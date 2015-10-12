#include "HiggsAnalysis/CombinedLimit/interface/BayesFactor.h"
#include "HiggsAnalysis/CombinedLimit/interface/CascadeMinimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h" 
#include "HiggsAnalysis/CombinedLimit/interface/Combine.h"  
#include "RooStats/ModelConfig.h"
#include "RooProjectedPdf.h"
#include "RooFitResult.h"
#include <iostream>
#include <cmath>

using namespace RooStats;
using namespace std;

bool  BayesFactor::marginalizeNuisances_ = false;
bool  BayesFactor::bgOnly_ = false;
float BayesFactor::muMax_ = 2.0;
int   BayesFactor::minimizerStrategy_  = 1;
float BayesFactor::preFitValue_ = 1.0;
// std::vector < std::string > BayesFactor::addDataCards_;

BayesFactor::BayesFactor() :
    LimitAlgo("BayesFactor specific options")
{
  options_.add_options()
  ("marginalizeNuisances",  "Marginalize nuisances instead of using the cascade minimizer")
  ("bgOnly",  "Compute background likelihood only")
//  ("addDataCards",        boost::program_options::value< vector < string > >(&addDataCards_)->multitoken(), "additional datacards, should only vary in signal hypothesis" )
  ("minimizerStrategy",  boost::program_options::value<int>(&minimizerStrategy_)->default_value(minimizerStrategy_),      "Stragegy for minimizer")
  ("preFitValue",        boost::program_options::value<float>(&preFitValue_)->default_value(preFitValue_),  "Value of signal strength pre-fit")
  ("muMax",        boost::program_options::value<float>(&muMax_)->default_value(muMax_), "maximum value of signal strength mu that is scanned");
}

BayesFactor::~BayesFactor(){
}

void BayesFactor::applyOptions(const boost::program_options::variables_map &vm)
{
  marginalizeNuisances_ = vm.count("marginalizeNuisances");
  bgOnly_ = vm.count("bgOnly");
}

float BayesFactor::getLikelihood ( RooWorkspace * w,
    const RooStats::ModelConfig * mc, RooAbsData & data, float r ) const
{
  const RooArgSet * nuisances = w->set( "nuisances" );
  const RooArgSet * observables = w->set( "observables" );
//  w->var("overall3")->setMin(-9.);
//  w->var("overall3")->setMax(9.);
  if ( r >= 0. )
  {
    RooRealVar * rrv = w->var("r");
    rrv->setVal ( r );
    rrv->setConstant ( true );
  }
  RooAbsPdf * bgpdf = mc->GetPdf();
  if ( marginalizeNuisances_ )
  {
    // cout << "[BayesFactor] marginalizeNuisances!" << endl;
    RooAbsPdf * projectedpdf = bgpdf->createProjection ( *nuisances );
    float ret = projectedpdf->getVal( *observables);
    //cout << "[BayesFactor] projectedpdf= " << ret << endl;
    delete projectedpdf;
    return ret;
  }
  RooAbsReal * nll = bgpdf->createNLL ( data );
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
  minim.setNuisanceParameters ( nuisances );
  minim.minimize( 0 );
  float ret = exp ( -nll->getVal() );
  delete nll;
  return ret;
}

pair < float, float > BayesFactor::signalIntegralOverMu ( RooWorkspace *w, RooStats::ModelConfig *mc_s, RooAbsData & data, float max ) const
{
  RooRealVar * rrv = w->var("r");
  const RooArgSet * nuisances = w->set( "nuisances" );
  rrv->setMin ( 0. );
  rrv->setMax ( max );
  RooAbsPdf * sigpdf = mc_s->GetPdf();
  RooAbsReal * nll = sigpdf->createNLL ( data, RooFit::Constrain ( *nuisances ) );
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained );
  minim.setNuisanceParameters ( nuisances );
  minim.minimize( 0 );
  float sum = 0.;
  float mmaxllhd=-1.;
  int nStepsMu_ = 20;
  float min = 0.;
  float lastllhd = 0.;
  for ( int i = 0; i < nStepsMu_; i++ )
  {
    float mu = min + (max-min) * float(i) / (nStepsMu_-1);
    lastllhd = this->getLikelihood ( w, mc_s, data, mu );
    sum += lastllhd;
    if ( lastllhd > mmaxllhd )
    {
      mmaxllhd=lastllhd;
    }
  }
  float sigllhd = sum / nStepsMu_;
  float maxl=0.01 * sigllhd;

  if ( lastllhd > maxl )
  {
    cout << "[BayesFactor:error] the llhd of last mu bin is very large: "
         << lastllhd << ">" << maxl << ". Consider choosing a larger mumax value."
         << endl;
  }
  float minl=1e-10 * sigllhd;
  if ( lastllhd < minl )
  {
    cout << "[BayesFactor:error] the llhd of the last mu bin is very small: "
         << lastllhd << "<" << minl << ". Consider choosing a smaller mumax value."
         << endl;
  }

  return pair < float, float > ( sigllhd, mmaxllhd );
}

RooRealVar * BayesFactor::computeSignalStrength ( RooWorkspace *w, 
    RooStats::ModelConfig *mc_s, RooAbsData &data ) const
{
  // code stolen from MaxLikelihoodFit, maximize the signal strength
  RooRealVar * r = w->var("r");
  r->setVal( preFitValue_ );
  r->setConstant(false);
  const RooCmdArg & constCmdArg_s = RooFit::Constrain(*mc_s->GetNuisanceParameters());
  const RooCmdArg &minosCmdArg = RooFit::Minos(*mc_s->GetParametersOfInterest());
  RooFitResult * res_s = 0;
  {
     CloseCoutSentry sentry(verbose < 5 );
     res_s = mc_s->GetPdf()->fitTo(data,
          RooFit::Save(1),
          RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()),
          RooFit::Strategy(minimizerStrategy_),
          RooFit::Extended(mc_s->GetPdf()->canBeExtended()), 
          constCmdArg_s, minosCmdArg);
  }
  if (!res_s)
  {
    cout << "[BayesFactor:error] res_s not defined." << endl;
    exit(-1);
  }
  RooRealVar *rf = dynamic_cast<RooRealVar*>(res_s->floatParsFinal().find(r->GetName()));
  if (!rf)
  {
    cout << "[BayesFactor:error] rf not defined." << endl;
    exit(-1);
  }
  if ( verbose > 5 ) res_s->Print("V");
  return rf;
};

bool BayesFactor::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint)
{
  float bgllhd = this->getLikelihood ( w, mc_b, data, -1 );
  if (bgOnly_)
  {
    cout << "[BayesFactor:result] bgllhd=" << bgllhd << endl;
    limit=bgllhd;
    return true;
  }
  pair <float,float> sig = this->signalIntegralOverMu ( w, mc_s, data, muMax_ );
  float sigllhd = sig.first;
  float maxllhd = sig.second;
  limit=log ( sigllhd/bgllhd );
  float significance = sqrt ( 2*log ( maxllhd / bgllhd ) );
  RooRealVar * mu = this->computeSignalStrength ( w, mc_s, data );
  cout << "[BayesFactor:result] bgllhd=" << bgllhd << ", signal=" << sigllhd
       << ", maxsigllhd=" << maxllhd << ", lnK=" << limit
       << ", sig=" << significance << ", mu=" << mu->getVal()
       << ", muerr=" << mu->getError()
       << endl;
  limitErr=0;
  hint=0;
  return true;
}
