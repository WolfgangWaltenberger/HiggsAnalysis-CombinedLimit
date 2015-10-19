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
double BayesFactor::muMax_ = 2.0;
int   BayesFactor::minimizerStrategy_  = 1;
double BayesFactor::preFitValue_ = 1.0;
// std::vector < std::string > BayesFactor::addDataCards_;

BayesFactor::BayesFactor() :
    LimitAlgo("BayesFactor specific options")
{
  options_.add_options()
  ("marginalizeNuisances",  "Marginalize nuisances instead of using the cascade minimizer")
  ("bgOnly",  "Compute background likelihood only")
//  ("addDataCards",        boost::program_options::value< vector < string > >(&addDataCards_)->multitoken(), "additional datacards, should only vary in signal hypothesis" )
  ("minimizerStrategy",  boost::program_options::value<int>(&minimizerStrategy_)->default_value(minimizerStrategy_),      "Stragegy for minimizer")
  ("preFitValue",        boost::program_options::value<double>(&preFitValue_)->default_value(preFitValue_),  "Value of signal strength pre-fit")
  ("muMax",        boost::program_options::value<double>(&muMax_)->default_value(muMax_), "maximum value of signal strength mu that is scanned");
}

BayesFactor::~BayesFactor(){
}

void BayesFactor::applyOptions(const boost::program_options::variables_map &vm)
{
  marginalizeNuisances_ = vm.count("marginalizeNuisances");
  bgOnly_ = vm.count("bgOnly");
}

double BayesFactor::getNLL ( RooWorkspace * w,
    const RooStats::ModelConfig * mc, RooAbsData & data, double r ) const
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
    double ret = projectedpdf->getVal( *observables);
    //cout << "[BayesFactor] projectedpdf= " << ret << endl;
    delete projectedpdf;
    return ret;
  }
  RooAbsReal * nll = bgpdf->createNLL ( data );
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
  minim.setNuisanceParameters ( nuisances );
  minim.minimize( 0 );
  double ret = nll->getVal();
  //double ret = exp ( -nll->getVal() );
  delete nll;
  return ret;
}

pair < long double, long double > BayesFactor::signalIntegralOverMu ( RooWorkspace *w, RooStats::ModelConfig *mc_s, RooAbsData & data, double max ) const
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
  long double sum = 0.;
  long double mmaxllhd=-1.;
  // double maxnll = -9e99;
  int nStepsMu_ = 20;
  long double min = 0.;
  long double lastllhd = 0.;
  for ( int i = 0; i < nStepsMu_; i++ )
  {
    long double mu = min + (max-min) * double(i) / (nStepsMu_-1);
    long double nll = this->getNLL ( w, mc_s, data, mu );
    /*
    if ( nll > maxnll )
    {
      maxnll=nll;
    } */
    lastllhd = expl ( - nll );
    sum += lastllhd;
    if ( lastllhd > mmaxllhd )
    {
      mmaxllhd=lastllhd;
    }
  }
  long double sigllhd = sum / nStepsMu_;
  long double maxl=0.01 * sigllhd;

  if ( lastllhd > maxl )
  {
    cout << "[BayesFactor:error] the llhd of last mu bin is very large: "
         << lastllhd << ">" << maxl << ". Consider choosing a larger mumax value."
         << endl;
  }
  long double minl=1e-10 * sigllhd;
  if ( lastllhd < minl )
  {
    cout << "[BayesFactor:error] the llhd of the last mu bin is very small: "
         << lastllhd << "<" << minl << ". Consider choosing a smaller mumax value."
         << endl;
  }

  return pair < long double, long double > ( sigllhd, mmaxllhd );
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
  double bgNLL = this->getNLL ( w, mc_b, data, -1 );
  if (bgOnly_)
  {
    cout << "[BayesFactor:result] bgNLL=" << bgNLL << endl;
    limit=bgNLL;
    return true;
  }
  pair < long double, long double> sig = this->signalIntegralOverMu ( w, mc_s, data, muMax_ );
  long double sigllhd = sig.first;
  long double maxllhd = sig.second;
  // long double maxllhd = exp ( - maxnll );
  limit=log ( sigllhd ) + bgNLL;
  double significance = sqrt ( 2*(log ( maxllhd ) + bgNLL ) );
  RooRealVar * mu = this->computeSignalStrength ( w, mc_s, data );
  cout << "[BayesFactor:result] bgNLL=" << bgNLL << ", signal=" << sigllhd
       << ", maxsigllhd=" << maxllhd << ", lnK=" << limit
       << ", sig=" << significance << ", mu=" << mu->getVal()
       << ", muerr=" << mu->getError()
       << endl;
  limitErr=0;
  hint=0;
  return true;
}
