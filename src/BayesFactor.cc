#include "HiggsAnalysis/CombinedLimit/interface/BayesFactor.h"
#include "HiggsAnalysis/CombinedLimit/interface/CascadeMinimizer.h"
#include "RooStats/ModelConfig.h"
#include "RooProjectedPdf.h"
#include "RooRealVar.h"
#include <iostream>
#include <cmath>

using namespace RooStats;
using namespace std;

bool  BayesFactor::marginalizeNuisances_ = false;
bool  BayesFactor::bgOnly_ = false;
float BayesFactor::muMax_ = 2.0;

BayesFactor::BayesFactor() :
    LimitAlgo("BayesFactor specific options")
{
  options_.add_options()
  ("marginalizeNuisances",  "Marginalize nuisances instead of using the cascade minimizer")
  ("bgOnly",  "Compute background likelihood only")
//  ("addDatacard",        boost::program_options::value<string>(&muMax_)->default_value(muMax_), "maximum value of signal strength mu that is scanned");
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
  /*
  RooAbsPdf * projectedpdf = sigpdf->createProjection ( *rrv );
  const RooArgSet * observables = w->set( "observables" );
  float ret = projectedpdf->getVal( *observables);
  float llhdedge = this->getLikelihood ( w, mc_s, data, max );
  float minllhd = 1e-10 * ret;
  if ( llhdedge < minllhd )
  {
    // the last bin doesnt contribute anymore to the mean likelihood!
    cout << "[BayesFactor::signalIntegralOverMu] llhd of last mu bin is very small:"
         << llhdedge <<"<"<< minllhd << endl;
    cout << "[BayesFactor::signalIntegralOverMu] consider a smaller value for mumax: " << muMax_ << endl;
  }
  float maxllhd = 0.01 * ret;
  if ( llhdedge > maxllhd )
  {
    // the last bin contributes still significantly to the
    // mean likelihood!
    cout << "[BayesFactor::signalIntegralOverMu] llhd of last mu bin is still large: "
         << llhdedge << ">" << maxllhd << "." << endl;
    cout << "[BayesFactor::signalIntegralOverMu] consider a higher value for mumax: " << muMax_ << endl;
  }
  */
  float sum = 0.;
  float mmaxllhd=-1.;
  // int maxllhdindex=-1;
  int nStepsMu_ = 20;
  float min = 0.;
  float lastllhd = 0.;
  for ( int i = 0; i < nStepsMu_; i++ )
  {
    float mu = min + (max-min) * float(i) / (nStepsMu_-1);
    lastllhd = this->getLikelihood ( w, mc_s, data, mu );
    // cout << "[BayesFactor] mu=" << mu << " llhd=" << sigllhd << endl;
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

  // cout << "[BayesFactor] maximum at " << maxllhdindex << endl;
  // cout << "[BayesFactor] sigllhd=" << sigllhd << endl;
  return pair < float, float > ( sigllhd, mmaxllhd );
}

/*
void BayesFactor::computeSignalStrength ( RooStats::ModelConfig *mc_s, RooAbsData &data ) const
{
    RooFitResult * res_s = mc_s->GetPdf()->fitTo(data,
    RooFit::Save(1),
    RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()),
    RooFit::Strategy(minimizerStrategy_),
    RooFit::Extended(mc_s->GetPdf()->canBeExtended()),
    constCmdArg_s, minosCmdArg);
       
};*/

bool BayesFactor::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
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
  cout << "[BayesFactor:result] bgllhd=" << bgllhd << ", signal=" << sigllhd 
       << ", maxsigllhd=" << maxllhd << ", lnK=" << limit 
       << ", sig=" << significance << endl;
  limitErr=0;
  hint=0;
  return true;
}
