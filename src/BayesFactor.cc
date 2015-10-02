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
float BayesFactor::muMax_ = 2.0;

BayesFactor::BayesFactor() :
    LimitAlgo("BayesFactor specific options")
{
  options_.add_options()
  ("marginalizeNuisances",  "Marginalize nuisances instead of using the cascade minimizer")
  ("muMax",        boost::program_options::value<float>(&muMax_)->default_value(muMax_), "maximum value of signal strength mu that is scanned");
}

BayesFactor::~BayesFactor(){
}

void BayesFactor::applyOptions(const boost::program_options::variables_map &vm) 
{
  marginalizeNuisances_ = vm.count("marginalizeNuisances"); 
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

float BayesFactor::signalIntegralOverMu ( RooWorkspace *w, RooStats::ModelConfig *mc_s, RooAbsData & data, float max ) const
{
  RooRealVar * rrv = w->var("r"); 
  const RooArgSet * nuisances = w->set( "nuisances" );
  rrv->setMin ( 0. );
  rrv->setMax ( max );
  RooAbsPdf * sigpdf = mc_s->GetPdf();
  RooAbsReal * nll = sigpdf->createNLL ( data );
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
  minim.setNuisanceParameters ( nuisances );
  minim.minimize( 0 );
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
    cout << "[BayesFactor::signalIntegralOverMu] consider a smaller value for " << muMax_ << endl;
  }
  float maxllhd = 0.01 * ret;
  if ( llhdedge > maxllhd )
  {
    // the last bin contributes still significantly to the
    // mean likelihood!
    cout << "[BayesFactor::signalIntegralOverMu] llhd of last mu bin is still large: "
         << llhdedge << ">" << maxllhd << "." << endl;
    cout << "[BayesFactor::signalIntegralOverMu] consider a higher value for " << muMax_ << endl;
  }
  return ret;

  /*
  float sum = 0.;
  float maxllhd=-1.;
  int maxllhdindex=-1;
  vector < float > llhds ( nStepsMu_ );
  for ( int i = 0; i < nStepsMu_; i++ )
  {
    float mu = min + (max-min) * float(i) / (nStepsMu_-1);
    float sigllhd = this->getLikelihood ( w, mc_s, data, mu );
    cout << "[BayesFactor] mu=" << mu << " llhd=" << sigllhd << endl;
    sum += sigllhd;
    llhds.push_back ( sigllhd );
    if ( sigllhd > maxllhd )
    {
      maxllhd=sigllhd;
      maxllhdindex=i;
    }
  }
  cout << "[BayesFactor] maximum at " << maxllhdindex << endl;
  float sigllhd = sum / nStepsMu_;
  cout << "[BayesFactor] projected = " << ret << " min=" << minrange << endl;
  */
}

bool BayesFactor::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
  float bgllhd = this->getLikelihood ( w, mc_b, data, -1 );
  float sigllhd = this->signalIntegralOverMu ( w, mc_s, data, muMax_ );
  limit=log ( sigllhd/bgllhd );
  cout << "[BayesFactor] bgllhd=" << bgllhd << ", signal=" << sigllhd << ", lnK=" << limit << endl;
  limitErr=0;
  hint=0;
  return true;
}
