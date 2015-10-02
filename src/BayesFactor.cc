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
  //float mmaxllhd=-1.;
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
    /*
    if ( lastllhd > maxllhd )
    {
      //mmaxllhd=sigllhd;
      maxllhdindex=i;
    }*/
  }
  float sigllhd = sum / nStepsMu_;
  float maxl=0.01 * sigllhd;

  if ( lastllhd > maxl )
  {
    cout << "[BayesFactor] the llhd of last mu bin is very large: " << lastllhd << ">"
         << maxl << ". Consider choosing a larger mumax value." << endl;
  }
  float minl=1e-10 * sigllhd;
  if ( lastllhd < minl )
  {
    cout << "[BayesFactor] the llhd of the last mu bin is very small: " << lastllhd 
         << "<" << minl << ". Consider choosing a smaller mumax value." << endl;
  }

 // cout << "[BayesFactor] maximum at " << maxllhdindex << endl;
  cout << "[BayesFactor] sigllhd=" << sigllhd << endl;
  return sigllhd;
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
