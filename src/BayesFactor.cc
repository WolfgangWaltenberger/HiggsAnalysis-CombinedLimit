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
float BayesFactor::mumax_ = 2.0;

BayesFactor::BayesFactor() :
    LimitAlgo("BayesFactor specific options")
{
  options_.add_options()
  ("marginalizeNuisances",  "Marginalize nuisances instead of using the cascade minimizer");
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

bool BayesFactor::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
  float bgllhd = this->getLikelihood ( w, mc_b, data, -1 );
  float sigllhd1 = this->getLikelihood ( w, mc_s, data, 1 );
  float sigllhd05 = this->getLikelihood ( w, mc_s, data, 0.5 );
  cout << "[BayesFactor] bgllhd=" << bgllhd << ", siglllhd(1)=" << sigllhd1 << ", sigllhd05=" << sigllhd05 << endl;
  limit=log ( sigllhd1/bgllhd );
  limitErr=0;
  hint=0;
  return true;
}
