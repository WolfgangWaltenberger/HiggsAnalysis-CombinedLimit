#ifndef HiggsAnalysis_CombinedLimit_BayesFactor_h
#define HiggsAnalysis_CombinedLimit_BayesFactor_h
/** \class BayesFactor
 *
 * Compute Bayes factor for a given set of signals, versus 
 * background theory only
 *
 * \author Wolfgang Waltenberger
 *
 */

#include "HiggsAnalysis/CombinedLimit/interface/LimitAlgo.h"
//#include <RooArgList.h>
//#include <RooFitResult.h>
//#include <boost/utility.hpp>

class BayesFactor : public LimitAlgo {
public:
  BayesFactor() ;
  virtual const std::string & name() const {
    static const std::string name("BayesFactor");
    return name;
  }
  ~BayesFactor();
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;

protected:
  virtual bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);

private:
  float getLikelihood ( RooWorkspace *, const RooStats::ModelConfig *, 
                 RooAbsData &, float r ) const ;

  //true: nuisances are marginalized, false: nuisances are treated w/ CascadeMinimizer
  static bool marginalizeNuisances_;
  //maximum signal strength
  static float mumax_;
};


#endif
