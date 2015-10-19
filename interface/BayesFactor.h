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
#include "RooRealVar.h"
#include <utility>
#include <vector>
#include <string>

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
  /// maximize the signal strength, taken from MaxLikelihoodfit
  RooRealVar * computeSignalStrength ( RooWorkspace *, RooStats::ModelConfig *mc_s, 
                               RooAbsData & ) const;
  double getNLL ( RooWorkspace *, const RooStats::ModelConfig *, 
                 RooAbsData &, double r ) const ;
  /// integrate out signal strength mu,
  /// report integral (.first) *and* maximum (.second)
  std::pair < long double, long double > signalIntegralOverMu ( RooWorkspace *w, 
      RooStats::ModelConfig *mc_s, RooAbsData &, double max ) const;

  //true: nuisances are marginalized, false: nuisances are treated w/ CascadeMinimizer
  static bool  marginalizeNuisances_;
  static bool  bgOnly_; // only compute likelihood of background
  static double muMax_; //maximum signal strength that we scan
  static int   minimizerStrategy_;
  static double preFitValue_; // value of mu before the fit
  // static std::vector < std::string > addDataCards_; // add datacards with more signals
};


#endif
