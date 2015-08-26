#ifndef LEVENBERGMARQUARDT
#define LEVENBERGMARQUARDT
#include "Parameters.h"
#include <cmath>

class ABC_Freq_obs
{
public:
  virtual const std::vector<std::vector<double>>& n_obs()const=0;

  virtual const std::vector<double>& n_obs(unsigned i)const=0;
  virtual const std::vector<double>& ntot_obs()const=0;
  unsigned numSamples()const;

  unsigned numCells()const;

  std::size_t numDF()const
  {
    std::size_t sum=0;
    for (std::size_t i =0; i<ntot_obs().size(); ++i)
      {
        if (ntot_obs()[i]>0)
          sum+=n_obs(i).size()-1;
      }
    return sum;
  }

  virtual ~ABC_Freq_obs(){}
  virtual std::ostream& put(std::ostream& s)const=0;
};


class ABC_Multinomial_Model
{
public:

  virtual void setPrior(const Parameters& parameters)=0;
  virtual const ABC_Freq_obs& getData()const=0;

  virtual const Parameters& getPrior()const=0;
  virtual  std::vector<std::vector<double>> p_exp(const Parameters& parameters)const=0;

  virtual  std::vector<std::vector<double>> J(const Parameters& p,
                                              const std::vector<std::vector<double>>& p_ex ,
                                              double delta=1e-4)const;

  std::vector<std::vector<double>>
  logLikCells(const std::vector<std::vector<double> > &p) const;

  std::vector<double>
  logLikSamples(const std::vector<std::vector<double> > &p) const;

  double logLik(const std::vector<std::vector<double>>& p)const;

  double logLik(const Parameters& p)const;

  std::vector<double> epsilon(const std::vector<std::vector<double>>& p)const;

  const std::vector<double> weight(const std::vector<std::vector<double>>& p)const;

  double logPrior(const Parameters& p) const;

  std::vector<double> PriorGradient(const Parameters& p) const;




  virtual ~ABC_Multinomial_Model();
};



class LevenbergMarquardtMultinomial
{
public:

  Parameters OptimParameters()const;

  std::size_t numEval()const;
  std::size_t numIter()const;
  double LogLik()const;
  double PostLogLik()const;
  std::vector<double> Gradient()const;


  LevenbergMarquardtMultinomial& optimize();


  LevenbergMarquardtMultinomial(ABC_Multinomial_Model* f,
                                const Parameters& initialParam,
                                std::size_t numIterations);

  double getEvidence()const;

  double getLogPostLik()const;

  double logDetPriorCov()const;
  double logDetPostCov()const;
  double logDetPostStd()const;


  LevenbergMarquardtMultinomial(const LevenbergMarquardtMultinomial& other);

  friend void swap(LevenbergMarquardtMultinomial& one, LevenbergMarquardtMultinomial& other);

  LevenbergMarquardtMultinomial& operator=(const LevenbergMarquardtMultinomial& other);

  LevenbergMarquardtMultinomial();

  ~LevenbergMarquardtMultinomial(){}
  std::string report()const;

  // void reset(const SimParameters& sp,const Treatment& tr);


  friend std::ostream& operator<<(std::ostream& s, LevenbergMarquardtMultinomial& LM);

private:
  ABC_Multinomial_Model* f_;

  std::vector<double> w_;

  Parameters ParamInitial_;

  std::size_t nPar_;
  std::size_t nData_;

 // parameters of the optimization
  /// delta x used for Jacobian approximation
  double dx_;
  std::size_t maxIter_;
  std::size_t maxFeval_;

  double ParamChangeMin_;
  double PostLogLikChangeMin_;
  double GradientNormPostMin_;

  double maxLanda_;

  // variables that change on each iteration

  double landa_;
  double landa0_;
  double v_;

  std::size_t nIter_;
  std::size_t nFeval_;
  std::size_t nDF_;

  // logLiks
  double logLikCurr_;
  double logLikNew_;
  double logLikNew0_;


//logPriors
  double logPriorCurr_;
  double logPriorNew_;
  double logPriorNew0_;



  double logPostLikCurr_;
  double logPostLogLikNew_;
  double logPostLikeNew0_;

  Parameters ParamCurr_;
  Parameters ParamNew_;
  Parameters ParamNew0_;

  std::vector<std::vector<double>> P_expCurr_;
  std::vector<std::vector<double>> P_exp_New_;
  std::vector<std::vector<double>> P_exp_New0_;

  std::vector< std::vector< double> > J_;
  std::vector<double> prior_G_;

  std::vector<double> epsilon_;


  std::vector<double> G_;
  std::vector< std::vector<double> > JTWJ_;
  std::vector< std::vector<double> > JTWJ_landa_;
  std::vector< std::vector<double> > JTWJinv_;

  // postLiks

  std::vector<double> d_;


  bool surpassIter_;
  bool surpassFeval_;
  bool surpassLanda_;

  double ParamChange_;
  double PostlogLikChange_;
  double GradNormPost_;

  bool smallParamChange_;

  bool smallSSChange_;

  bool smallGradient_;
  bool isNanLogPostLik_;

  bool meetConvergenceCriteria();

  void initialize();
  void iterate();

  void computeJacobian();
  void computeSearchDirection();
  void updateLanda();



};


#endif // LEVENBERGMARQUARDT

