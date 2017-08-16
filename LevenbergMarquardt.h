#ifndef LEVENBERGMARQUARDT
#define LEVENBERGMARQUARDT
#include "Parameters.h"
#include <cmath>



class CortexLikelihood;



class ABC_Freq_obs
{
public:

  virtual const std::vector<std::vector<double>>& n_obs()const=0;

  virtual const std::vector<double>& bin_dens()const=0;

  virtual const std::vector<double>& n_obs(unsigned i)const=0;
  virtual const std::vector<double>& ntot_obs()const=0;
  unsigned numSamples()const;

  unsigned numCells()const;

  std::size_t numDF()const;

  virtual ~ABC_Freq_obs(){}
  virtual std::ostream& put(std::ostream& s)const=0;
};





class ABC_Distribution_Model
{
public:

  virtual void setPrior(const Parameters& parameters)=0;
  virtual const ABC_Freq_obs& getData()const=0;

  virtual const Parameters& getPrior()const=0;


  virtual  std::vector<std::vector<double>> f(const Parameters& parameters, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const=0;

  virtual  std::vector<std::vector<double>> f(const std::vector<double>& p)const;

  virtual  std::vector<std::vector<double>> f(const std::vector<double>& p, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const;

  std::vector<std::vector<double>>
  J(const Parameters& p,
    const std::vector<std::vector<double>>& f0 ,
    double delta=1e-4)const;
  virtual std::vector<std::vector<double>>
  logLikCells(const std::vector<std::vector<double> > &p) const=0;

  virtual std::vector<double>
  logLikSamples(const std::vector<std::vector<double> > &p) const=0;

  virtual double logLik(const std::vector<std::vector<double>>& p)const=0;

  virtual double logLik(const Parameters& p)const;

  virtual double logLik(const std::vector<double>& o)const;



  virtual std::vector<double> epsilon(const std::vector<std::vector<double>>& p)const=0;

  virtual const std::vector<double> weight(const std::vector<std::vector<double>>& p)const=0;

  virtual double logPrior(const Parameters& p) const;


  virtual std::vector<double> PriorGradient(const Parameters& p) const;


  virtual ~ABC_Distribution_Model(){}
  double logPrior(const std::vector<double> &o) const;
};



class ABC_Multinomial_Model: virtual public ABC_Distribution_Model
{
public:

  virtual void setPrior(const Parameters& parameters)=0;
  virtual const ABC_Freq_obs& getData()const=0;

  virtual const Parameters& getPrior()const=0;


  virtual  std::vector<std::vector<double>> f(const Parameters& parameters, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const=0;

  virtual

  std::vector<std::vector<double>>
  logLikCells(const std::vector<std::vector<double> > &p) const;

  std::vector<double>
  logLikSamples(const std::vector<std::vector<double> > &p) const;

  double logLik(const std::vector<std::vector<double>>& p)const;


  std::vector<double> epsilon(const std::vector<std::vector<double>>& p)const;

  const std::vector<double> weight(const std::vector<std::vector<double>>& p)const;


  virtual ~ABC_Multinomial_Model(){}
};





class ABC_MultiPoison_Model:virtual public ABC_Distribution_Model

{
public:

//  virtual void setPrior(const Parameters& parameters)=0;
//  virtual const ABC_Freq_obs& getData()const=0;

//  virtual const Parameters& getPrior()const=0;


  //virtual  std::vector<std::vector<double>> f(const Parameters& parameters)const=0;


  std::vector<std::vector<double>>
  logLikCells(const std::vector<std::vector<double> > &landa) const override;

  std::vector<double>
  logLikSamples(const std::vector<std::vector<double> > &landa) const override;

  double logLik(const std::vector<std::vector<double>>& landa)const override;

  double logLik(const Parameters& p)const override;

  std::vector<double> epsilon(const std::vector<std::vector<double>>& landa)const override;

  const std::vector<double> weight(const std::vector<std::vector<double>>& landa)const override;



  virtual ~ABC_MultiPoison_Model(){}
};



class LevenbergMarquardtDistribution: public BaseObject
{
public:

  Parameters OptimParameters()const;

  std::size_t numEval()const;
  std::size_t numIter()const;
  double LogLik()const;
  double PostLogLik()const;
  std::vector<double> Gradient()const;


  LevenbergMarquardtDistribution& optimize();


  LevenbergMarquardtDistribution& optimize(std::string optname,
                                           double factor,
                                           std::size_t numSeeds,
                                           std::mt19937_64::result_type initseed=0);


  LevenbergMarquardtDistribution(const CortexLikelihood* f,
                                 const Parameters& initialParam,
                                 std::size_t numIterations,
                                 double maxDuration_mins,
                                 const std::string&  name);

  double getEvidence()const;

  double getLogPostLik()const;

  double logDetPriorCov()const;
  double logDetPostCov()const;
  double logDetPostStd()const;


  LevenbergMarquardtDistribution(const LevenbergMarquardtDistribution& other);

  friend void swap(LevenbergMarquardtDistribution& one, LevenbergMarquardtDistribution& other);

  LevenbergMarquardtDistribution& operator=(const LevenbergMarquardtDistribution& other);

  LevenbergMarquardtDistribution();

  ~LevenbergMarquardtDistribution(){}
  std::string report()const;

  // void reset(const SimParameters& sp,const Treatment& tr);


  friend std::ostream& operator<<(std::ostream& s, LevenbergMarquardtDistribution& LM);

private:
  std::chrono::steady_clock::time_point startTime_;

  std::string fname_;
  std::ofstream os_;
  const CortexLikelihood* CL_;
  std::vector<double> w_;
  Parameters ParamInitial_;

  std::size_t nPar_;
  std::size_t nData_;

  // parameters of the optimization
  /// delta x used for Jacobian approximation
  double dx_;
  double maxDur_in_min_;
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

  double timeOpt_;
  double timeIter_;
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

  bool surpassDuration_;
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

  bool isJacobianInvalid_;

  bool meetConvergenceCriteria();

  void initialize();
  void iterate();

  void computeJacobian();
  void computeSearchDirection();
  void updateLanda();




  // BaseClass interface
public:
  static std::string ClassName(){return "LevenbergMarquardtDistribution";}
  virtual std::string myClass() const override{return ClassName();}

  // BaseObject interface
public:
  virtual LevenbergMarquardtDistribution *create() const override{ return nullptr;}
  virtual std::ostream &writeBody(std::ostream &s) const override{ return s;}
  virtual void clear() override{}
  virtual bool readBody(std::string &, std::istream &, std::ostream& ) override{ return false;}

protected:
  virtual void update() override{}
};










#endif // LEVENBERGMARQUARDT

