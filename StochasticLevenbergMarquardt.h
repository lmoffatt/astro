#ifndef STOCHASTICLEVENBERGMARQUARDT_H
#define STOCHASTICLEVENBERGMARQUARDT_H
#include "Parameters.h"
#include <cmath>



class CortexLikelihood;


class StochasticLevenbergMarquardtDistribution: public BaseObject
{
public:

  StochasticLevenbergMarquardtDistribution& optimize();




  StochasticLevenbergMarquardtDistribution(const CortexLikelihood *f
                                           , double beta,
                                           const Parameters& initialParam
                                           , std::size_t numSamples
                                           , std::size_t n_skip_iter
                                           , double maxDuration_mins
                                           , const std::string& name
                                           ,bool storeParameters
                                           ,bool storeData);



  ~StochasticLevenbergMarquardtDistribution(){}


  friend std::ostream& operator<<(std::ostream& s, StochasticLevenbergMarquardtDistribution& LM);



private:

  struct state
  {

    state(std::size_t nData_, std::size_t nPar_);
    bool update(const CortexLikelihood* CL_
                , double beta,
                Parameters &&param
                , double dx
                , double initLanda_
                ,double multLanda
                ,double maxLanda
                ,double tol);


    bool update_landa(const CortexLikelihood* CL_,double landa_,double tol);

    std::vector<double> w_;
    double logLik;
    double logPrior;
    double beta_;
    double beta_LogLik_;
    double beta_LogLik_Next_exp_;

    double landa_;
    Parameters Param_;
    Parameters ParamNext_;
    std::vector<std::vector<double>> P_exp_;
    std::vector<std::vector<double>> P_exp_Next_;
    double logLik_Next_;
    double logPrior_Next_;
    double beta_LogLik_Next_;

    double epsilon_beta_log_Lik_;


    bool isJacobianInvalid_;
    std::vector< std::vector< double> > J_;
    std::vector<double> prior_G_;
    std::vector<double> beta_G_;
    std::vector<double> epsilon_;
    std::vector< std::vector<double> > beta_JTWJ_;
    std::vector< std::vector<double> > beta_JTWJ_landa_;
    std::vector< std::vector<double> > beta_JTWJinv_;
    std::vector<double> beta_d_;


  public:
  };

  struct MetropolisHastings_test
  {
    MetropolisHastings_test()=default;

    bool operator()(std::mt19937_64 &mt
                    ,const state& current
                    ,const state &candidate);


    double logPcurrent;
    double logPcandidate;

    double logQforward;

    double logQbackward;

    double logA;

    double A;
    double r;
    bool accept_;

  };




  struct step
  {
    step(std::size_t nData, std::size_t nPar)
      :
        s_Curr_(nData,nPar),s_Test_(nData,nPar),test_(),nTryNext_(0),nIter_(0)
      ,timeOpt_(0),timeIter_(0),fail_(false)
    {}
    state s_Curr_;
    state s_Test_;
    MetropolisHastings_test test_;

    std::size_t nTryNext_;
    std::size_t nIter_;
    double timeOpt_;
    double timeIter_;
    bool fail_;

  };

  struct stochastic_iter
  {
    void operator()(const CortexLikelihood* CL_,double beta,std::mt19937_64 &mt,
                    step& s)const;

    stochastic_iter(double tol=1.0,
                    double dx=1e-6,
                    double initLanda=1,
                    double multLanda=2,
                    double maxLanda=1e5,std::size_t maxTries=10): maxTries_(maxTries),dx_(dx),initLanda_(initLanda),multLanda_(multLanda),maxLanda_(maxLanda),tol_(tol){}

    std::size_t maxTries_;
    double dx_;
    double initLanda_;
    double multLanda_;
    double maxLanda_;
    double tol_;

  };
  struct step_report
  {
    std::string title()const ;
    std::string operator()(const step& s)const;
  };

  struct state_store;
  struct mcmc
  {
    mcmc(std::size_t nsamples_,const state_store& ss);
    bool isFull_;
    std::size_t nsamples_;
    std::vector<double> logLiks;
    std::vector<Parameters> parameters;
    std::vector<std::vector<std::vector<double>>> data;
    double sumLogLik;
    double sumSqrLogLik;
    double varLogLik(){return sumSqrLogLik/nsamples_-meanLogLik()*meanLogLik();}
    double meanLogLik(){return sumLogLik/nsamples_;}

    bool isFull()const {return nsamples_>=(logLiks.size()-1);}
  };

  struct state_store
  {
    std::size_t num_skip_samples_;
    bool storeParameters;
    bool storeData;


    state_store(std::size_t num_skipped_samples,
                bool storeParam=false, bool storeDat=false)
      : num_skip_samples_(num_skipped_samples),
        storeParameters(storeParam),storeData(storeDat){}


    bool operator()(mcmc& m, const step& s)
    {
      std::size_t isample=(s.nIter_)/(num_skip_samples_+1);
      if ((isample*(num_skip_samples_+1)==s.nIter_))
        {
        this->operator()(m,s.s_Curr_,isample);
          return true;
        }
      else return false;
    }

    void operator()(mcmc& m, const state& s, std::size_t isample)
    {
      if (!m.isFull())
        {
          m.nsamples_++;
          m.logLiks[isample]=s.logLik;
          m.sumLogLik+=s.logLik;
          m.sumSqrLogLik+=s.logLik*s.logLik;
          if (storeParameters)
            m.parameters[isample]=s.Param_;
          if (storeData)
            m.data[isample]=s.P_exp_;
        }
    }
  };

  struct Finalization_report

  {
    Finalization_report(): surpassDuration_(false),full_mcmc_(false),failStep_(false){}
    bool surpassDuration_;
    bool full_mcmc_;
    bool failStep_;

  };

  struct
      Finalization_test
  {
    Finalization_test(double maxDur_in_min);

    bool operator()(const mcmc& m,const step& s, Finalization_report& r)const;
    double maxDur_in_min_;


  public:
   };



  void initialize();
  bool iterate();
  void report_step();
  void report_Optimization();



private:
  std::chrono::steady_clock::time_point startTime_;

  std::string fname_;
  std::ofstream os_;
  const CortexLikelihood* CL_;
  Parameters ParamInitial_;

  std::size_t nPar_;
  std::size_t nData_;

  step step_;
  stochastic_iter stoc_iter_;
  step_report step_report_;
  state_store storeIt_;
  mcmc mcmc_;
  Finalization_test f_test_;
  Finalization_report f_report_;












  // BaseClass interface
public:
  static std::string ClassName(){return "LevenbergMarquardtDistribution";}
  virtual std::string myClass() const override{return ClassName();}

  // BaseObject interface
public:
  virtual StochasticLevenbergMarquardtDistribution *create() const override{ return nullptr;}
  virtual std::ostream &writeBody(std::ostream &s) const override{ return s;}
  virtual void clear() override{}
  virtual bool readBody(std::string &, std::istream &) override{ return false;}

protected:
  virtual void update() override{}
};

#endif // STOCHASTICLEVENBERGMARQUARDT_H
