#ifndef BAYESITERATION
#define BAYESITERATION


#include <vector>
#include <string>
#include <map>
#include "LevenbergMarquardt.h"

#include <random>


inline double mean(const std::vector<double>& x)
{
  double m=0;
  for (std::size_t i=0; i<x.size(); ++i)
    m+=x[i];
  return m/x.size();
}


inline std::vector<double> removeMean(const std::vector<double>& x,double m)
{
  std::vector<double> o(x.size());
  for (std::size_t i=0; i<x.size(); ++i)
    o[i]=x[i]-m;
  return o;
}



inline double AutoCorrelationMean(const std::vector<double>& d, std::size_t k)
{
  double sum=0;
  for (std::size_t i=0; i+k<d.size(); ++i)
    sum+=d[i]*d[i+k];
  sum/=d.size()-k;
  return sum;
}

inline double var(const std::vector<double>& d)
{
  double sum=0;
  for (std::size_t i=0; i<d.size(); ++i)
    sum+=d[i]*d[i];
  sum/=d.size();
  return sum;
}

inline std::vector<double> reduce(const std::vector<double>& x)
{
  std::vector<double> o(x.size()/2);
  for (std::size_t i=0; i<o.size(); ++i)
    {
      o[i]=x[2*i]+x[2*i+1];
    }
  return o;
}


inline double integratedAutoCorrelationTimed(const std::vector<double>& d,std::size_t n)
{
  double v=var(d);
  double tau=v;
  for (std::size_t i=1; i<n; ++i)
    {
      tau+=AutoCorrelationMean(d,i)*2.0;
    }
  tau/=v;
  if ((tau>n)&&(n*n<d.size()))
    {
      auto r=reduce(d);
      return integratedAutoCorrelationTimed(r,n);
    }
  else
    return tau;
}


inline double integratedAutoCorrelationTime(const std::vector<double>& x,std::size_t n)
{
  double m=mean(x);
  std::vector<double> d=removeMean(x,m);
  return integratedAutoCorrelationTimed(d,n);

}

class mcmc;
class BaseProposalDistribution;

class BaseMcmcState: public BaseObject
{

public:
  virtual  std::ostream &writeValues(std::ostream &s, std::size_t isample)const=0;
  virtual  std::ostream &writeMeans(std::ostream &s)const=0;
  virtual  std::ostream &writeValuesTitles(std::ostream &s)const=0;
  virtual   std::ostream &writeMeansTitles(std::ostream &s)const=0;

  virtual  std::ostream &writeYValuesTitles(std::ostream &s, const CortexLikelihood *CL)const=0;
  virtual  std::ostream &writeYValues(std::ostream &s, std::size_t isample, const CortexLikelihood *CL)const=0;


  virtual ~BaseMcmcState(){}
  static  BaseMcmcState* createChild(const std::string& className);
private:
  static std::map<std::string, BaseMcmcState*> childs_;
  static std::map<std::string, BaseMcmcState*> getChilds();



  // BaseObject interface
public:
  virtual BaseMcmcState *create() const =0;
};

class BaseBetaManagement: public BaseObject
{
public:
  virtual double beta(const mcmc* myMCMC)const=0;

  virtual ~BaseBetaManagement(){}
  static  BaseBetaManagement* createChild(const std::string& className);
private:
  static std::map<std::string, BaseBetaManagement*> childs_;
  static std::map<std::string, BaseBetaManagement*> getChilds();



  // BaseObject interface
public:
  virtual BaseBetaManagement *create() const =0;
};


class Beta_Ramp: public BaseBetaManagement
{


  // BaseClass interface
public:
  static std::string ClassName(){return "Beta_Ramp";}

  virtual std::__cxx11::string myClass() const override
  {
    return ClassName();
  }

  // BaseObject interface
public:
  Beta_Ramp(){}
  virtual Beta_Ramp *create() const override
  {
    return new Beta_Ramp;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
  }
  virtual void clear() override
  {
  }
  virtual bool readBody(std::__cxx11::string &line, std::istream &s) override
  {
  }

protected:
  virtual void update() override
  {
  }

  // BaseBetaManagement interface
public:
  virtual double beta(const mcmc *myMCMC) const override;
public:
  Beta_Ramp(std::size_t betaSamples,std::size_t quasiPriorSamples):
    quasiPriorSamples_(quasiPriorSamples), betaSamples_(betaSamples){}

protected:

  std::size_t quasiPriorSamples_;
  std::size_t betaSamples_;

};



class BaseInitializer: public BaseObject
{
public:
  virtual  BaseMcmcState* initialize(const mcmc* myMc,
                                     const CortexLikelihood* CL,
                                     const Parameters& initial,
                                     std::mt19937_64 mt,
                                     std::size_t* ifeval,
                                     BaseBetaManagement* beta)=0;

  virtual ~BaseInitializer(){}

  virtual BaseInitializer* create()const =0;

  static  BaseInitializer* createChild(const std::string& className);
private:
  static std::map<std::string, BaseInitializer*> childs_;
  static std::map<std::string, BaseInitializer *> getChilds();


};


class Initialize_Ensamble: public BaseInitializer
{


  // BaseClass interface
public:
  static std::string ClassName(){return "Initialize_Ensamble";}
  virtual std::__cxx11::string myClass() const override
  {
    return ClassName();
  }

  // BaseObject interface
public:

  virtual Initialize_Ensamble *create() const override
  {
    return new Initialize_Ensamble;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"radius",radius_);
    writeField(s,"numWalkers",numWalkers_);
  }
  virtual void clear() override
  {
    radius_=0;
    numWalkers_=0;
  }
  virtual bool readBody(std::__cxx11::string &line, std::istream &s) override
  {
    if (!readField(line,s,"radius",radius_))
      return false;
    else if (!readField(line,s,"numWalkers",numWalkers_))
      return false;
    else return true;
  }

protected:
  virtual void update() override
  {
  }

  // BaseInitializer interface
public:

  virtual BaseMcmcState *initialize(const mcmc* myMc,
                                    const CortexLikelihood *CL, const Parameters &initial, std::mt19937_64 mt, std::size_t *ifeval, BaseBetaManagement *beta) override;



  Initialize_Ensamble(double radius, std::size_t numWalkers):
    radius_(radius),numWalkers_(numWalkers)
  {}

  Initialize_Ensamble(){}
protected:
  double radius_;
  std::size_t numWalkers_;

private:

};




class BaseFinalize: public BaseObject
{
public:
  virtual bool meetConditions(const mcmc* myMcmc)const=0;

  virtual ~BaseFinalize(){}
  virtual BaseFinalize* create()const =0;

  static  BaseFinalize* createChild(const std::string& className);
private:
  static std::map<std::string, BaseFinalize*> childs_;
  static std::map<std::string, BaseFinalize*> getChilds();


};



class Finalize: public BaseFinalize
{


  // BaseClass interface
public:
  static std::string ClassName(){return "Finalize";}
  virtual std::__cxx11::string myClass() const override
  {
    return ClassName();
  }

  // BaseObject interface
public:
  virtual Finalize *create() const override
  {
    return new Finalize;

  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"numIter", numIter_);

    writeField(s,"maxTime", maxtime_);

  }
  virtual void clear() override
  {
    maxtime_=0;
    numIter_=0;
  }
  virtual bool readBody(std::__cxx11::string &line, std::istream &s) override
  {
    if (!readField(line,s,"numIter",numIter_)) return false;
    else if (!readField(line,s,"maxTime",maxtime_)) return false;
    return true;

  }

protected:
  virtual void update() override
  {
  }

  // BaseFinalize interface
public:
  virtual bool meetConditions(const mcmc *myMcmc) const override;

  Finalize(std::size_t numIter, double maxtime):
    numIter_(numIter), maxtime_(maxtime){}

  Finalize(){}
private:
  std::size_t numIter_;
  double maxtime_;
};



class BaseSampler: public BaseObject
{
public:
  virtual void InitializePrint(const mcmc* myMcmc)=0;

  virtual void IterateTitlesPrint(const mcmc* myMcmc
                                  , const BaseProposalDistribution* prop
                                  , const BaseMcmcState * state
                                  , const CortexLikelihood *CL)=0;
  virtual void IteratePrint(const mcmc* myMcmc)=0;

  virtual ~BaseSampler(){}
  virtual BaseSampler* create()const =0;

  static  BaseSampler* createChild(const std::string& className);
private:
  static std::map<std::string, BaseSampler*> childs_;
  static std::map<std::string, BaseSampler*> getChilds();


};


class Sampler:public BaseSampler
{


  // BaseClass interface
public:
  static std::string ClassName(){return "Sampler";}
  virtual std::__cxx11::string myClass() const override
  {
    return ClassName();

  }

  // BaseObject interface
public:
  Sampler(){}
  virtual Sampler *create() const override
  {
    return new Sampler;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
  }
  virtual void clear() override
  {
  }
  virtual bool readBody(std::__cxx11::string &line, std::istream &s) override
  {
  }

protected:
  virtual void update() override
  {
  }

  // BaseSampler interface
public:
  virtual void InitializePrint(const mcmc *myMcmc) override;

  virtual void IterateTitlesPrint(const mcmc* myMcmc
                                  , const BaseProposalDistribution* prop
                                  , const BaseMcmcState * state
                                  , const CortexLikelihood *CL) override;

  virtual void IteratePrint(const mcmc *myMcmc);
  Sampler(std::ostream* os,const std::string& fname, bool includePredictions ):
    os_(os),
    fname_(fname),
    includePredictions_(includePredictions),
    has_os_(os!=nullptr)
  {
    os_par_.open(fname_+"_par.txt",std::ofstream::out);
    os_mean_.open(fname_+"_mean.txt",std::ofstream::out);
    if (includePredictions_)
      os_pred_.open(fname_+"_pred.txt",std::ofstream::out);

  }
protected:
  std::ostream* os_;
  std::string fname_;
  bool includePredictions_;
  bool has_os_;
  std::ofstream os_par_;
  std::ofstream os_mean_;
  std::ofstream os_pred_;



};




class BaseProposalDistribution:public BaseObject
{
public:
  virtual void next(const mcmc* myMc,
                    const CortexLikelihood* CL,
                    const BaseBetaManagement* beta,
                    BaseMcmcState* s,
                    std::mt19937_64& mt_,
                    std::size_t* ifeval)=0;

  virtual void TitlesPrint(std::ostream& os)const=0;
  virtual void IteratePrint(std::ostream& os)const=0;

  static constexpr double acc_ratio_opt=0.234;
  static constexpr double acc_ratio_opt_mala=0.574;

  virtual BaseProposalDistribution* create()const =0;

  static  BaseProposalDistribution* createChild(const std::string& className);
private:
  static std::map<std::string, BaseProposalDistribution*> childs_;
  static std::map<std::string, BaseProposalDistribution*> getChilds();


};



class StretchProposal: public BaseProposalDistribution
{



  // BaseClass interface
public:
  static std::string ClassName(){return "StretchProposal";}
  virtual std::__cxx11::string myClass() const override
  {
    return ClassName();
  }

  // BaseObject interface
public:
  virtual StretchProposal *create() const override
  {
    return new StretchProposal;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"a0",a0_);
    writeField(s,"delta_0",delta_0_);
    writeField(s,"n_skip",n_skip_);
  }

  StretchProposal(){}
  virtual void clear() override
  {

  }
  virtual bool readBody(std::__cxx11::string &line, std::istream &s) override
  {
    if (!readField(line,s,"a0",a0_))
      return false;
    else if (!readField(line,s,"delta_0",delta_0_))
      return false;
    else if (!readField(line,s,"n_skip",n_skip_))
      return false;
    else return true;

  }

protected:
  virtual void update() override
  {
    a_=a0_;
    delta_=delta_0_;
    acc_ratio_min_=0;
    acc_ratio_max_=0;
  }

  // BaseProposalDistribution interface
public:
  virtual void next(const mcmc *myMc
                    , const CortexLikelihood *CL
                    , const BaseBetaManagement *beta
                    , BaseMcmcState *s
                    , std::mt19937_64& mt
                    , std::size_t *ifeval) override;


  StretchProposal(double a0,
  double delta_0,
  std::size_t n_skip):
  a0_(a0),delta_0_(delta_0),a_(a0_),delta_(delta_0),n_skip_(n_skip),acc_count_(),acc_ratio_min_(0),acc_ratio_max_(0)
  {}

  // BaseProposalDistribution interface
public:
  virtual void TitlesPrint(std::ostream &os)const override
  {
    os<<"\t"<<"a_value"<<"\t"<<"accep_ratio_min"<<"\t"<<"accep_ratio_max";
  }
  virtual void IteratePrint(std::ostream &os)const override
  {
    os<<"\t"<<a_<<"\t"<<acc_ratio_min_<<"\t"<<acc_ratio_max_;
  }


private:
double a0_;
double delta_0_;
std::size_t n_skip_;
double a_;
double delta_;
std::vector<std::size_t> acc_count_;
double acc_ratio_min_;
double acc_ratio_max_;


};


class mcmcEnsambleState: public BaseMcmcState
{
public:
  //  mcmcWalkerState(const Parameters& prior
  //                  , std::vector<std::vector<double>> trMeans
  //                  , double beta
  //                  , std::vector<double> dataLiks
  //                  , std::vector<double> logPrios);

  mcmcEnsambleState(const Parameters &prior
                    , std::size_t numWalkers
                    , double beta);



  static mcmcEnsambleState create(const CortexLikelihood* f,
                                  const Parameters& initialParam,
                                  std::size_t numWalkers,
                                  double radiusWalkers,
                                  std::mt19937_64& mt,
                                  std::size_t& ifeval,
                                  double beta);


  mcmcEnsambleState();

  const Parameters& getMeanState()const;

  std::size_t numWalkers() const;

  std::size_t numParameters()const;

  std::vector<double>& operator[](std::size_t i);

  const std::vector<double>& operator[](std::size_t i)const;

  double logBetaLik(std::size_t i)const;

  double& logPrior(std::size_t i);
  const double& logPrior(std::size_t i) const;

  void update_Mean();

  void update_Covariance();

  std::ostream& writeTrValues(std::ostream& s,std::size_t isample);




private:
  Parameters mean_;
  double dataLikMean_;
  double dataLikStd_;
  double dataLikMin_;
  double dataLikMax_;

  double priorLikMean_;
  double priorLikStd_;
  double priorLikMin_;
  double priorLikMax_;

  double beta_;
  std::vector<std::vector<std::vector<double>>> f_;
  std::vector<std::vector<double>> trMeans_;
  std::vector<double> dataLiks_;
  std::vector<double> priorLiks_;


  // BaseClass interface
public:
  static std::string ClassName(){ return "mcmcWalkerState";}
  virtual std::__cxx11::string myClass() const override{ return ClassName();}

  // BaseObject interface
public:
  virtual mcmcEnsambleState *create() const override{ return new mcmcEnsambleState;}
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override;
  virtual bool readBody(std::string &line, std::istream &s) override;

  const double &logDataLikMean() const;
  const double &logDataLikStd() const;
  const double &logPriorLikMean() const;
  const double &logPriorLikStd() const;
  std::ostream &writeValues(std::ostream &s, std::size_t isample)const;
  std::ostream &writeMeans(std::ostream &s)const;
  std::ostream &writeValuesTitles(std::ostream &s)const;
  std::ostream &writeMeansTitles(std::ostream &s)const;
  const double &logDataLikMax() const;
  const double &logDataLikMin() const;
  const double &logPriorLikMax() const;
  const double &logPriorLikMin() const;
  const double &logDataLik(std::size_t i) const;
  double &logDataLik(std::size_t i);
  double beta() const;
  std::size_t numMeasures() const;
  std::vector<std::vector<double> > &f(std::size_t i);
  std::ostream &writeYValuesTitles(std::ostream &s, const CortexLikelihood *CL)const;
  std::ostream &writeYValues(std::ostream &s, std::size_t isample, const CortexLikelihood *CL)const;
protected:
  virtual void update() override{}
};



class mcmcWrun: public BaseAgent
{
public:


  mcmcWrun(const mcmc* e,
           std::size_t numIterations):
    e_(e),
    run_(numIterations){}

  mcmcWrun(){}

  void push_back(mcmcEnsambleState w);



private:
  const mcmc* e_;
  std::vector<mcmcEnsambleState> run_;
  std::size_t iend_;

  // BaseClass interface
public:
  static std::string ClassName(){return "mcmcWRun";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual BaseObject *create() const override{

    return new mcmcWrun;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override{}
  virtual bool readBody(std::string &line, std::istream &s) override;

protected:
  virtual void update() override{}
};

class mcmc: public BaseAgent
{

public:


  mcmc(const std::string name,
       const CortexLikelihood* f,
       const Parameters& initialParam,
       BaseInitializer * initializer,
       BaseProposalDistribution* prop,
       BaseFinalize* final,
       BaseBetaManagement* betaManag,
       BaseSampler* sampler,
       std::mt19937_64::result_type initseed)
    :
      startTime_(std::chrono::steady_clock::now()),
      mt_()
    , fname_(name),
      CL_(f)
    ,initial_(initialParam)
    , initializer_(initializer),
      prop_(prop),
      final_(final),
      betaManag_(betaManag),
      sampler_(sampler)
  {

    if (initseed!=0)
      {
        mt_.seed(initseed);
        std::cout<<"Seed for random generator provided="<<initseed<<std::endl;
        if (os_mean_.is_open())
          os_mean_<<"Seed for random generator provided="<<initseed<<std::endl;
      }
    else
      {
        std::random_device rd;
        std::mt19937_64::result_type seed=rd();
        mt_.seed(seed);
        std::cout<<"Seed of random generator="<<seed<<std::endl;
        if (os_mean_.is_open())
          os_mean_<<"Seed of random generator="<<seed<<std::endl;
      }


  }

  virtual void run();
  virtual void initialize();
  virtual void iterate();
  virtual void finalize(){}


  void IterateTitlesPrint(std::ostream& s) const;
  void IteratePrint(std::ostream& s)const;


  std::size_t iter()const{return i_;}

  double runTime()const {return runt_;}

protected:
  std::chrono::steady_clock::time_point startTime_;
  std::mt19937_64 mt_;
  bool seedIsProvided_;
  std::mt19937_64::result_type seed_;

  std::string fname_;
  std::size_t i_;
  double runt_;
  double timeIter_;

  const CortexLikelihood* CL_;
  Parameters initial_;
  double acc_ratio_;
  std::size_t ifeval_;
  BaseMcmcState* ws_0_;
  BaseInitializer * initializer_;
  BaseProposalDistribution* prop_;
  BaseFinalize* final_;
  BaseBetaManagement* betaManag_;
  BaseSampler* sampler_;


  mcmc(){}









  // BaseClass interface
public:
  static std::string ClassName(){return "mcmc"; }
  virtual std::__cxx11::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual mcmc *create() const override{ return new mcmc;}
  virtual std::ostream &writeBody(std::ostream &s) const override;




  virtual void clear() override{}

  virtual bool readBody(std::string &line, std::istream &s) override;

protected:
  virtual void update() override{}

private:
  std::ofstream os_mean_;

  // Base_mcmc interface
};





#endif // BAYESITERATION

