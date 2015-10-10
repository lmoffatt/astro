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



class mcmcWalkerState: public BaseObject
{
public:
//  mcmcWalkerState(const Parameters& prior
//                  , std::vector<std::vector<double>> trMeans
//                  , double beta
//                  , std::vector<double> dataLiks
//                  , std::vector<double> logPrios);

  mcmcWalkerState(const Parameters &prior
                  , std::size_t numWalkers
                  , double beta);



  static mcmcWalkerState create(const CortexLikelihood* f,
                                const Parameters& initialParam,
                                std::size_t numWalkers,
                                double radiusWalkers,
                                std::mt19937& mt,
                                std::size_t& ifeval,
                                double beta);


  mcmcWalkerState();

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
  virtual BaseObject *create() const override{ return new mcmcWalkerState;}
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override;
  virtual bool readBody(std::string &line, std::istream &s) override;

  const double &logDataLikMean() const;
  const double &logDataLikStd() const;
  const double &logPriorLikMean() const;
  const double &logPriorLikStd() const;
  std::ostream &writeValues(std::ostream &s, std::size_t isample);
  std::ostream &writeMeans(std::ostream &s);
  std::ostream &writeValuesTitles(std::ostream &s);
  std::ostream &writeMeansTitles(std::ostream &s);
  const double &logDataLikMax() const;
  const double &logDataLikMin() const;
  const double &logPriorLikMax() const;
  const double &logPriorLikMin() const;
  const double &logDataLik(std::size_t i) const;
  double &logDataLik(std::size_t i);
  double beta() const;
  std::size_t numMeasures() const;
  std::vector<std::vector<double> > &f(std::size_t i);
  std::ostream &writeYValuesTitles(std::ostream &s, const CortexLikelihood *CL);
  std::ostream &writeYValues(std::ostream &s, std::size_t isample, const CortexLikelihood *CL);
protected:
  virtual void update() override{}
};


class mcmc;
class mcmcWrun: public BaseAgent
{
public:


  mcmcWrun(const mcmc* e,
           std::size_t numIterations):
    e_(e),
    run_(numIterations){}

  mcmcWrun(){}

  void push_back(mcmcWalkerState w);



private:
  const mcmc* e_;
  std::vector<mcmcWalkerState> run_;
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
  enum method {STRETCH,WALK};


  mcmc(const CortexLikelihood* f,
       const Parameters& initialParam,
       double maxDuration_minutes,
       std::size_t betaSamples,
       std::size_t eqSamples,
       std::size_t numWalkers,
       std::size_t nSkip,
       double radiusWalkers,
       const std::string&  name,
       double amin, double amax,
       std::size_t N_for_Walk,
       double rWalk,
       method m, bool includePredictions,
       std::mt19937::result_type initseed);

  mcmc(){}

  void run();




  void next();




  // BaseClass interface
public:
  static std::string ClassName(){return "mcmc"; }
  virtual std::__cxx11::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual BaseObject *create() const override{ return new mcmc;}
  virtual std::ostream &writeBody(std::ostream &s) const override;




  virtual void clear() override{}

  virtual bool readBody(std::string &line, std::istream &s) override;

protected:
  virtual void update() override{}

private:
  std::chrono::steady_clock::time_point startTime_;
  std::mt19937 mt_;
  std::string fname_;
  std::ofstream os_par_;
  std::ofstream os_mean_;
  std::ofstream os_pred_;

  const CortexLikelihood* CL_;
  Parameters initial_;
  double betarun_;
  double maxDuration_minutes_;
  std::size_t betaSamples_;
  std::size_t numSamples_;
  std::size_t numWalkers_;
  std::size_t n_skip_;
  double acc_ratio_;
  std::size_t ifeval_;
  double radiusWalkers_;
  double amin_;
  double amax_;

  std::size_t N_for_Walk_;
  double rWalk_;
  bool includePredictions_;
  method method_;
};


#endif // BAYESITERATION

