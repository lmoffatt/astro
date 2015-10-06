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
  mcmcWalkerState(const Parameters& prior
                  , std::vector<std::vector<double>> trMeans
                  ,std::vector<double> postLiks
                  , std::vector<double> logPrios);

  mcmcWalkerState(const Parameters &prior
                  , std::size_t numWalkers);



  static mcmcWalkerState create(const CortexLikelihood* f,
                                const Parameters& initialParam,
                                std::size_t numWalkers,
                                double radiusWalkers,
                                std::mt19937& mt,
                                std::size_t& ifeval);


  mcmcWalkerState();

  const Parameters& getMeanState()const;

  std::size_t numWalkers() const;

  std::size_t numParameters()const;

  std::vector<double>& operator[](std::size_t i);

  const std::vector<double>& operator[](std::size_t i)const;

  double& logPostLik(std::size_t i);
  const double& logPostLik(std::size_t i)const;

  double& logPrior(std::size_t i);
  const double& logPrior(std::size_t i) const;

  void update_Mean();

  void update_Covariance();

  std::ostream& writeTrValues(std::ostream& s,std::size_t isample);




private:
  Parameters mean_;
  double postLikMean_;
  double postLikStd_;
  double priorLikMean_;
  double priorLikStd_;
  std::vector<std::vector<double>> trMeans_;
  std::vector<double> postLiks_;
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

  const double &logPostLikMean() const;
  const double &logPostLikStd() const;
  const double &logPriorLikMean() const;
  const double &logPriorLikStd() const;
  std::ostream &writeValues(std::ostream &s, std::size_t isample);
  std::ostream &writeMeans(std::ostream &s);
  std::ostream &writeValuesTitles(std::ostream &s);
  std::ostream &writeMeansTitles(std::ostream &s);
protected:
  virtual void update() override{}
};


class emcee_mcmc;
class mcmcWrun: public BaseAgent
{
public:


  mcmcWrun(const emcee_mcmc* e,
           std::size_t numIterations):
    e_(e),
    run_(numIterations){}

  mcmcWrun(){}

  void push_back(mcmcWalkerState w);



private:
  const emcee_mcmc* e_;
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

class emcee_mcmc: public BaseAgent
{
public:
  emcee_mcmc(const CortexLikelihood* f,
             const Parameters& initialParam,
             double maxDuration_minutes,
             std::size_t numIterations,
             std::size_t numWalkers,
             std::size_t nSkip,
             double radiusWalkers,
             const std::string&  name,
             double a,
             std::mt19937::result_type initseed);

  emcee_mcmc(){}

  mcmcWrun run();


  void next();


  // BaseClass interface
public:
  static std::string ClassName(){return "emcee_mcmc"; }
  virtual std::__cxx11::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual BaseObject *create() const override{ return new emcee_mcmc;}
  virtual std::ostream &writeBody(std::ostream &s) const override;




  virtual void clear() override{}

  virtual bool readBody(std::string &line, std::istream &s) override;

protected:
  virtual void update() override{}

private:
  std::chrono::steady_clock::time_point startTime_;
  std::mt19937 mt_;
  std::string fname_;
  std::ofstream os_val_;
  std::ofstream os_mean_;
  const CortexLikelihood* CL_;
  Parameters initial_;
  double maxDuration_minutes_;
  std::size_t numSamples_;
  std::size_t numWalkers_;
  std::size_t n_skip_;
  std::size_t ifeval_;
  double radiusWalkers_;
  double a_;
};


#endif // BAYESITERATION

