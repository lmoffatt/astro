#include <vector>
#include <fstream>
#include <cmath>
#include <map>
#include "BayesIteration.h"
#include "MatrixInverse.h"
#include "CortexLikelihood.h"
#include "CommandManager.h"

#include <tuple>

mcmcWalkerState
nextState(const CortexLikelihood* CL,
          std::mt19937& mt,
          double a,
          mcmcWalkerState& x
          ,std::size_t& ifeval)
{

  mcmcWalkerState o(x.getMeanState(),x.numWalkers());
  std::size_t K=x.numWalkers();
  std::size_t N=x.numParameters();
#pragma omp parallel for
  for (std::size_t k=0; k<K; ++k)
    {
       std::uniform_int_distribution<std::size_t> i(0,K-2);
      auto j=i(mt);
      if (j>=k) ++j;
      std::uniform_real_distribution<double> u(0,1);
      double z=std::pow(u(mt)*(std::sqrt(a)-std::sqrt(1.0/a))+std::sqrt(1/a),2);
      std::vector<double> y=x[j]+(x[k]-x[j])*z;
      double logprior=CL->logPrior(y);
      double logpostlik=CL->logLik(y)+logprior;
      ++ifeval;
      double logq=log(z)*(N-1)+logpostlik-x.logPostLik(k);
      double logr=log(u(mt));
      if (logr<=logq)
        {
          o[k]=y;
          o.logPostLik(k)=logpostlik;
          o.logPrior(k)=logprior;
        }
      else
        {
          o[k]=x[k];
          o.logPostLik(k)=x.logPostLik(k);
          o.logPrior(k)=x.logPrior(k);

        }
    }

  o.update_Mean();
  o.update_Covariance();
  return o;
}

mcmcWalkerState::mcmcWalkerState(const Parameters &prior, std::vector<std::vector<double> > trMeans, std::vector<double> postLiks, std::vector<double> logPrios):
  mean_(prior)
,trMeans_(trMeans)
,postLiks_(postLiks)
,priorLiks_(logPrios){
  update_Mean();
  update_Covariance();

}

mcmcWalkerState::mcmcWalkerState(const Parameters& prior, std::size_t numWalkers):
  mean_(prior)
,trMeans_(numWalkers)
,postLiks_(numWalkers)
,priorLiks_(numWalkers){}

mcmcWalkerState mcmcWalkerState::create(const CortexLikelihood *f
                                        , const Parameters &initialParam
                                        , std::size_t numWalkers
                                        , double radiusWalkers
                                        , std::mt19937 &mt
                                        , std::size_t &ifeval)
{
  mcmcWalkerState o(f->getPrior(),numWalkers);
  std::size_t i=0;
  while    (i<numWalkers)
    {
      o[i]=initialParam.randomSampleValues(mt,f->getPrior(),radiusWalkers);
      double logprior=f->logPrior(o[i]);
      double logPostLik=logprior+f->logLik(o[i]);
      ++ifeval;
      o.logPostLik(i)=logPostLik;
      o.logPrior(i)=logprior;
      if (!std::isnan(logPostLik))
        ++i;

    }
  o.update_Mean();
  o.update_Covariance();
  return o;
}

mcmcWalkerState::mcmcWalkerState(){}

const Parameters &mcmcWalkerState::getMeanState() const
{
  return mean_;
}

std::size_t mcmcWalkerState::numWalkers() const
{
  return trMeans_.size();
}

std::size_t mcmcWalkerState::numParameters() const
{
  return mean_.size();
}

std::vector<double> &mcmcWalkerState::operator[](std::size_t i)
{
  return trMeans_[i];
}

const std::vector<double> &mcmcWalkerState::operator[](std::size_t i) const
{
  return trMeans_[i];
}

double &mcmcWalkerState::logPostLik(std::size_t i)
{
  return postLiks_[i];
}

const double &mcmcWalkerState::logPostLik(std::size_t i) const
{
  return postLiks_[i];
}
const double &mcmcWalkerState::logPostLikMean() const
{
  return postLikMean_;
}
const double &mcmcWalkerState::logPostLikStd() const
{
  return postLikStd_;
}
const double &mcmcWalkerState::logPriorLikMean() const
{
  return priorLikMean_;
}
const double &mcmcWalkerState::logPriorLikStd() const
{
  return priorLikStd_;
}



double &mcmcWalkerState::logPrior(std::size_t i)
{
  return priorLiks_[i];
}

const double &mcmcWalkerState::logPrior(std::size_t i) const
{
  return priorLiks_[i];
}

void mcmcWalkerState::update_Mean()
{
  std::vector<double> mean(numParameters(),0);

  for (std::size_t j=0; j<numParameters(); ++j)
    {
      for (std::size_t i=0; i<numWalkers(); ++i)
        mean[j]+=trMeans_[i][j];
      mean[j]/=numWalkers();
    }
  postLikMean_=0;
  priorLikMean_=0;
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      postLikMean_+=postLiks_[i];
      priorLikMean_+=priorLiks_[i];
    }
  postLikMean_/=numWalkers();
  priorLikMean_/=numWalkers();
  mean_.settMeans(mean);
}

void mcmcWalkerState::update_Covariance()
{
  std::vector<std::vector<double>> cov(numParameters(),std::vector<double>(numParameters(),0));
  for (std::size_t k=0; k<numWalkers(); ++k)
    {
      for (std::size_t i=0; i<numParameters(); ++i)
        for (std::size_t j=i; j<numParameters(); ++j)
          cov[i][j]+=(trMeans_[k][i]-mean_[i])*(trMeans_[k][j]-mean_[j]);
    }
  for (std::size_t i=0; i<numParameters(); ++i)
    for (std::size_t j=i; j<numParameters(); ++j)
      {
        cov[i][j]/=numWalkers();
        cov[j][i]=cov[i][j];
      }

  postLikStd_=0;
  priorLikStd_=0;
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      postLikStd_+=std::pow(postLiks_[i]-postLikMean_,2);
      priorLikStd_+=std::pow(priorLiks_[i]-priorLikMean_,2);
    }

  postLikStd_/=numWalkers();
  postLikStd_=std::sqrt(postLikStd_);
  priorLikStd_/=numWalkers();
  priorLikStd_=std::sqrt(priorLikStd_);



  mean_.setCovariance(cov);
}

std::ostream &mcmcWalkerState::writeTrValues(std::ostream &s, std::size_t isample)
{
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      s<<isample<<"\t"<<i<<"\t"<<logPostLik(i)<<"\t"<<logPrior(i);
      for (std::size_t k=0; k<numParameters(); ++k)
        s<<"\t"<<trMeans_[i][k];
      s<<"\n";
    }
  return s;
}

std::ostream &mcmcWalkerState::writeValues(std::ostream &s, std::size_t isample)
{
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      s<<isample<<"\t"<<i<<"\t"<<logPostLik(i)<<"\t"<<logPrior(i);
      for (std::size_t k=0; k<numParameters(); ++k)
        s<<"\t"<<mean_.getTransform(k)->inverse(trMeans_[i][k]);
      s<<"\n";
    }
  return s;
}


std::ostream &mcmcWalkerState::writeValuesTitles(std::ostream &s)
{
  s<<"isample"<<"\t"<<"i"<<"\t"<<"logPostLik(i)"<<"\t"<<"logPrior(i)";
  for (std::size_t k=0; k<numParameters(); ++k)
    s<<"\t"<<mean_.indexToName(k);
  s<<"\n";
  return s;
}

std::ostream &mcmcWalkerState::writeMeans(std::ostream &s)
{
    s<<"\t"<<logPostLikMean()<<"\t"<<logPostLikStd();
    s<<"\t"<<logPriorLikMean()<<"\t"<<logPriorLikStd();

      for (std::size_t k=0; k<numParameters(); ++k)
      {
          s<<"\t"<<mean_.mean(k);
          s<<"\t"<<mean_.dBStd(k);
      }
      s<<"\n";
  return s;
}

std::ostream &mcmcWalkerState::writeMeansTitles(std::ostream &s)
{
    s<<"\t"<<"logPostLikMean"<<"\t"<<"logPostLikStd";
    s<<"\t"<<"logPriorLikMean"<<"\t"<<"logPriorLikStd";

      for (std::size_t k=0; k<numParameters(); ++k)
      {
          s<<"\t"<<mean_.indexToName(k)<<"_mean";
          s<<"\t"<<mean_.indexToName(k)<<"_dB";
        }
      s<<"\n";
  return s;
}




std::ostream &mcmcWalkerState::writeBody(std::ostream &s) const
{
  writeField(s,"mean_state",mean_);
  writeField(s,"transformed_Means",trMeans_);
  writeField(s,"postLikelihood",postLiks_);
  writeField(s,"priorLikelihood",priorLiks_);
  return s;
}

void mcmcWalkerState::clear() {
  priorLiks_.clear();
  postLiks_.clear();
  trMeans_.clear();
}

bool mcmcWalkerState::readBody(std::string &line, std::istream &s)
{
  if (!readField(line,s,"mean_state",mean_)) return false;
  else if (!readField(line,s,"transformed_Means",trMeans_)) return false;
  else if (!readField(line,s,"postLikelihood",postLiks_)) return false;
  else if (!readField(line,s,"priorLikelihood",priorLiks_)) return false;
  else return true;
}








emcee_mcmc::emcee_mcmc(const CortexLikelihood *f
                       , const Parameters &initialParam
                       , double maxDuration_minutes
                       , std::size_t numIterations
                       , std::size_t numWalkers
                       , std::size_t nSkip
                       , double radiusWalkers
                       , const std::string &name
                       , double a
                       , std::mt19937::result_type initseed)
  :
    startTime_(std::chrono::steady_clock::now()),
    mt_()
  , fname_(name),
    os_val_(),
    os_mean_(),
    CL_(f)
  ,initial_(initialParam)
  , maxDuration_minutes_(maxDuration_minutes)
  , numSamples_(numIterations)
  , numWalkers_(numWalkers)
  ,n_skip_(nSkip)
  , radiusWalkers_(radiusWalkers)
  , a_(a)
{

  fname_=getSaveName(fname_);
  os_val_.open(fname_+"_val.txt",std::ofstream::out);
  os_mean_.open(fname_+"_mean.txt",std::ofstream::out);


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
      std::mt19937::result_type seed=rd();
      mt_.seed(seed);
      std::cout<<"Seed of random generator="<<seed<<std::endl;
      if (os_mean_.is_open())
        os_mean_<<"Seed of random generator="<<seed<<std::endl;
    }

}

mcmcWrun emcee_mcmc::run()
{
  startTime_=std::chrono::steady_clock::now();
  ifeval_=0;

  mcmcWrun o(this,numSamples_);
  mcmcWalkerState ws_0=mcmcWalkerState::create(CL_,initial_,numWalkers_,radiusWalkers_,mt_,ifeval_);
  mcmcWalkerState ws_1;
  std::size_t i=0;
  double runt=0;
  os_mean_.precision(10);
  os_val_.precision(10);
  std::cout<<"i"<<"\t"<<"ifeval"<<"\t"<<"runt"<<"\t"<<"timeIter"<<"\t";
  std::cout<<"logPostLikMean()"<<"\t"<<"logPostLikStd()";
  std::cout<<"\t"<<"logPriorLikMean()"<<"\t"<<"logPriorLikStd()"<<std::endl;

  os_mean_<<"i"<<"\t"<<"ifeval"<<"\t"<<"runt"<<"\t"<<"timeIter";
  ws_0.writeMeansTitles(os_mean_);

  while ((runt<maxDuration_minutes_)&&(i<numSamples_))
    {
      auto tnow=std::chrono::steady_clock::now();
      auto d=tnow-startTime_;
      double t0=runt;
      runt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
      double timeIter=60*(runt-t0);
      for (std::size_t ii=0; ii<n_skip_; ++ii)
        {
          ws_1=nextState(CL_,mt_,a_,ws_0,ifeval_);
          std::swap(ws_0,ws_1);
        }

      std::cout<<i<<"\t"<<ifeval_<<"\t"<<runt<<"\t"<<timeIter<<"\t";
      std::cout<<ws_0.logPostLikMean()<<"\t"<<ws_0.logPostLikStd();
      std::cout<<"\t"<<ws_0.logPriorLikMean()<<"\t"<<ws_0.logPriorLikStd()<<std::endl;
      ws_0.writeValuesTitles(os_val_);
      ws_0.writeValues(os_val_,i);
      os_mean_<<i<<"\t"<<ifeval_<<"\t"<<runt<<"\t"<<timeIter;
      ws_0.writeMeans(os_mean_);
      os_mean_.flush();
      os_val_.flush();
      ++i;
    }
}

std::ostream &emcee_mcmc::writeBody(std::ostream &s) const
{
  writeField(s,"Likelihood_Model",CL_->id());
  writeField(s,"Initial_Parameter",initial_);
  writeField(s,"Maximal_duration",maxDuration_minutes_);
  writeField(s,"Number_Samples",numSamples_);
  writeField(s,"Number_Walkers",numWalkers_);
  writeField(s,"Initial_Radius",radiusWalkers_);
  writeField(s,"eemcee_a",a_);

  return s;
}

bool emcee_mcmc::readBody(std::string &line, std::istream &s)
{
  std::string likelihoodName;
  if (!readField(line,s,"Likelihood_Model",likelihoodName)) return false;
  else
    {
      CL_=cm_->getLikelihood(likelihoodName);
      if (!readField(line,s,"Initial_Parameter",initial_))
        return false;
      else if (!readField(line,s,"Maximal_duration",maxDuration_minutes_))
        return false;
      else if (!readField(line,s,"Initial_Parameter",initial_))
        return false;
      else if (!readField(line,s,"Number_Samples",numSamples_))
        return false;
      else if (!readField(line,s,"Number_Walkers",numWalkers_))
        return false;
      else if (!readField(line,s,"Initial_Radius",radiusWalkers_))
        return false;
      else if (!readField(line,s,"eemcee_a",a_))
        return false;
      else return true;
    }
}



void mcmcWrun::push_back(mcmcWalkerState w)
{
  run_[iend_]=w;
  iend_++;
}

std::ostream &mcmcWrun::writeBody(std::ostream &s) const
{

  writeField(s,"mcmc_Algorithm",e_->id());
  writeField(s,"mcmc_run",run_);

}

bool mcmcWrun::readBody(std::string &line, std::istream &s)
{
  std::string algorithm;
  if (!readField(line,s,"mcmc_Algorithm",algorithm)) return false;
  else {
      auto e=cm_->getMcmc(algorithm);
      if (e==nullptr) return false;
      else
        {
          e_=e;
          if (!readField(line,s,"mcmc_run",run_)) return false;
          else return true;

        }
    }
}
