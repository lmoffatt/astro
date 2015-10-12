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
nextStateStretch(const CortexLikelihood* CL,
                 std::mt19937& mt,
                 double amin,
                 double amax,
                 double beta,
                 mcmcWalkerState& x
                 ,std::size_t& ifeval
                 ,std::size_t& accept_count)
{

 // mcmcWalkerState o(x.getMeanState(),x.numWalkers(),beta);
  std::size_t K=x.numWalkers();
  std::size_t N=x.numParameters();

  std::vector<std::size_t> js(K);

  std::uniform_int_distribution<std::size_t> is(0,K/2-1);
  for (std::size_t k=0; k<K/2; ++k)
    {
      js[k]=is(mt)+K/2;
      js[k+K/2]=is(mt);

    }
  std::uniform_real_distribution<double> u(0,1);
  std::vector<double> rs(K);
  for (std::size_t k=0; k<K; ++k)
    {
      rs[k]=u(mt);
    }

  std::vector<double> a;
  if (amin==amax)
    a.push_back(amin);
  else
    {
      a.resize(K);
      std::uniform_real_distribution<double> arand(amin,amax);
      for (std::size_t k=0; k<K; ++k)
        a[k]=arand(mt);
  }


// first half
for (std::size_t ih=0; ih<2; ++ih)
{
#pragma omp parallel for default (none) shared (ih,K,CL,rs,a,amin,amax,x,N,beta,js) reduction (+ : ifeval, accept_count)
  for (std::size_t k=ih*K/2; k<K/2+ih*K/2; ++k)
    {
      std::size_t j=js[k];
      double av;
      if (amin==amax)
        av=a[0];
      else
        av=a[k];
      double z=std::pow(rs[k]*(std::sqrt(av)-std::sqrt(1.0/av))+std::sqrt(1/av),2);
      std::vector<double> y=x[j]+(x[k]-x[j])*z;
      double logprior=CL->logPrior(y);
      std::vector<std::vector<double>> f=CL->f(y);
      double logdatalik=CL->logLik(f);
      ++ifeval;
      double logq=log(z)*(N-1)+beta*logdatalik+logprior-beta*x.logDataLik(k)-x.logPrior(k);
      double logr=log(rs[k]);
      if (logr<=logq)
        {
          x[k]=y;
          x.f(k)=f;
          x.logDataLik(k)=logdatalik;
          x.logPrior(k)=logprior;
          accept_count++;
        }
      }

}
  x.update_Mean();
  x.update_Covariance();
  return x;
}



mcmcWalkerState
nextStateStretch_2(const CortexLikelihood* CL,
                   std::mt19937& mt,
                   double a,
                   double beta,
                   mcmcWalkerState& x
                   ,std::size_t& ifeval
                   ,std::size_t& accept_count)
{

  mcmcWalkerState o(x.getMeanState(),x.numWalkers(),beta);
  std::size_t K=x.numWalkers();
  std::size_t N=x.numParameters();

  std::vector<std::size_t> js(K);
  std::uniform_int_distribution<std::size_t> is(0,K-2);
  for (std::size_t k=0; k<K; ++k)
    {
      js[k]=is(mt);
      if (js[k]>=k) ++js[k];
    }
  std::uniform_real_distribution<double> u(0,1);
  std::vector<double> rs(K);
  for (std::size_t k=0; k<K; ++k)
    {
      rs[k]=u(mt);
    }



  std::vector<double> z(K);
  std::vector<std::vector<double>> y(K);
  std::vector<double> logprior(K);

  std::vector<Parameters> yp(K);
  for (std::size_t k=0; k<K; ++k)
    {
      std::size_t j=js[k];
      z[k]=std::pow(rs[k]*(std::sqrt(a)-std::sqrt(1.0/a))+std::sqrt(1/a),2);
      y[k]=x[j]+(x[k]-x[j])*z[k];
      yp[k]=CL->getPrior().toParameters(y[k]);
    }


#pragma omp parallel for default (none) shared (o,K,CL,yp,logprior)
  for (std::size_t k=0; k<K; ++k)
    o.f(k)=CL->f(yp[k]);


  for (std::size_t k=0; k<K; ++k)
    {
      logprior[k]=CL->logPrior(y[k]);
      double logdatalik=CL->logLik(o.f(k));
      ++ifeval;
      double logq=log(z[k])*(N-1)+beta*logdatalik+logprior[k]-beta*x.logDataLik(k)-x.logPrior(k);
      double logr=log(rs[k]);
      if (logr<=logq)
        {
          o[k]=y[k];
          o.logDataLik(k)=logdatalik;
          o.logPrior(k)=logprior[k];
          accept_count++;
        }
      else
        {
          o[k]=x[k];
          o.f(k)=x.f(k);
          o.logDataLik(k)=x.logDataLik(k);
          o.logPrior(k)=x.logPrior(k);
        }
    }

  o.update_Mean();
  o.update_Covariance();
  return o;
}






mcmcWalkerState
nextStateWalk(const CortexLikelihood* CL,
              std::mt19937& mt,
              std::size_t N,
              double rWalk,
              double beta,

              mcmcWalkerState& x
              ,std::size_t& ifeval
              ,std::size_t& accept_count)
{
  mcmcWalkerState o(x.getMeanState(),x.numWalkers(),beta);
  std::size_t K=x.numWalkers();

  std::vector<std::vector<std::size_t>> js(K,std::vector<std::size_t>(N));

  std::uniform_int_distribution<std::size_t> is(0,K-1);
  for (std::size_t k=0; k<K; ++k)
    {
      std::vector<bool> sel(k,false);
      sel[k]=true;
      std::size_t jsel=k;
      std::size_t nsel=0;
      while (nsel<N)
        {
          while (sel[jsel])
            jsel=is(mt);
          js[k][nsel]=jsel;
          sel[jsel]=true;
          ++nsel;
        }
    }
  std::vector<std::vector<double>> zs(K,std::vector<double>(N));

  std::normal_distribution<double> n;
  for (std::size_t k=0; k<K; ++k)
    for (std::size_t i=0; i<N; ++i)
      zs[k][i]=n(mt);

  std::uniform_real_distribution<double> u(0,1);
  std::vector<double> rs(K);
  for (std::size_t k=0; k<K; ++k)
    {
      rs[k]=u(mt);
    }


#pragma omp parallel for default (none) shared (o,K,CL,rs,x,N,beta,js,rWalk,zs) reduction (+ : ifeval,accept_count)
  for (std::size_t k=0; k<K; ++k)
    {
      std::vector<double> xmean=x[js[k][0]];
      for (std::size_t jj=1; jj<N; ++jj)
        {
          xmean+=x[js[k][jj]];
        }
      xmean*=1.0/N;
      std::vector<double> y=x[k];
      for (std::size_t jj=0; jj<N; ++jj)
        {
          y+=(x[js[k][jj]]-xmean)*zs[k][jj]*rWalk;
        }
      double logprior=CL->logPrior(y);
      double logdatalik=CL->logLik(y);
      ++ifeval;
      double logq=beta*logdatalik+logprior-beta*x.logDataLik(k)-x.logPrior(k);
      double logr=log(rs[k]);
      if (logr<=logq)
        {
          o[k]=y;
          o.logDataLik(k)=logdatalik;
          o.logPrior(k)=logprior;
          accept_count++;
        }
      else
        {
          o[k]=x[k];
          o.logDataLik(k)=x.logDataLik(k);
          o.logPrior(k)=x.logPrior(k);

        }
    }
  o.update_Mean();
  o.update_Covariance();
  return o;
}



//mcmcWalkerState::mcmcWalkerState(const Parameters &prior,
//                                 std::vector<std::vector<double> > trMeans,
//                                 double beta,
//                                 std::vector<double> dataLiks,
//                                 std::vector<double> logPrios):
//  mean_(prior)
//,beta_(beta)
//,trMeans_(trMeans)
//,dataLiks_(dataLiks)
//,priorLiks_(logPrios)
//{
//  update_Mean();
//  update_Covariance();
//}

mcmcWalkerState::mcmcWalkerState(const Parameters& prior, std::size_t numWalkers, double beta):
  mean_(prior)
,beta_(beta)
,f_(numWalkers)
,trMeans_(numWalkers)
,dataLiks_(numWalkers)
,priorLiks_(numWalkers)
{}

mcmcWalkerState mcmcWalkerState::create(const CortexLikelihood *f
                                        , const Parameters &initialParam
                                        , std::size_t numWalkers
                                        , double radiusWalkers
                                        , std::mt19937 &mt
                                        , std::size_t &ifeval,
                                        double beta)
{
  mcmcWalkerState o(f->getPrior(),numWalkers,beta);
  std::size_t i=0;
  while    (i<numWalkers)
    {
      o[i]=initialParam.randomSampleValues(mt,f->getPrior(),radiusWalkers);
      double logprior=f->logPrior(o[i]);
      double logDataLik=f->logLik(o[i]);
      ++ifeval;
      o.logDataLik(i)=logDataLik;
      o.logPrior(i)=logprior;
      if (!std::isnan(logDataLik))
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

std::size_t mcmcWalkerState::numMeasures() const
{
  return f_[0].size();
}


std::vector<double> &mcmcWalkerState::operator[](std::size_t i)
{
  return trMeans_[i];
}

std::vector<std::vector<double>> &mcmcWalkerState::f(std::size_t i)
{
  return f_[i];
}


const std::vector<double> &mcmcWalkerState::operator[](std::size_t i) const
{
  return trMeans_[i];
}

double &mcmcWalkerState::logDataLik(std::size_t i)
{
  return dataLiks_[i];
}

const double &mcmcWalkerState::logDataLik(std::size_t i)const
{
  return dataLiks_[i];
}

double mcmcWalkerState::logBetaLik(std::size_t i) const
{
  return beta_*logDataLik(i)+logPrior(i);
}
const double &mcmcWalkerState::logDataLikMean() const
{
  return dataLikMean_;
}

const double &mcmcWalkerState::logPriorLikMax() const
{
  return priorLikMax_;
}

const double &mcmcWalkerState::logPriorLikMin() const
{
  return priorLikMin_;
}


const double &mcmcWalkerState::logDataLikMax() const
{
  return dataLikMax_;
}

const double &mcmcWalkerState::logDataLikMin() const
{
  return dataLikMin_;
}

double mcmcWalkerState::beta()const
{
  return beta_;
}


const double &mcmcWalkerState::logDataLikStd() const
{
  return dataLikStd_;
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
  dataLikMean_=0;
  priorLikMean_=0;
  dataLikMax_=dataLiks_[0];
  dataLikMin_=dataLikMax_;
  priorLikMax_=priorLiks_[0];
  priorLikMin_=priorLikMax_;

  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      dataLikMean_+=dataLiks_[i];
      priorLikMean_+=priorLiks_[i];
      if (dataLiks_[i]>dataLikMax_) dataLikMax_=dataLiks_[i];
      if (dataLiks_[i]<dataLikMin_) dataLikMin_=dataLiks_[i];
      if (priorLiks_[i]>priorLikMax_) priorLikMax_=priorLiks_[i];
      if (priorLiks_[i]<priorLikMin_) priorLikMin_=priorLiks_[i];
    }
  dataLikMean_/=numWalkers();
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

  dataLikStd_=0;
  priorLikStd_=0;
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      dataLikStd_+=std::pow(dataLiks_[i]-dataLikMean_,2);
      priorLikStd_+=std::pow(priorLiks_[i]-priorLikMean_,2);
    }

  dataLikStd_/=numWalkers();
  dataLikStd_=std::sqrt(dataLikStd_);
  priorLikStd_/=numWalkers();
  priorLikStd_=std::sqrt(priorLikStd_);



  mean_.setCovariance(cov);
}

std::ostream &mcmcWalkerState::writeTrValues(std::ostream &s, std::size_t isample)
{
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      s<<"isamp="<<isample<<"\t"<<beta()<<"\t"<<i<<"\t"<<logDataLik(i)<<"\t"<<logPrior(i);
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
      s<<"isamp="<<isample<<"\t"<<beta()<<"\t"<<i<<"\t"<<logDataLik(i)<<"\t"<<logPrior(i);
      for (std::size_t k=0; k<numParameters(); ++k)
        s<<"\t"<<mean_.getTransform(k)->inverse(trMeans_[i][k]);
      s<<"\n";
    }
  return s;
}

std::ostream &mcmcWalkerState::writeYValues(std::ostream &s, std::size_t isample, const CortexLikelihood* CL)
{
  for (std::size_t i=0; i<numWalkers(); ++i)
    {
      s<<isample<<"\t"<<beta_<<"\t"<<i<<"\t"<<logDataLik(i)<<"\t"<<logPrior(i);

      std::size_t numMeasures=CL->getExperiment()->numMeasures();
      std::size_t js=0;

      for (std::size_t im=0; im<numMeasures; ++im)
        {
          s<<"\t"<<"\t"<<CL->getExperiment()->getMeasure(im)->id();
          if (CL->getExperiment()->getMeasure(im)->inj_Width()>0)
            {
              s<<"\t"<<"injury";
              for (std::size_t is=0; is<CL->getModel()->getNumberOfSimulatedStates(); ++is)
                {
                  s<<"\t"<<f_[i][js][0];
                  ++js;
                }
            }

          for (std::size_t ix=0; ix<CL->getExperiment()->getMeasure(im)->meanAstro().size(); ++ix)
            {
              s<<"\t";
              for (std::size_t is=0; is<f_[i][0].size(); ++is)
                {
                  s<<"\t"<<f_[i][js][is];
                }
              ++js;
            }
        }
      s<<"\n";
    }
  return s;
}
std::ostream &mcmcWalkerState::writeYValuesTitles(std::ostream &s, const CortexLikelihood * CL)
{
  s<<"isample"<<"\t"<<"beta"<<"\t"<<"i"<<"\t"<<"logDataLik(i)"<<"\t"<<"logPrior(i)";

  std::size_t numMeasures=CL->getExperiment()->numMeasures();
  for (std::size_t i=0; i<numMeasures; ++i)
    {
      s<<"\t"<<"\t"<<CL->getExperiment()->getMeasure(i)->id();
      if (CL->getExperiment()->getMeasure(i)->inj_Width()>0)
        {
          s<<"\t"<<"injury";
          for (std::size_t is=0; is<CL->getModel()->getNumberOfSimulatedStates(); ++is)
            s<<"\t"<<"S_INJ_"<<is;
        }

      for (std::size_t ix=0; ix<CL->getExperiment()->getMeasure(i)->meanAstro().size(); ++ix)
        {
          s<<"\t";
          for (std::size_t is=0; is<f_[i][0].size(); ++is)
            {
              s<<"\t"<<"n_pred_"<<is<<"_at_"<<CL->getExperiment()->getMeasure(i)->xpos()[ix];
            }
        }
    }
  s<<"\n";
  return s;
}



std::ostream &mcmcWalkerState::writeValuesTitles(std::ostream &s)
{
  s<<"isample"<<"\t"<<"beta"<<"\t"<<"i"<<"\t"<<"logDataLik(i)"<<"\t"<<"logPrior(i)";
  for (std::size_t k=0; k<numParameters(); ++k)
    s<<"\t"<<mean_.indexToName(k);
  s<<"\n";
  return s;
}

std::ostream &mcmcWalkerState::writeMeans(std::ostream &s)
{
  s<<"\t"<<beta();
  s<<"\t"<<logDataLikMean()<<"\t"<<logDataLikStd()<<"\t"<<logDataLikMax()<<"\t"<<logDataLikMin();
  s<<"\t"<<logPriorLikMean()<<"\t"<<logPriorLikStd()<<"\t"<<logPriorLikMax()<<"\t"<<logPriorLikMin();

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
  s<<"\t"<<"beta";
  s<<"\t"<<"DataLikMean"<<"\t"<<"DataLikStd"<<"\t"<<"DataLikMax"<<"\t"<<"DataLikMin";
  s<<"\t"<<"PriorLikMean"<<"\t"<<"PriorLikStd"<<"\t"<<"PriorLikMax"<<"\t"<<"PriorLikMin";

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
  writeField(s,"beta",beta_);
  writeField(s,"transformed_Means",trMeans_);
  writeField(s,"dataLikelihood",dataLiks_);
  writeField(s,"priorLikelihood",priorLiks_);
  return s;
}

void mcmcWalkerState::clear() {
  priorLiks_.clear();
  dataLiks_.clear();
  trMeans_.clear();
}

bool mcmcWalkerState::readBody(std::string &line, std::istream &s)
{
  if (!readField(line,s,"mean_state",mean_)) return false;
  else if (!readField(line,s,"beta",beta_)) return false;
  else if (!readField(line,s,"transformed_Means",trMeans_)) return false;
  else if (!readField(line,s,"dataLikelihood",dataLiks_)) return false;
  else if (!readField(line,s,"priorLikelihood",priorLiks_)) return false;
  else return true;
}








mcmc::mcmc(const CortexLikelihood *f
           , const Parameters &initialParam
           , double maxDuration_minutes
           , std::size_t betaSamples
           , std::size_t eqSamples
           , std::size_t numWalkers
           , std::size_t nSkip
           , double radiusWalkers
           , const std::string &name
           , double amin
           , double amax
           ,std::size_t N_for_Walk
           ,double rWalk
           ,method m
           , bool includePredictions
           , std::mt19937::result_type initseed)
  :
    startTime_(std::chrono::steady_clock::now()),
    mt_()
  , fname_(name),
    os_par_(),
    os_mean_(),
    os_pred_(),
    CL_(f)
  ,initial_(initialParam)
  , maxDuration_minutes_(maxDuration_minutes)
  , betaSamples_(betaSamples)
  , numSamples_(betaSamples+eqSamples)
  , numWalkers_(numWalkers)
  ,n_skip_(nSkip)
  ,acc_ratio_()
  , radiusWalkers_(radiusWalkers)
  ,amin_(amin)
  ,amax_(amax)
  ,N_for_Walk_(N_for_Walk)
  ,rWalk_(rWalk)
  ,includePredictions_(includePredictions)
  ,method_(m)
{

  fname_=getSaveName(fname_);
  os_par_.open(fname_+"_par.txt",std::ofstream::out);
  os_mean_.open(fname_+"_mean.txt",std::ofstream::out);
  if (includePredictions_)
    os_pred_.open(fname_+"_pred.txt",std::ofstream::out);


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

void mcmc::run()
{
  startTime_=std::chrono::steady_clock::now();
  ifeval_=0;
  betarun_=0;
  mcmcWalkerState ws_0=mcmcWalkerState::create(CL_,initial_,numWalkers_,radiusWalkers_,mt_,ifeval_,betarun_);
  std::size_t i=0;
  double runt=0;
  os_mean_.precision(10);
  os_par_.precision(10);
  if (includePredictions_)
    os_pred_.precision(10);

  std::cout<<"mcmc run. method=";
  os_mean_<<"mcmc run. method=";
  switch (method_) {
    case WALK:
      std::cout
          <<"walk. Number of points="<<N_for_Walk_<<" walk ratio="<<rWalk_<<std::endl;
      os_mean_
          <<"walk. Number of points="<<N_for_Walk_<<" walk ratio="<<rWalk_<<std::endl;

      break;
    case STRETCH:
      std::cout
          <<"stretch. value of parameter amin="<<amin_<<" amax="<<amax_<<std::endl;
      os_mean_
          <<"stretch. value of parameter amin="<<amin_<<" amax="<<amax_<<std::endl;
      break;

    default:
      break;
    }


  std::cout<<"i"<<"\t"<<"beta"<<"\t"<<"ifeval"<<"\t"<<"runt"<<"\t"<<"timeIter"<<"\t"<<"acc_r"<<"\t";
  std::cout<<"DataLikMean()"<<"\t"<<"DataLikStd()"<<"\t";
  std::cout<<"DataLikMax()"<<"\t"<<"DataLikMin()";
  std::cout<<"\t"<<"PriorLikMean()"<<"\t"<<"PriorLikStd()";
  std::cout<<"PriorLikMax()"<<"\t"<<"PriorLikMin()"<<std::endl;

  os_mean_<<"i"<<"\t"<<"ifeval"<<"\t"<<"runt"<<"\t"<<"timeIter"<<"\t"<<"acc_r";
  ws_0.writeMeansTitles(os_mean_);

  while ((runt<maxDuration_minutes_)&&(i<numSamples_))
    {
      betarun_=std::min(1.0*i/betaSamples_,1.0);
      auto tnow=std::chrono::steady_clock::now();
      auto d=tnow-startTime_;
      double t0=runt;
      runt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
      double timeIter=60*(runt-t0);
      std::size_t acc_count=0;
      for (std::size_t ii=0; ii<n_skip_; ++ii)
        {
          //          ws_1=nextState(CL_,mt_,a_,betarun_,ws_0,ifeval_,acc_count);
          //          std::swap(ws_0,ws_1);

          if (method_==STRETCH)
            ws_0=nextStateStretch(CL_,mt_,amin_,amax_,betarun_,ws_0,ifeval_,acc_count);
          else
            ws_0=nextStateWalk(CL_,mt_,N_for_Walk_,rWalk_,betarun_,ws_0,ifeval_,acc_count);

        }
      acc_ratio_=1.0*acc_count/n_skip_/numWalkers_;


      std::cout<<"isamp="<<i<<"\t"<<betarun_<<"\t"<<ifeval_;
      std::cout<<"\t"<<runt<<"\t"<<timeIter<<"\t"<<acc_ratio_<<"\t";
      std::cout<<ws_0.logDataLikMean()<<"\t"<<ws_0.logDataLikStd();
      std::cout<<"\t"<<ws_0.logDataLikMax()<<"\t"<<ws_0.logDataLikMin();
      std::cout<<"\t"<<ws_0.logPriorLikMean()<<"\t"<<ws_0.logPriorLikStd();
      std::cout<<"\t"<<ws_0.logPriorLikMax()<<"\t"<<ws_0.logPriorLikMin()
              <<std::endl;
      ws_0.writeValuesTitles(os_par_);
      if (includePredictions_)
        ws_0.writeYValuesTitles(os_pred_,CL_);

      ws_0.writeValues(os_par_,i);
      if (includePredictions_)
        ws_0.writeYValues(os_pred_,i,CL_);

      os_mean_<<"isamp="<<i<<"\t"<<ifeval_<<"\t"<<runt<<"\t"<<timeIter<<"\t"<<acc_ratio_;
      ws_0.writeMeans(os_mean_);
      os_mean_.flush();
      os_par_.flush();
      if (includePredictions_)
        os_pred_.flush();
      ++i;
    }
}



std::ostream &mcmc::writeBody(std::ostream &s) const
{
  writeField(s,"Likelihood_Model",CL_->id());
  writeField(s,"Initial_Parameter",initial_);
  writeField(s,"Maximal_duration",maxDuration_minutes_);
  writeField(s,"Beta_Samples",betaSamples_);
  writeField(s,"Number_Samples",numSamples_);
  writeField(s,"Number_Walkers",numWalkers_);
  writeField(s,"Initial_Radius",radiusWalkers_);
  writeField(s,"eemcee_amin",amin_);
  writeField(s,"eemcee_amax",amax_);

  return s;
}

bool mcmc::readBody(std::string &line, std::istream &s)
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
      else if (!readField(line,s,"Beta_Samples",betaSamples_))
        return false;
      else if (!readField(line,s,"Number_Samples",numSamples_))
        return false;
      else if (!readField(line,s,"Number_Walkers",numWalkers_))
        return false;
      else if (!readField(line,s,"Initial_Radius",radiusWalkers_))
        return false;
      else if (!readField(line,s,"eemcee_amin",amin_))
        return false;
      else if (!readField(line,s,"eemcee_amax",amax_))
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
  return s;
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
