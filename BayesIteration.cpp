#include <vector>
#include <fstream>
#include <cmath>
#include <map>
#include "BayesIteration.h"
#include "MatrixInverse.h"
#include "CortexLikelihood.h"
#include "CommandManager.h"
#include "CortexMeasure.h"

#include <tuple>

mcmcEnsambleState
nextStateStretch(const CortexLikelihood* CL,
                 std::mt19937_64& mt,
                 double a,
                 double beta,
                 mcmcEnsambleState& x
                 ,std::size_t& ifeval
                 ,std::vector<std::size_t>& accept_count)
{

  // mcmcWalkerState o(x.getMeanState(),x.numWalkers(),beta);
  std::size_t K=x.numWalkers();
  std::size_t N=x.numParameters();
  accept_count.resize(K,0);

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

  // first half
  for (std::size_t ih=0; ih<2; ++ih)
    {
#pragma omp parallel for default (none) shared (ih,K,CL,rs,a,x,N,beta,js,accept_count) reduction (+ : ifeval)
      for (std::size_t k=ih*K/2; k<K/2+ih*K/2; ++k)
        {
          std::size_t j=js[k];
          double z=std::pow(rs[k]*(std::sqrt(a)-std::sqrt(1.0/a))+std::sqrt(1/a),2);
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
              accept_count[k]++;
            }
        }

    }
  x.update_Mean();
  x.update_Covariance();
  return x;
}



mcmcEnsambleState
nextStateStretch_2(const CortexLikelihood* CL,
                   std::mt19937_64& mt,
                   double a,
                   double beta,
                   mcmcEnsambleState& x
                   ,std::size_t& ifeval
                   ,std::size_t& accept_count)
{

  mcmcEnsambleState o(x.getMeanState(),x.numWalkers(),beta);
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






mcmcEnsambleState
nextStateWalk(const CortexLikelihood* CL,
              std::mt19937_64& mt,
              std::size_t N,
              double rWalk,
              double beta,

              mcmcEnsambleState& x
              ,std::size_t& ifeval
              ,std::size_t& accept_count)
{
  mcmcEnsambleState o(x.getMeanState(),x.numWalkers(),beta);
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

mcmcEnsambleState::mcmcEnsambleState(const Parameters& prior, std::size_t numWalkers, double beta):
  mean_(prior)
,beta_(beta)
,f_(numWalkers)
,trMeans_(numWalkers)
,dataLiks_(numWalkers)
,priorLiks_(numWalkers)
{}

mcmcEnsambleState mcmcEnsambleState::create(const CortexLikelihood *f
                                            , const Parameters &initialParam
                                            , std::size_t numWalkers
                                            , double radiusWalkers
                                            , std::mt19937_64 &mt
                                            , std::size_t &ifeval,
                                            double beta)
{
  mcmcEnsambleState o(f->getPrior(),numWalkers,beta);
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

mcmcEnsambleState::mcmcEnsambleState(){}

const Parameters &mcmcEnsambleState::getMeanState() const
{
  return mean_;
}

std::size_t mcmcEnsambleState::numWalkers() const
{
  return trMeans_.size();
}

std::size_t mcmcEnsambleState::numParameters() const
{
  return mean_.size();
}

std::size_t mcmcEnsambleState::numMeasures() const
{
  return f_[0].size();
}


std::vector<double> &mcmcEnsambleState::operator[](std::size_t i)
{
  return trMeans_[i];
}

std::vector<std::vector<double>> &mcmcEnsambleState::f(std::size_t i)
{
  return f_[i];
}


const std::vector<double> &mcmcEnsambleState::operator[](std::size_t i) const
{
  return trMeans_[i];
}

double &mcmcEnsambleState::logDataLik(std::size_t i)
{
  return dataLiks_[i];
}

const double &mcmcEnsambleState::logDataLik(std::size_t i)const
{
  return dataLiks_[i];
}

double mcmcEnsambleState::logBetaLik(std::size_t i) const
{
  return beta_*logDataLik(i)+logPrior(i);
}
const double &mcmcEnsambleState::logDataLikMean() const
{
  return dataLikMean_;
}

const double &mcmcEnsambleState::logPriorLikMax() const
{
  return priorLikMax_;
}

const double &mcmcEnsambleState::logPriorLikMin() const
{
  return priorLikMin_;
}


const double &mcmcEnsambleState::logDataLikMax() const
{
  return dataLikMax_;
}

const double &mcmcEnsambleState::logDataLikMin() const
{
  return dataLikMin_;
}

double mcmcEnsambleState::beta()const
{
  return beta_;
}


const double &mcmcEnsambleState::logDataLikStd() const
{
  return dataLikStd_;
}
const double &mcmcEnsambleState::logPriorLikMean() const
{
  return priorLikMean_;
}
const double &mcmcEnsambleState::logPriorLikStd() const
{
  return priorLikStd_;
}



double &mcmcEnsambleState::logPrior(std::size_t i)
{
  return priorLiks_[i];
}

const double &mcmcEnsambleState::logPrior(std::size_t i) const
{
  return priorLiks_[i];
}

void mcmcEnsambleState::update_Mean()
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

void mcmcEnsambleState::update_Covariance()
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

std::ostream &mcmcEnsambleState::writeTrValues(std::ostream &s, std::size_t isample)
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

std::ostream &mcmcEnsambleState::writeValues(std::ostream &s, std::size_t isample) const
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

std::ostream &mcmcEnsambleState::writeYValues(std::ostream &s, std::size_t isample, const CortexLikelihood* CL) const
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
std::ostream &mcmcEnsambleState::writeYValuesTitles(std::ostream &s, const CortexLikelihood * CL) const
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



std::ostream &mcmcEnsambleState::writeValuesTitles(std::ostream &s) const
{
  s<<"isample"<<"\t"<<"beta"<<"\t"<<"i"<<"\t"<<"logDataLik(i)"<<"\t"<<"logPrior(i)";
  for (std::size_t k=0; k<numParameters(); ++k)
    s<<"\t"<<mean_.indexToName(k);
  s<<"\n";
  return s;
}

std::ostream &mcmcEnsambleState::writeMeans(std::ostream &s) const
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

std::ostream &mcmcEnsambleState::writeMeansTitles(std::ostream &s) const
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




std::ostream &mcmcEnsambleState::writeBody(std::ostream &s) const
{
  writeField(s,"mean_state",mean_);
  writeField(s,"beta",beta_);
  writeField(s,"transformed_Means",trMeans_);
  writeField(s,"dataLikelihood",dataLiks_);
  writeField(s,"priorLikelihood",priorLiks_);
  return s;
}

void mcmcEnsambleState::clear() {
  priorLiks_.clear();
  dataLiks_.clear();
  trMeans_.clear();
}

bool mcmcEnsambleState::readBody(std::string &line, std::istream &s)
{
  if (!readField(line,s,"mean_state",mean_)) return false;
  else if (!readField(line,s,"beta",beta_)) return false;
  else if (!readField(line,s,"transformed_Means",trMeans_)) return false;
  else if (!readField(line,s,"dataLikelihood",dataLiks_)) return false;
  else if (!readField(line,s,"priorLikelihood",priorLiks_)) return false;
  else return true;
}








//mcmc::mcmc(
//    const CortexLikelihood *f
//           , const Parameters &initialParam
//           , double maxDuration_minutes
//           , std::size_t quasiPriorSamples
//           , std::size_t betaSamples
//           , std::size_t eqSamples
//           , std::size_t numWalkers
//           , std::size_t nSkip
//           , double radiusWalkers
//           , const std::string &name
//           , double amin
//           , double amax
//           , double a_b
//           ,std::size_t N_for_Walk
//           ,double rWalk
//           , bool includePredictions
//           , std::mt19937_64::result_type initseed)
//  :
//    startTime_(std::chrono::steady_clock::now()),
//    mt_()
//  , fname_(name),
//    os_par_(),
//    os_mean_(),
//    os_pred_(),
//    CL_(f)
//  ,initial_(initialParam)
//  , maxDuration_minutes_(maxDuration_minutes)
//  ,quasiPriorSamples_(quasiPriorSamples)
//  , betaSamples_(betaSamples)
//  , numSamples_(quasiPriorSamples+betaSamples+eqSamples)
//  , numWalkers_(numWalkers)
//  ,n_skip_(nSkip)
//  ,acc_ratio_()
//  , radiusWalkers_(radiusWalkers)
//  ,amin_(amin)
//  ,amax_(amax)
//  ,a_b_(a_b)
//  ,N_for_Walk_(N_for_Walk)
//  ,rWalk_(rWalk)
//  ,includePredictions_(includePredictions)
//  ,method_(m)
//{

//  fname_=getSaveName(fname_);
//  os_par_.open(fname_+"_par.txt",std::ofstream::out);
//  os_mean_.open(fname_+"_mean.txt",std::ofstream::out);
//  if (includePredictions_)
//    os_pred_.open(fname_+"_pred.txt",std::ofstream::out);


//  if (initseed!=0)
//    {
//      mt_.seed(initseed);
//      std::cout<<"Seed for random generator provided="<<initseed<<std::endl;
//      if (os_mean_.is_open())
//        os_mean_<<"Seed for random generator provided="<<initseed<<std::endl;
//    }
//  else
//    {
//      std::random_device rd;
//      std::mt19937_64::result_type seed=rd();
//      mt_.seed(seed);
//      std::cout<<"Seed of random generator="<<seed<<std::endl;
//      if (os_mean_.is_open())
//        os_mean_<<"Seed of random generator="<<seed<<std::endl;
//    }

//}


void mcmc::initialize()
{

}


//bool mcmc::meetConditions()
//{
//  return (runt_<maxDuration_minutes_)&&(i_<numSamples_);
//}



void mcmc::iterate()
{

  prop_->next(this,CL_,betaManag_,ws_0_,mt_,&ifeval_);

  ++i_;
  auto tnow=std::chrono::steady_clock::now();
  auto d=tnow-startTime_;
  double t0=runt_;
  runt_=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
  timeIter_=60*(runt_-t0);

  sampler_->IteratePrint(this);


  //  std::size_t acc_count=0;
  //  betarun_=std::min(((1.0*std::max(i_,quasiPriorSamples_+1)-quasiPriorSamples_)/betaSamples_),1.0);


  //  for (std::size_t ii=0; ii<n_skip_; ++ii)
  //    {
  //      //          ws_1=nextState(CL_,mt_,a_,betarun_,ws_0,ifeval_,acc_count);
  //      //          std::swap(ws_0,ws_1);

  //      if (method_==STRETCH)
  //        ws_0_=nextStateStretch(CL_,mt_,amin_,amax_,a_b_,betarun_,ws_0_,ifeval_,acc_count);
  //      else
  //        ws_0_=nextStateWalk(CL_,mt_,N_for_Walk_,rWalk_,betarun_,ws_0_,ifeval_,acc_count);

  //    }
  //  acc_ratio_=1.0*acc_count/n_skip_/numWalkers_;

}

void mcmc::IterateTitlesPrint(std::ostream &s)const
{
  s<<"\t"<<"isamp="<<"\t"<<"betarun"<<"\t"<<"ifeval";
  s<<"\t"<<"runt"<<"\t"<<"timeIter";
}



//void mcmc::iteratePrint()
//{



//  write_titles(os_out_);

//  *os_out_<<"isamp="<<i_<<"\t"<<betarun_<<"\t"<<ifeval_;
//  *os_out_<<"\t"<<runt_<<"\t"<<timeIter<<"\t"<<acc_ratio_<<"\t";
//  *os_out_<<ws_0_.logDataLikMean()<<"\t"<<ws_0_.logDataLikStd();
//  *os_out_<<"\t"<<ws_0_.logDataLikMax()<<"\t"<<ws_0_.logDataLikMin();
//  *os_out_<<"\t"<<ws_0_.logPriorLikMean()<<"\t"<<ws_0_.logPriorLikStd();
//  *os_out_<<"\t"<<ws_0_.logPriorLikMax()<<"\t"<<ws_0_.logPriorLikMin()
//          <<std::endl;

//  ws_0_.writeValuesTitles(os_par_);
//  if (includePredictions_)
//    ws_0_.writeYValuesTitles(os_pred_,CL_);

//  ws_0_.writeValues(os_par_,i_);
//  if (includePredictions_)
//    ws_0_.writeYValues(os_pred_,i_,CL_);

//  os_mean_<<"isamp="<<i_<<"\t"<<ifeval_<<"\t"<<runt_<<"\t"<<timeIter<<"\t"<<acc_ratio_;
//  ws_0_.writeMeans(os_mean_);
//  os_mean_.flush();
//  os_par_.flush();
//  if (includePredictions_)
//    os_pred_.flush();

//}


//void mcmc::initializePrint()
//{
//  os_par_.precision(10);
//  if (includePredictions_)
//    os_pred_.precision(10);

//  *os_out_<<"mcmc run. method=";
//  initialPrint(*os_out_);



//  *os_out_<<"i"<<"\t"<<"beta"<<"\t"<<"ifeval"<<"\t"<<"runt"<<"\t"<<"timeIter";
//  *os_out_<<"\t"<<"acc_r"<<"\t"<<std::endl;
//   ws_0_->writeMeansTitles(os_out_);


//}



//void mcmc::initializePrint()
//{

//  BaseSampler
//  os_mean_.precision(10);

//  os_mean_<<"mcmc run. method=";
//  initialPrint(*os_mean_);
//  os_mean_<<"i"<<"\t"<<"ifeval"<<"\t"<<"runt"<<"\t"<<"timeIter"<<"\t"<<"acc_r";
//  ws_0_->writeMeansTitles(os_mean_);

//}


void mcmc::run()
{
  startTime_=std::chrono::steady_clock::now();
  ifeval_=0;
  i_=0;
  runt_=0;

  ws_0_=initializer_->initialize(this,CL_,initial_,mt_,&ifeval_,betaManag_);

  sampler_->InitializePrint(this);

  sampler_->IterateTitlesPrint(this,prop_,ws_0_,CL_);

  while (final_->meetConditions(this))
    {
      iterate();
    }

  finalize();

}



std::ostream &mcmc::writeBody(std::ostream &s) const
{
  writeField(s,"Likelihood_Model",CL_->id());

  writeField(s,"Initial_Parameter",initial_);

  writeField(s,"seedIsProvided",seedIsProvided_);

  writeField(s,"seed",seed_);

  writePtrField(s,"Initializer",initializer_);

  writePtrField(s,"proposedDistribution",prop_);

  writePtrField(s,"Finalizer",final_);

  writePtrField(s,"BetaCourse",betaManag_);

  writePtrField(s,"Sampler",sampler_);

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
      else if (!readField(line,s,"seedIsProvided",seedIsProvided_))
        return false;
      else if (!readField(line,s,"seed",seed_))
        return false;
      else if (!readPtrField(line,s,"Initializer",initializer_))
        return false;
      else if (!readPtrField(line,s,"proposedDistribution",prop_))
        return false;

      else if (!readPtrField(line,s,"Finalizer",final_))
        return false;
      else if (!readPtrField(line,s,"BetaCourse",betaManag_))
        return false;
      else if (!readPtrField(line,s,"Sampler",sampler_))
        return false;
      else return true;
    }
}



void mcmcWrun::push_back(mcmcEnsambleState w)
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


BaseMcmcState *Initialize_Ensamble::initialize(const mcmc *myMc, const CortexLikelihood *f, const Parameters &initialParam, std::mt19937_64 mt, std::size_t *ifeval, BaseBetaManagement *beta)
{

  mcmcEnsambleState*  o=new
      mcmcEnsambleState(f->getPrior(),numWalkers_,beta->beta(myMc));
  std::size_t i=0;
  while    (i<numWalkers_)
    {
      (*o)[i]=initialParam.randomSampleValues(mt,f->getPrior(),radius_);
      double logprior=f->logPrior(o->operator[](i));
      double logDataLik=f->logLik((*o)[i]);
      ++ifeval;
      o->logDataLik(i)=logDataLik;
      o->logPrior(i)=logprior;
      if (!std::isnan(logDataLik))
        ++i;

    }
  o->update_Mean();
  o->update_Covariance();
  return o;
}

double Beta_Ramp::beta(const mcmc *myMCMC) const
{
  return std::min(((1.0*std::max(myMCMC->iter(),quasiPriorSamples_+1)-quasiPriorSamples_)/betaSamples_),1.0);

}

void Sampler::InitializePrint(const mcmc *myMcmc)
{
  if (has_os_)
    myMcmc->write(*os_);
  myMcmc->write(os_mean_);
  myMcmc->write(os_par_);
  if (includePredictions_)
    {
      myMcmc->write(os_pred_);
    }

}

void Sampler::
IterateTitlesPrint(const mcmc *myMcmc, const BaseProposalDistribution *prop, const BaseMcmcState *state, const CortexLikelihood* CL)
{


  myMcmc->IterateTitlesPrint(*os_);
  *os_<<"\t";
  prop->TitlesPrint(*os_);
  *os_<<"\t";
  state->writeMeansTitles(*os_);
  *os_<<std::endl;

  myMcmc->IterateTitlesPrint(os_mean_);
  os_mean_<<"\t";
  prop->TitlesPrint(os_mean_);
  os_mean_<<"\t";
  state->writeMeansTitles(os_mean_);
  os_mean_<<std::endl;
  os_mean_.flush();

  myMcmc->IterateTitlesPrint(os_par_);
  os_par_<<"\t";
  prop->TitlesPrint(os_par_);
  os_par_<<"\t";
  state->writeValuesTitles(os_par_);
  os_par_<<std::endl;
  os_par_.flush();
  if (this->includePredictions_)
    {
      myMcmc->IterateTitlesPrint(os_pred_);
      os_pred_<<"\t";
      prop->TitlesPrint(os_pred_);
      os_pred_<<"\t";
      state->writeYValuesTitles(os_pred_,CL);
      os_pred_<<std::endl;
    }

}

void Sampler::IteratePrint(const mcmc *myMcmc)
{
  myMcmc->write(os_mean_);
  myMcmc->write(os_par_);
  if (includePredictions_)
    {
      myMcmc->write(os_pred_);
    }


}

BaseMcmcState *BaseMcmcState::createChild(const std::__cxx11::string &className)
{
  auto it=childs_.find(className);
  if (it!=childs_.end())
    {
      return it->second->create();
    }
  else
    return nullptr;

}

std::map<std::__cxx11::string, BaseMcmcState *> BaseMcmcState::getChilds()
{
  std::map<std::__cxx11::string, BaseMcmcState *> o;
  o[mcmcEnsambleState::ClassName()]=new mcmcEnsambleState;

  return o;
}
std::map<std::__cxx11::string, BaseMcmcState *> BaseMcmcState::childs_=BaseMcmcState::getChilds();

BaseBetaManagement *BaseBetaManagement::createChild(const std::__cxx11::string &className)
{
  auto it=childs_.find(className);
  if (it!=childs_.end())
    {
      return it->second->create();
    }
  else
    return nullptr;

}

std::map<std::__cxx11::string, BaseBetaManagement *> BaseBetaManagement::childs_=BaseBetaManagement::getChilds();

std::map<std::__cxx11::string, BaseBetaManagement *> BaseBetaManagement::getChilds()
{
  std::map<std::__cxx11::string, BaseBetaManagement *> o;
  o[Beta_Ramp::ClassName()]=new Beta_Ramp;
}

BaseInitializer *BaseInitializer::createChild(const std::__cxx11::string &className)
{
  auto it=childs_.find(className);
  if (it!=childs_.end())
    {
      return it->second->create();
    }
  else
    return nullptr;

}

std::map<std::__cxx11::string, BaseInitializer *> BaseInitializer::childs_=BaseInitializer::getChilds();

std::map<std::__cxx11::string, BaseInitializer *> BaseInitializer::getChilds()
{
  std::map<std::__cxx11::string, BaseInitializer *> o;
  o[Initialize_Ensamble::ClassName()]=new Initialize_Ensamble;
  return o;
}

BaseProposalDistribution *BaseProposalDistribution::createChild(const std::__cxx11::string &className)
{
  auto it=childs_.find(className);
  if (it!=childs_.end())
    {
      return it->second->create();
    }
  else
    return nullptr;

}

std::map<std::__cxx11::string, BaseProposalDistribution *> BaseProposalDistribution::getChilds()
{
  std::map<std::__cxx11::string, BaseProposalDistribution *>  o;
  o[StretchProposal::ClassName()]=new StretchProposal;
  return o;

}
std::map<std::__cxx11::string, BaseProposalDistribution *> BaseProposalDistribution::childs_=BaseProposalDistribution::getChilds();

BaseFinalize *BaseFinalize::createChild(const std::__cxx11::string &className)
{
  auto it=childs_.find(className);
  if (it!=childs_.end())
    {
      return it->second->create();
    }
  else
    return nullptr;

}

std::map<std::__cxx11::string, BaseFinalize *> BaseFinalize
::getChilds()
{
  std::map<std::__cxx11::string, BaseFinalize *>  o;
  o[Finalize::ClassName()]=new Finalize;
  return o;
}
std::map<std::__cxx11::string, BaseFinalize *> BaseFinalize::childs_=BaseFinalize::getChilds();

BaseSampler *BaseSampler::createChild(const std::__cxx11::string &className)
{
  auto it=childs_.find(className);
  if (it!=childs_.end())
    {
      return it->second->create();
    }
  else
    return nullptr;

}

std::map<std::__cxx11::string, BaseSampler *>
BaseSampler::getChilds()
{
  std::map<std::__cxx11::string, BaseSampler *> o;
  o[Sampler::ClassName()]=new Sampler;
  return o;
}
std::map<std::__cxx11::string, BaseSampler *>
BaseSampler::childs_=BaseSampler::getChilds();

bool Finalize::meetConditions(const mcmc *myMcmc) const
{
  return (myMcmc->iter()<numIter_)&&(myMcmc->runTime()<maxtime_);
}

void StretchProposal::next(const mcmc *myMc, const CortexLikelihood *CL, const BaseBetaManagement *beta, BaseMcmcState *s, std::mt19937_64 &mt, std::size_t *ifeval)
{


  if (s->myClass()==mcmcEnsambleState::ClassName())
    {
      delta_=std::min(delta_0_,std::pow(myMc->iter(),-2));
      if (acc_ratio_min_<acc_ratio_opt)
          a_/=1.0-delta_;
      else
        a_*=1+delta_;
      mcmcEnsambleState *e=reinterpret_cast<mcmcEnsambleState*> (s);
      double betarun=beta->beta(myMc);
      for (std::size_t ii=0; ii<n_skip_; ++ii)
        {
          *e=nextStateStretch(CL,mt,a_,betarun,*e,*ifeval,acc_count_);
        }
      acc_ratio_min_=1.0*min(acc_count_)/n_skip_;
      acc_ratio_max_=1.0*max(acc_count_)/n_skip_;
    }

}
