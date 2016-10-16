#include "StochasticLevenbergMarquardt.h"
#include "MatrixInverse.h"
#include "CortexLikelihood.h"

#include <vector>
#include <limits>
#include <cmath>
#include <sstream>


StochasticLevenbergMarquardtDistribution::StochasticLevenbergMarquardtDistribution
(const CortexLikelihood *f, double beta,
 const Parameters& initialParam
 , std::size_t numSamples, std::size_t n_skip_iter
 , double maxDuration_mins
 , const std::string& name,bool storeParameters,bool storeData):
  startTime_(std::chrono::steady_clock::now()),
  fname_(name),
  os_(),
  CL_(f),
  ParamInitial_(initialParam),
  nPar_(ParamInitial_.size()),
  nData_(f->getData().numCells())
, step_(nData_,nPar_),
  stoc_iter_(),
  step_report_(),
  storeIt_(n_skip_iter,storeParameters,storeData),
  mcmc_(numSamples,storeIt_),
  f_test_(maxDuration_mins),
  f_report_()

{
  fname_=getSaveName(fname_);
  os_.open(fname_,std::ofstream::out);
}






void update_Jacobian(const ABC_Distribution_Model* f,
                     const Parameters& ParamCurr,
                     const std::vector<std::vector<double>>&P_expCurr,
                     std::vector<std::vector<double>>& J,
                     bool& invalidJ,
                     std::vector<double>& w,
                     std::vector<double>& epsilon,
                     std::vector<double>& priorG,
                     double dx
                     )
{
  auto J0=f->J(ParamCurr, P_expCurr,dx);
  if (J0.empty())
    {
      invalidJ=true;
    }
  else
    {
      J=J0;
      invalidJ=false;
      w=f->weight(P_expCurr);
      priorG=f->PriorGradient(ParamCurr);
      epsilon=f->epsilon(P_expCurr);
    }
}
void update_Gradient(double beta,
                     std::vector<double>& beta_G,
                     const std::vector<double>& prior_G,
                     const std::vector<double>& epsilon,
                     const std::vector<std::vector<double>>& J)
{
  for (std::size_t i=0; i<beta_G.size(); ++i)
    {
      beta_G[i]=prior_G[i];
      for (std::size_t n=0; n<epsilon.size();++n)
        {
          beta_G[i]+=epsilon[n]*J[n][i]*beta;
        }
    }
}


void update_JTWJ_landa(const ABC_Distribution_Model* f,
                       double beta,
                       std::vector<std::vector<double>>& beta_JTWJ
                       ,std::vector<std::vector<double>>& beta_JTWJ_landa
                       ,const std::vector<std::vector<double>>& J
                       ,const std::vector<double>& w
                       ,std::size_t nPar,std::size_t nData)
{
  for (std::size_t i=0; i<nPar; ++i)
    for (std::size_t j=0; j<nPar; ++j)
      {
        beta_JTWJ[i][j]=f->getPrior().getInvCovariance()[i][j];
        for (std::size_t n=0; n<nData; ++n)
          {
            beta_JTWJ[i][j]+=J[n][i]*J[n][j]*w[n]*beta;
          }
        beta_JTWJ_landa[i][j]=beta_JTWJ[i][j];

      }
}




void update_Likelihoods(const ABC_Distribution_Model* f
                        ,double beta
                        ,const Parameters& ParamNew
                        ,std::vector<std::vector<double>>& P_exp_New
                        ,double& logLikNew,
                        double& logPriorNew,
                        double& beta_LogLik)
{
  P_exp_New=f->f(ParamNew);
  if (P_exp_New.empty())
    logLikNew=std::numeric_limits<double>::quiet_NaN();
  else
    logLikNew=f->logLik(P_exp_New);
  logPriorNew=f->logPrior(ParamNew);
  beta_LogLik=logPriorNew+beta*logLikNew;
}



inline
void landaMult_JTWJ(std::vector<std::vector<double>>& JTWJ_landa
                    ,const std::vector<std::vector<double>>& JTWJ
                    ,double landa)
{
  for (std::size_t i=0; i<JTWJ.size(); ++i)
    {
      JTWJ_landa[i][i]=JTWJ[i][i]*(1+landa);
    }
}



void update_Search_Direction(std::vector<double>& d
                             ,double& beta_LogLik_Next_exp
                             ,const std::vector<std::vector<double>>& JTWJinv
                             ,const std::vector<double>& G
                             ,double beta_LogLik
                             ,std::size_t nPar)
{
  for (std::size_t i=0; i<nPar; ++i)
    {
      d[i]=0;
      for (std::size_t j=0; j<nPar;++j)
        {
          d[i]+=JTWJinv[i][j]*G[j];
        }
    }

  beta_LogLik_Next_exp=beta_LogLik;
  for (std::size_t i=0; i<nPar; ++i)
    beta_LogLik_Next_exp+=d[i]*G[i];

}





















StochasticLevenbergMarquardtDistribution::state::state(std::size_t nData_, std::size_t nPar_)
  :
    w_(), logLik(),logPrior(),beta_(),beta_LogLik_(),beta_LogLik_Next_exp_(),
    landa_(), Param_(),ParamNext_(), P_exp_(),P_exp_Next_(),
    logLik_Next_(),logPrior_Next_(),beta_LogLik_Next_(), epsilon_beta_log_Lik_(),
    isJacobianInvalid_(false),
    J_(std::vector< std::vector< double> >(nData_,std::vector<double>(nPar_))),
    prior_G_(nPar_),
    beta_G_(nPar_)
  ,epsilon_(),
    beta_JTWJ_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),

    beta_JTWJ_landa_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
    beta_JTWJinv_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
    beta_d_(nPar_){
}

bool StochasticLevenbergMarquardtDistribution::state::update(const CortexLikelihood* CL_
                                                             ,double beta
                                                             ,Parameters &&param
                                                             , double dx
                                                             , double initLanda_
                                                             ,double multLanda
                                                             ,double maxLanda
                                                             ,double tol)
{
  Param_=param;
  beta_=beta;
  std::size_t nPar_=Param_.size();
  update_Likelihoods(CL_,beta_,Param_,P_exp_,logLik,logPrior,beta_LogLik_);
  update_Jacobian(CL_,Param_,P_exp_,J_
                  ,isJacobianInvalid_,w_,epsilon_,prior_G_,dx);
  std::size_t nData_=w_.size();
  if (!isJacobianInvalid_)
    {
      update_Gradient(beta_,beta_G_,prior_G_,epsilon_,J_);

      update_JTWJ_landa(CL_,beta_,beta_JTWJ_,beta_JTWJ_landa_,J_,w_,nPar_,nData_);

      double landa=initLanda_;
      while (!update_landa(CL_,landa,tol)&& (landa<maxLanda))
        landa*=multLanda;
      return landa<maxLanda;
    }
  else return false;

}

bool StochasticLevenbergMarquardtDistribution::state::update_landa
(const CortexLikelihood* CL_,double landa, double tol)
{
  landa_=landa;
  std::size_t nPar_=Param_.size();
  landaMult_JTWJ(beta_JTWJ_landa_,beta_JTWJ_,landa_);

  beta_JTWJinv_=inv(beta_JTWJ_landa_);


  update_Search_Direction(beta_d_,beta_LogLik_Next_exp_,beta_JTWJinv_,beta_G_,beta_LogLik_,nPar_);

  ParamNext_=Param_;
  for (std::size_t i=0; i<nPar_; ++i)
    ParamNext_[i]-=beta_d_[i];
  ParamNext_.setCovariance(beta_JTWJ_landa_);
  update_Likelihoods(CL_,beta_,ParamNext_,P_exp_Next_,logLik_Next_,logPrior_Next_,beta_LogLik_Next_);
  epsilon_beta_log_Lik_=beta_LogLik_Next_-beta_LogLik_Next_exp_;
  return std::abs(epsilon_beta_log_Lik_)<tol*std::sqrt(2*nPar_);


}




bool StochasticLevenbergMarquardtDistribution::MetropolisHastings_test::operator()(std::mt19937_64 &mt, const StochasticLevenbergMarquardtDistribution::state &current, const StochasticLevenbergMarquardtDistribution::state &candidate)
{
  logPcurrent=current.beta_LogLik_;
  logPcandidate=candidate.beta_LogLik_;

  logQforward=current.ParamNext_.logProb(candidate.Param_);

  logQbackward=candidate.ParamNext_.logProb(current.Param_);

  logA=logPcandidate-logQforward-logPcurrent+logQbackward;

  A=std::min(1.0,exp(logA));

  std::uniform_real_distribution<double> u(0,1);
  r=u(mt);
  if (candidate.isJacobianInvalid_)
    accept_=false;
  else
    accept_=r<A;
  return accept_;

}


void StochasticLevenbergMarquardtDistribution::stochastic_iter::operator()
(const CortexLikelihood* CL_,double beta,std::mt19937_64 &mt,
 StochasticLevenbergMarquardtDistribution::step &s)const
{
  s.nTryNext_=0;
  while (!s.s_Test_.update(CL_,beta,s.s_Curr_.ParamNext_.randomSample(mt),dx_,initLanda_,multLanda_,maxLanda_,tol_)&&
         s.nTryNext_<maxTries_)
    ++s.nTryNext_;
  if (s.nTryNext_==maxTries_)
    s.fail_=true;
  else
    {
      if (s.test_(mt,s.s_Curr_,s.s_Test_))
        {
          state tmp=std::move(s.s_Curr_);
          s.s_Curr_=std::move(s.s_Test_);
          s.s_Test_=std::move(tmp);
        }


    }
}
StochasticLevenbergMarquardtDistribution::Finalization_test::Finalization_test(double maxDur_in_min): maxDur_in_min_(maxDur_in_min){}

bool StochasticLevenbergMarquardtDistribution::Finalization_test::operator()(const mcmc& m,
    const StochasticLevenbergMarquardtDistribution::step &s, StochasticLevenbergMarquardtDistribution::Finalization_report &r) const
{
  if (s.fail_)
    {
      r.failStep_=true;
      return true;
    }
  else if (m.isFull())
    {
      r.full_mcmc_=true;
      return true;
    }
  else if (s.timeOpt_>=maxDur_in_min_)
    {
      r.surpassDuration_=true;
      return true;
    }
  else return false;
}



StochasticLevenbergMarquardtDistribution::mcmc::mcmc(std::size_t nsamples, const StochasticLevenbergMarquardtDistribution::state_store &ss)
  :nsamples_(0),logLiks(nsamples),parameters(),data(),sumLogLik(0),sumSqrLogLik(0)
{
  if (ss.storeParameters)
    parameters=std::vector<Parameters>(nsamples);
  if (ss.storeData)
    data=std::vector<std::vector<std::vector<double>>>(nsamples);
}
