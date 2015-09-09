#include "LevenbergMarquardt.h"
#include "MatrixInverse.h"
#include "CortexLikelihood.h"

#include <vector>
#include <limits>
#include <cmath>
#include <sstream>


LevenbergMarquardtDistribution::LevenbergMarquardtDistribution(const CortexLikelihood *f,
                                                               const Parameters& initialParam
                                                               , std::size_t numIterations
                                                               , const std::string& name):
  fname_(name),
  os_(),
  CL_(f),
  w_(),
  ParamInitial_(initialParam),
  nPar_(ParamInitial_.size()),
  nData_(f->getData().numCells()),
  dx_(std::sqrt(std::numeric_limits<double>::epsilon())),
  maxIter_(numIterations),
  maxFeval_(10000),
  ParamChangeMin_(1e-9),
  PostLogLikChangeMin_(1e-9),
  GradientNormPostMin_(0.0000001),
  maxLanda_(1e11),
  landa_(1000),
  v_(3),
  nIter_(0),
  nFeval_(0),
  nDF_(0),
  logPostLikCurr_(std::numeric_limits<double>::infinity()),
  logPostLikeNew0_(std::numeric_limits<double>::infinity()),
  ParamCurr_(ParamInitial_),
  ParamNew_(ParamCurr_),
  ParamNew0_(ParamCurr_),
  P_expCurr_(),
  P_exp_New_(),
  P_exp_New0_(),
  J_(std::vector< std::vector< double> >(nData_,std::vector<double>(nPar_))),
  G_(ParamCurr_.size()),
  JTWJ_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
  JTWJ_landa_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
  JTWJinv_(JTWJ_),
  d_(ParamCurr_.size()),
  surpassIter_(false),
  surpassFeval_(false),
  ParamChange_(std::numeric_limits<double>::infinity()),
  PostlogLikChange_(ParamChange_),
  GradNormPost_(ParamChange_),
  smallParamChange_(false),
  smallSSChange_(false),
  smallGradient_(false)
{
  fname_=getSaveName(fname_);
  os_.open(fname_,std::ofstream::out);


}



LevenbergMarquardtDistribution::LevenbergMarquardtDistribution(){}

LevenbergMarquardtDistribution::LevenbergMarquardtDistribution (const LevenbergMarquardtDistribution& other):
  os_(),
  fname_(other.fname_),
  w_(other.w_),
  ParamInitial_(other.ParamInitial_),
  nPar_(other.nPar_),
  nData_(other.nData_),
  dx_(other.dx_),
  maxIter_(other.maxIter_),
  maxFeval_(other.maxFeval_),
  ParamChangeMin_(other.ParamChangeMin_),
  PostLogLikChangeMin_(other.PostLogLikChangeMin_),
  GradientNormPostMin_(other.GradientNormPostMin_),
  maxLanda_(other.maxLanda_),
  landa_(other.landa_),
  landa0_(other.landa0_),
  v_(other.v_),
  nIter_(other.nIter_),
  nFeval_(other.nFeval_),
  nDF_(other.nDF_),


  logPostLikCurr_(other.logPostLikCurr_),
  logPostLikeNew0_ (other.logPostLikeNew0_),
  ParamCurr_ (other.ParamCurr_),
  ParamNew_ (other.ParamNew_),
  ParamNew0_ (other.ParamNew0_),
  P_expCurr_(other.P_expCurr_),
  P_exp_New_(other.P_exp_New_),
  P_exp_New0_(other.P_exp_New0_),

  J_(other.J_),
  G_(other.G_),
  JTWJ_(other.JTWJ_),
  JTWJinv_(other.JTWJinv_),

  d_(other.d_),


  surpassIter_(other.surpassIter_),
  surpassFeval_(other.surpassFeval_),

  ParamChange_(other.ParamChange_),
  PostlogLikChange_(other.PostlogLikChange_),
  GradNormPost_(other.GradNormPost_),

  smallParamChange_(other.smallParamChange_),

  smallSSChange_(other.smallSSChange_),

  smallGradient_(other.smallGradient_)

{}



LevenbergMarquardtDistribution&
LevenbergMarquardtDistribution::operator=(const LevenbergMarquardtDistribution& other)
{
  if (this!=&other)
    {
      LevenbergMarquardtDistribution tmp(other);
      swap(*this,tmp);
    }
  return *this;
}

void swap(LevenbergMarquardtDistribution& one, LevenbergMarquardtDistribution& other)
{
  std::swap(one.os_,other.os_);
  std::swap(one.w_,other.w_);

  std::swap(one.ParamInitial_,other.ParamInitial_);
  std::swap(one.nPar_,other.nPar_);
  std::swap(one.nData_,other.nData_);
  std::swap(one.dx_,other.dx_);
  std::swap(one.maxIter_,other.maxIter_);
  std::swap(one.maxFeval_,other.maxFeval_);
  std::swap(one.ParamChangeMin_,other.ParamChangeMin_);
  std::swap(one.PostLogLikChangeMin_,other.PostLogLikChangeMin_);
  std::swap(one.GradientNormPostMin_,other.GradientNormPostMin_);
  std::swap(one.maxLanda_,other.maxLanda_);
  std::swap(one.landa_,other.landa_);
  std::swap(one.landa0_,other.landa0_);
  std::swap(one.v_,other.v_);
  std::swap(one.nIter_,other.nIter_);
  std::swap(one.nFeval_,other.nFeval_);

  std::swap(one.nDF_,other.nDF_);


  std::swap(one.logPostLikCurr_,other.logPostLikCurr_);
  std::swap(one.logPostLikeNew0_,other.logPostLikeNew0_);
  std::swap(one.ParamCurr_,other.ParamCurr_);
  std::swap(one.ParamNew_,other.ParamNew_);
  std::swap(one.ParamNew0_,other.ParamNew0_);
  std::swap(one.P_expCurr_,other.P_expCurr_);
  std::swap(one.P_exp_New_,other.P_exp_New_);
  std::swap(one.P_exp_New_,other.P_exp_New_);
  std::swap(one.J_,other.J_);
  std::swap(one.G_,other.G_);
  std::swap(one.JTWJ_,other.JTWJ_);
  std::swap(one.JTWJinv_,other.JTWJinv_);

  std::swap(one.d_,other.d_);

  std::swap(one.ParamCurr_,other.ParamCurr_);

  std::swap(one.surpassIter_,other.surpassIter_);
  std::swap(one.surpassFeval_,other.surpassFeval_);

  std::swap(one.ParamChange_,other.ParamChange_);
  std::swap(one.PostlogLikChange_,other.PostlogLikChange_);
  std::swap(one.GradNormPost_,other.GradNormPost_);

  std::swap(one.smallParamChange_,other.smallParamChange_);

  std::swap(one.smallSSChange_,other.smallSSChange_);

  std::swap(one.smallGradient_,other.smallGradient_);

}

LevenbergMarquardtDistribution& LevenbergMarquardtDistribution::optimize()
{
  std::cout<<"num DF="<<CL_->getData().numDF()<<"\n";
  std::cout<<"nIter"<<"\t"<<"PostLogLik"<<"\t"<<"LogPrior"<<"\t"<<"landa"<<"\t";
  std::cout<<"ParamChange"<<"\t"<<"PostLogLikChange"<<"\t"<<"NormGrad"<<"\n";
  if (os_.is_open())
    {
      os_<<"num DF="<<CL_->getData().numDF()<<"\n";
      os_<<"nIter"<<"\t"<<"PostLogLik"<<"\t"<<"LogPrior"<<"\t"<<"landa"<<"\t";
      os_<<"ParamChange"<<"\t"<<"PostLogLikChange"<<"\t"<<"NormGrad"<<"\n";
    }
  initialize();
  while (!meetConvergenceCriteria())
    iterate();

  std::cout<<*this;
  if (os_.is_open())
    {
      os_<<*this;
      os_.flush();
    }
  if (!isNanLogPostLik_)
    {
      JTWJinv_=inv(JTWJ_);
      ParamCurr_.setCovariance(JTWJinv_);
    }
  return *this;
}

LevenbergMarquardtDistribution &LevenbergMarquardtDistribution::
optimize(std::string optname,
         double factor,
         std::size_t numSeeds,
         std::mt19937::result_type initseed)
{
  std::mt19937 mt;
  std::random_device rd;

  if (initseed!=0)
    {
      mt.seed(initseed);
      std::cout<<"Seed for random generator provided="<<initseed<<std::endl;
      numSeeds=1;
      if (os_.is_open())
        os_<<"Seed for random generator provided="<<initseed<<std::endl;
    }
  for (std::size_t i=0;i<numSeeds;i++)
    {
      if (initseed==0)
        {
          std::mt19937::result_type seed=rd();
          mt.seed(seed);
          std::cout<<"Seed of random generator="<<seed<<std::endl;
          if (os_.is_open())
            os_<<"Seed of random generator="<<seed<<std::endl;
        }


      ParamCurr_=CL_->getPrior().randomSample(mt,ParamInitial_,factor);
      optimize();
      std::string optfname=OptimParameters().save(optname+"_"+std::to_string(PostLogLik()));
      CortexMultinomialLikelihoodEvaluation CE(*CL_,OptimParameters());
      std::ofstream fo;
      std::string fnameout=optfname.substr(0,optfname.size()-3)+"_lik.txt";
      fo.open(fnameout.c_str());
      CE.extract(fo);
      fo.close();

    }
  return *this;
}

double LevenbergMarquardtDistribution::getEvidence()const
{

  return logPostLikCurr_-0.5*ParamCurr_.logDetCov();
}
double LevenbergMarquardtDistribution::getLogPostLik()const
{

  return logPostLikCurr_;
}

double LevenbergMarquardtDistribution::logDetPriorCov()const
{
  return CL_->getPrior().logDetCov();
}
double LevenbergMarquardtDistribution::logDetPostCov()const
{
  return ParamCurr_.logDetCov();
}

double LevenbergMarquardtDistribution::logDetPostStd()const
{
  double d=0;
  for (std::size_t i=0; i<ParamCurr_.size(); ++i)
    d+=log(ParamCurr_.pStds()[i]);
  return 2*d;
}

void LevenbergMarquardtDistribution::iterate()
{
  computeJacobian();
  computeSearchDirection();
  updateLanda();
  std::cout<<nIter_<<"\t"<<logPostLikCurr_<<"\t"<<logPriorCurr_<<"\t"<<landa_<<"\t";
  std::cout<<ParamChange_<<"\t"<<PostlogLikChange_<<"\t"<<GradNormPost_<<"\n";
  if (os_.is_open())
    {
      os_<<nIter_<<"\t"<<logPostLikCurr_<<"\t"<<logPriorCurr_<<"\t"<<landa_<<"\t";
      os_<<ParamChange_<<"\t"<<PostlogLikChange_<<"\t"<<GradNormPost_<<"\n";
      os_.flush();

    }
  nIter_++;
}


void update_Jacobian(const ABC_Distribution_Model* f,
                     const Parameters& ParamCurr,
                     const std::vector<std::vector<double>>&P_expCurr,std::size_t nPar,
                     std::vector<std::vector<double>>& J,
                     std::vector<double>& w,
                     std::vector<double>& epsilon,
                     std::vector<double>& priorG,
                     std::size_t& nFeval,
                     double dx
                     )
{
  J=f->J(ParamCurr, P_expCurr,dx);
  nFeval+=nPar;
  w=f->weight(P_expCurr);
  priorG=f->PriorGradient(ParamCurr);
  epsilon=f->epsilon(P_expCurr);
}

void update_Gradient(std::vector<double>& G,
                     const std::vector<double>& prior_G,
                     const std::vector<double>& epsilon,
                     const std::vector<std::vector<double>>& J)
{
  for (std::size_t i=0; i<G.size(); ++i)
    {
      G[i]=prior_G[i];
      for (std::size_t n=0; n<epsilon.size();++n)
        {
          G[i]+=epsilon[n]*J[n][i];
        }
    }
}


void update_JTWJ(std::vector<std::vector<double>>& JTWJ
                 ,std::vector<std::vector<double>>& JTWJ_landa
                 ,const Parameters& ParamCurr
                 ,const std::vector<std::vector<double>>& J
                 ,const std::vector<double>& w
                 ,std::size_t nPar,std::size_t nData)
{
  for (std::size_t i=0; i<nPar; ++i)
    for (std::size_t j=0; j<nPar; ++j)
      {
        JTWJ[i][j]=ParamCurr.getInvCovariance()[i][j];
        for (std::size_t n=0; n<nData; ++n)
          {
            JTWJ[i][j]+=J[n][i]*J[n][j]*w[n];
            JTWJ_landa[i][j]=JTWJ[i][j];
          }
      }
}



void LevenbergMarquardtDistribution::computeJacobian()
{

  update_Jacobian(CL_,ParamCurr_,P_expCurr_,nPar_,J_,w_,epsilon_,prior_G_,nFeval_,dx_);

  update_Gradient(G_,prior_G_,epsilon_,J_);

  update_JTWJ(JTWJ_,JTWJ_landa_,ParamCurr_,J_,w_,nPar_,nData_);

}

void update_Likelihoods(const ABC_Distribution_Model* f
                        ,const Parameters& ParamNew
                        ,std::vector<std::vector<double>>& P_exp_New
                        ,double& logLikNew,double& logPriorNew
                        , double& logPostLogLikNew
                        ,std::size_t& NFeval)
{
  P_exp_New=f->f(ParamNew);
  logLikNew=f->logLik(P_exp_New);
  logPriorNew=f->logPrior(ParamNew);
  logPostLogLikNew=logLikNew+logPriorNew;
  NFeval++;
}




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
                             ,const std::vector<std::vector<double>>& JTWJinv
                             ,const std::vector<double>& G
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

}
void LevenbergMarquardtDistribution::initialize()
{
  nIter_=0;
  nFeval_=0;
  nDF_=CL_->getData().numDF();
  ParamChange_=std::numeric_limits<double>::infinity();
  PostlogLikChange_=std::numeric_limits<double>::infinity();
  GradNormPost_=std::numeric_limits<double>::infinity();
  //ParamCurr_=ParamInitial_;

  update_Likelihoods(CL_,ParamCurr_,P_expCurr_,logLikCurr_,logPriorCurr_,logPostLikCurr_
                     ,nFeval_);
}


void LevenbergMarquardtDistribution::computeSearchDirection()
{

  landaMult_JTWJ(JTWJ_landa_,JTWJ_,landa_);

  JTWJinv_=inv(JTWJ_landa_);


  update_Search_Direction(d_,JTWJinv_,G_,nPar_);


  ParamNew_=ParamCurr_;
  for (std::size_t i=0; i<nPar_; ++i)
    ParamNew_[i]-=d_[i];

  update_Likelihoods(CL_,ParamNew_,P_exp_New_,logLikNew_,logPriorNew_,logPostLogLikNew_,nFeval_);
}






void LevenbergMarquardtDistribution::updateLanda()
{
  std::size_t ifevalLoop=0;

  /// no es mejor
  if ((logPostLogLikNew_<=logPostLikCurr_)
      ||std::isnan(logPostLogLikNew_))
    {
      /// mientras no sea mejor
      while(((logPostLogLikNew_<=logPostLikCurr_)
             &&(nFeval_<maxFeval_)
             )||isnan(logPostLogLikNew_))
        {
          /// si me paso freno...
          if (landa_*v_>=maxLanda_) break;
          landa0_=landa_;
          landa_=landa0_*v_;
          logLikNew0_=logLikNew_;
          logPriorNew0_=logPriorNew_;
          logPostLikeNew0_=logPostLogLikNew_;
          ParamNew0_=ParamNew_;
          P_exp_New0_=P_exp_New_;
          computeSearchDirection();
          ifevalLoop++;
          //   std::cerr<<landa_<<" ";
        }

      /// de aca sali porque cosegui algo mejor o porque el landa es muy grande
    }
  else
    {
      /// aqui ya tengo un candidato
      landa0_=landa_;
      landa_=landa_/v_;
      logLikNew0_=logLikNew_;
      logPriorNew0_=logPriorNew_;
      logPostLikeNew0_=logPostLogLikNew_;
      ParamNew0_=ParamNew_;
      P_exp_New0_=P_exp_New_;

      computeSearchDirection();
      ifevalLoop++;
      /// mientras pueda mejorar
      while((logPostLogLikNew_>logPostLikeNew0_)&&(!std::isnan(logPostLogLikNew_))
            &&(landa_>0.5))
        {
          landa0_=landa_;
          landa_=landa_/v_;
          logLikNew0_=logLikNew_;
          logPriorNew0_=logPriorNew_;
          logPostLikeNew0_=logPostLogLikNew_;
          ParamNew0_=ParamNew_;
          P_exp_New0_=P_exp_New_;
          computeSearchDirection();
          ifevalLoop++;
        }

      /// si me pase, recupero el anterior
      if ((logPostLogLikNew_<=logPostLikeNew0_)||std::isnan(logPostLogLikNew_))
        {
          landa_=landa0_;
          ParamNew_=ParamNew0_;
          logPostLogLikNew_=logPostLikeNew0_;
          logLikNew_=logLikNew0_;
          logPriorNew_=logPriorNew0_;
          P_exp_New_=P_exp_New0_;
        }
    }
  if ((logPostLogLikNew_<=logPostLikCurr_)|| std::isnan(logPostLogLikNew_))
    {
      ParamChange_=0;
      PostlogLikChange_=0;
    }
  else
    {
      ParamChange_=0;
      for (std::size_t i=0; i<nPar_; ++i)
        ParamChange_+=(ParamCurr_[i]-ParamNew_[i])*(ParamCurr_[i]-ParamNew_[i]);
      ParamChange_=sqrt(ParamChange_);
      PostlogLikChange_=-logPostLikCurr_+logPostLogLikNew_;
      ParamCurr_=ParamNew_;
      logPostLikCurr_=logPostLogLikNew_;
      logPriorCurr_=logPriorNew_;
      logLikCurr_=logLikNew_;
      P_expCurr_=P_exp_New_;
    }
  GradNormPost_=0;
  for (std::size_t i=0; i<nPar_; ++i)
    GradNormPost_+=G_[i]*G_[i];
  GradNormPost_=sqrt(GradNormPost_);
}


bool LevenbergMarquardtDistribution::meetConvergenceCriteria()
{
  surpassIter_=nIter_>=maxIter_;
  surpassFeval_=nFeval_>=maxFeval_;
  surpassLanda_=landa_>=maxLanda_;
  smallParamChange_=ParamChange_<ParamChangeMin_;
  smallSSChange_=(PostlogLikChange_)<PostLogLikChangeMin_;
  smallGradient_=GradNormPost_<GradientNormPostMin_;
  isNanLogPostLik_=std::isnan(logPostLikCurr_);


  return surpassIter_||
      surpassFeval_||
      smallParamChange_||
      smallSSChange_||
      smallGradient_||
      surpassLanda_||
      isNanLogPostLik_;
}



Parameters LevenbergMarquardtDistribution::OptimParameters()const
{
  return ParamCurr_;
}

std::size_t LevenbergMarquardtDistribution::numEval()const
{
  return nFeval_;
}
std::size_t LevenbergMarquardtDistribution::numIter()const
{
  return nIter_;
}
double LevenbergMarquardtDistribution::LogLik()const
{
  return logLikCurr_;
}

double LevenbergMarquardtDistribution::PostLogLik() const
{
  return logPostLikCurr_;
}
std::vector<double> LevenbergMarquardtDistribution::Gradient()const
{
  return G_;
}
std::string LevenbergMarquardtDistribution::report()const{
  std::stringstream output;

  output<<"Convergence critera: \t";
  if (surpassIter_)
    output<<nIter_<<" surpass number of iterations="<<maxIter_<<"\t";
  if (surpassFeval_)
    output<<nFeval_<<" surpass number of function evaluations="<<maxFeval_<<"\t";
  if (surpassLanda_)
    output<<landa_<<" surpass maximum landa value="<<maxLanda_<<"\t";
  if (smallParamChange_)
    output<<ParamChange_<<" surpass minium parameter change value="<<ParamChangeMin_<<"\t";
  if (smallSSChange_)
    output<<PostlogLikChange_<<" surpass minium SS change value="<<PostLogLikChangeMin_<<"\t";
  if (smallGradient_)
    output<<GradNormPost_<<" surpass minium gradient norm value="<<GradientNormPostMin_<<"\t";
  if (logPostLikCurr_!=logPostLikCurr_)
    output<<logPostLikCurr_<<" invalid value of the square sum"<<"\t";


  return output.str();
}


std::ostream& operator<<(std::ostream& s, LevenbergMarquardtDistribution& LM)
{
  s<<"LevenbergMarquardtParameters\n";
  s<<"Number of iterations \t"<<LM.numIter()<<"\t";
  s<<"Posterior Log Likelihood"<<LM.logPostLikCurr_<<"\n";
  s<<"Report \t"<<LM.report()<<"\n";
  s<<"Landa value \t"<<LM.landa_<<"\t";
  s<<"Parameter change \t"<<LM.ParamChangeMin_<<"\t";
  s<<"PostLogLik change \t"<<LM.PostLogLikChangeMin_<<"\t";
  s<<"Gradient norm \t"<<LM.GradNormPost_<<"\n";
  s<<"Evidence \t"<<LM.getEvidence()<<"\n";
  s<<"End\n";
  return s;
}




unsigned ABC_Freq_obs::numSamples() const
{
  return ntot_obs().size();
}

unsigned ABC_Freq_obs::numCells() const
{
  if (!ntot_obs().empty())
    return ntot_obs().size()*n_obs(0).size();
  else return 0;
}


double ABC_Distribution_Model::logPrior(const Parameters &p)const
{

  double sumL=0;
  if (!getPrior().hasCovariance())
    {
      for (unsigned i=0; i<p.size(); ++i)
        sumL+=std::pow(p[i]-getPrior()[i],2)/std::pow(getPrior().pStds()[i],2);
    }
  else
    {
      for (unsigned i=0; i<p.size(); ++i)
        for (unsigned j=0; j<p.size(); ++j)
          sumL+=(p[i]-getPrior()[i])*(p[j]-getPrior()[j])*getPrior().getInvCovariance()[i][j];
    }
  return 0.5*(-sumL);
  //+getPrior().logDetCov());
}

std::vector<double> ABC_Distribution_Model::PriorGradient(const Parameters &p)const
{
  std::vector<double> out(p.size(),0);
  if (!getPrior().hasCovariance())
    {
      for (unsigned i=0; i<p.size(); ++i)
        out[i]=(p[i]-getPrior()[i])/std::pow(getPrior().pStds()[i],2);
    }
  else
    {
      for (unsigned i=0; i<p.size(); ++i)
        for (unsigned j=0; j<p.size(); ++j)
          out[i]+=(p[j]-getPrior()[j])*getPrior().getInvCovariance()[i][j];
    }
  return out;
}





std::vector<std::vector<double> >
ABC_Distribution_Model::J(const Parameters &p,
                          const std::vector<std::vector<double> > &f0
                          , double delta) const
{


  std::vector<std::vector<double>> out(getData().numCells()
                                       ,std::vector<double>(p.size()));
  for (std::size_t i=0; i< p.size(); ++i)
    {
      unsigned ic=0;
      double dp=delta*(1/getPrior().pStds()[i]+p[i]);
      Parameters pr(p);
      pr[i]=p[i]+dp;
      std::vector<std::vector<double>> f_delta=f(pr);
      for (unsigned ii=0; ii<getData().ntot_obs().size(); ++ii)
        {
          for (unsigned j=0; j<getData().n_obs(ii).size(); ++j)
            {
              if (f0[ii][j]>0)
                {
                  double dpdbeta=(log(f_delta[ii][j])-log(f0[ii][j]))/dp;
                  out[ic][i]=dpdbeta;
                }
              else
                out[ic][i]=0;
              ++ic;
            }
        }
    }

  return out;


}

double ABC_Distribution_Model::logLik(const Parameters &p) const
{
  return logLik(f(p));
}



std::vector<std::vector<double> >
ABC_Multinomial_Model::logLikCells(const std::vector<std::vector<double> > &p) const
{
  std::vector<std::vector<double> > cL(p.size(),std::vector<double>(p[0].size()));
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      double n=getData().ntot_obs()[i];
      double pL=0;
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double nj=getData().n_obs(i)[j];
          if (nj>0)
            cL[i][j]=nj*log(n*p[i][j]/nj);
          else
            cL[i][j]=0;
        }
    }
  return cL;

}

std::vector<double> ABC_Multinomial_Model::logLikSamples(const std::vector<std::vector<double> > &p) const
{
  std::vector<double>  cL(p.size(),0);
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      double n=getData().ntot_obs()[i];
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double nj=getData().n_obs(i)[j];
          if (nj>0)
            cL[i]+=nj*log(n*p[i][j]/nj);
        }
    }
  return cL;
}

double ABC_Multinomial_Model::logLik(const std::vector<std::vector<double> > &p) const
{
  double sumL=0;
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      double n=getData().ntot_obs()[i];
      double pL=0;
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double nj=getData().n_obs(i)[j];
          if (nj>0)
            pL+=nj*log(n*p[i][j]/nj);
        }
      sumL+=pL;
    }
  return sumL;
}


std::vector<double> ABC_Multinomial_Model::epsilon(const std::vector<std::vector<double> > &p) const
{
  std::vector<double> out(getData().numCells());
  unsigned ic=0;
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      double n=getData().ntot_obs()[i];
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double nj=getData().n_obs(i)[j];
          out[ic]=n*p[i][j]-nj;
          ++ic;
        }
    }
  return out;

}

const std::vector<double> ABC_Multinomial_Model::weight(const std::vector<std::vector<double> > &p) const
{
  std::vector<double> out(getData().numCells());
  unsigned ic=0;
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      double n=getData().ntot_obs()[i];
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          out[ic]=n*p[i][j];
          ++ic;
        }
    }
  return out;

}


std::vector<std::vector<double> >
ABC_MultiPoison_Model::logLikCells(const std::vector<std::vector<double> > &landa) const
{
  std::vector<std::vector<double> > cL(landa.size(),std::vector<double>(landa[0].size(),0));
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double k=getData().n_obs()[i][j];
          double landav=landa[i][j];
          if (k>0)
            cL[i][j]=+k*log(landav)-landav-lgamma(k+1);
        }
    }
  return cL;

}

std::vector<double> ABC_MultiPoison_Model::logLikSamples(const std::vector<std::vector<double> > &landa) const
{
  std::vector<double>  cL(landa.size(),0);
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double k=getData().n_obs()[i][j];
          double landav=landa[i][j];
          if (landav>0)
            cL[i]+=k*log(landav)-landav-lgamma(k+1);
        }
    }
  return cL;


}


double ABC_MultiPoison_Model::logLik(const std::vector<std::vector<double> > &landa) const
{
  double sumL=0;
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double k=getData().n_obs(i)[j];
          double landav=landa[i][j];
          if (std::isnan(landav))
            return landav;
          else if (landav!=0)
            {
              double logL=k*log(landav)-landav-lgamma(k+1);

              sumL+=logL;
            }
        }
    }
  return sumL;
}


double ABC_MultiPoison_Model::logLik(const Parameters &p) const
{
  return logLik(f(p));
}

std::vector<double> ABC_MultiPoison_Model::epsilon(const std::vector<std::vector<double> > &landa) const
{
  std::vector<double> out(getData().numCells());
  unsigned ic=0;
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      double n=getData().ntot_obs()[i];
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          double k=getData().n_obs(i)[j];
          double landav=landa[i][j];
          out[ic]=landav-k;
          ++ic;
        }
    }
  return out;

}

const std::vector<double> ABC_MultiPoison_Model::weight(const std::vector<std::vector<double> > &landa) const
{
  std::vector<double> out(getData().numCells());
  unsigned ic=0;
  for (unsigned i=0; i<getData().ntot_obs().size(); ++i)
    {
      for (unsigned j=0; j<getData().n_obs(i).size(); ++j)
        {
          out[ic]=landa[i][j];
          ++ic;
        }
    }
  return out;

}

