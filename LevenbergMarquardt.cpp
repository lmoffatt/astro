#include "LevenbergMarquardt.h"
#include "MatrixInverse.h"

#include <vector>
#include <limits>
#include <cmath>
#include <sstream>


LevenbergMarquardtMultinomial::LevenbergMarquardtMultinomial(ABC_Multinomial_Model *f,
                                                             const Parameters& initialParam, std::size_t numIterations):
  f_(f),
  w_(),
  ParamInitial_(initialParam),
  nPar_(ParamInitial_.size()),
  nData_(f->getData().numCells()),
  dx_(1e-9),
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
  smallGradient_(false){}



LevenbergMarquardtMultinomial::LevenbergMarquardtMultinomial(){}

LevenbergMarquardtMultinomial::LevenbergMarquardtMultinomial (const LevenbergMarquardtMultinomial& other):
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



LevenbergMarquardtMultinomial&
LevenbergMarquardtMultinomial::operator=(const LevenbergMarquardtMultinomial& other)
{
  if (this!=&other)
    {
      LevenbergMarquardtMultinomial tmp(other);
      swap(*this,tmp);
    }
  return *this;
}

void swap(LevenbergMarquardtMultinomial& one, LevenbergMarquardtMultinomial& other)
{
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

LevenbergMarquardtMultinomial& LevenbergMarquardtMultinomial::optimize()
{
  std::cout<<"num DF="<<f_->getData().numDF()<<"\n";
  std::cout<<"nIter"<<"\t"<<"PostLogLik"<<"\t"<<"LogPrior"<<"\t"<<"landa"<<"\t";
  std::cout<<"ParamChange"<<"\t"<<"PostLogLikChange"<<"\t"<<"NormGrad"<<"\n";
  initialize();
  while (!meetConvergenceCriteria())
    iterate();

  std::cout<<report();
  JTWJinv_=inv(JTWJ_);
  ParamCurr_.setCovariance(JTWJinv_);
  return *this;
}

double LevenbergMarquardtMultinomial::getEvidence()const
{

  return logPostLikCurr_-0.5*ParamCurr_.logDetCov();
}
double LevenbergMarquardtMultinomial::getLogPostLik()const
{

  return logPostLikCurr_;
}

double LevenbergMarquardtMultinomial::logDetPriorCov()const
{
  return f_->getPrior().logDetCov();
}
double LevenbergMarquardtMultinomial::logDetPostCov()const
{
  return ParamCurr_.logDetCov();
}

double LevenbergMarquardtMultinomial::logDetPostStd()const
{
  double d=0;
  for (std::size_t i=0; i<ParamCurr_.size(); ++i)
    d+=log(ParamCurr_.pStds()[i]);
  return 2*d;
}

void LevenbergMarquardtMultinomial::iterate()
{
  computeJacobian();
  computeSearchDirection();
  updateLanda();
  std::cout<<nIter_<<"\t"<<logPostLikCurr_<<"\t"<<logPriorCurr_<<"\t"<<landa_<<"\t";
  std::cout<<ParamChange_<<"\t"<<PostlogLikChange_<<"\t"<<GradNormPost_<<"\n";
  nIter_++;
}


void update_Jacobian(const ABC_Multinomial_Model* f,
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



void LevenbergMarquardtMultinomial::computeJacobian()
{

  update_Jacobian(f_,ParamCurr_,P_expCurr_,nPar_,J_,w_,epsilon_,prior_G_,nFeval_,dx_);

  update_Gradient(G_,prior_G_,epsilon_,J_);

  update_JTWJ(JTWJ_,JTWJ_landa_,ParamCurr_,J_,w_,nPar_,nData_);

}

void update_Likelihoods(const ABC_Multinomial_Model* f
                        ,const Parameters& ParamNew
                        ,std::vector<std::vector<double>>& P_exp_New
                        ,double& logLikNew,double& logPriorNew
                        , double& logPostLogLikNew
                        ,std::size_t& NFeval)
{
  P_exp_New=f->p_exp(ParamNew);
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
void LevenbergMarquardtMultinomial::initialize()
{
  nIter_=0;
  nFeval_=0;
  nDF_=f_->getData().numDF();
  ParamChange_=std::numeric_limits<double>::infinity();
  PostlogLikChange_=std::numeric_limits<double>::infinity();
  GradNormPost_=std::numeric_limits<double>::infinity();
  ParamCurr_=ParamInitial_;

  update_Likelihoods(f_,ParamCurr_,P_expCurr_,logLikCurr_,logPriorCurr_,logPostLikCurr_
                     ,nFeval_);
}


void LevenbergMarquardtMultinomial::computeSearchDirection()
{

  landaMult_JTWJ(JTWJ_landa_,JTWJ_,landa_);

  JTWJinv_=inv(JTWJ_landa_);


  update_Search_Direction(d_,JTWJinv_,G_,nPar_);


  ParamNew_=ParamCurr_;
  for (std::size_t i=0; i<nPar_; ++i)
    ParamNew_[i]-=d_[i];

  update_Likelihoods(f_,ParamNew_,P_exp_New_,logLikNew_,logPriorNew_,logPostLogLikNew_,nFeval_);
}






void LevenbergMarquardtMultinomial::updateLanda()
{
  std::size_t ifevalLoop=0;
  if ((logPostLogLikNew_<=logPostLikCurr_)
      ||(logPostLogLikNew_!=logPostLogLikNew_))
    {
      while(((logPostLogLikNew_<=logPostLikCurr_)
             &&(nFeval_<maxFeval_)
             )||(logPostLogLikNew_!=logPostLogLikNew_))
        {
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

    }
  else
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
      while((logPostLogLikNew_>logPostLikeNew0_)&&(!logPostLogLikNew_!=logPostLogLikNew_)
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

      if ((logPostLogLikNew_<=logPostLikeNew0_)||(logPostLogLikNew_!=logPostLogLikNew_))
        {
          landa_=landa0_;
          ParamNew_=ParamNew0_;
          logPostLogLikNew_=logPostLikeNew0_;
          logLikNew_=logLikNew0_;
          logPriorNew_=logPriorNew0_;
          P_exp_New_=P_exp_New0_;
        }
    }
  if (logPostLogLikNew_<=logPostLikCurr_)
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


bool LevenbergMarquardtMultinomial::meetConvergenceCriteria()
{
  surpassIter_=nIter_>=maxIter_;
  surpassFeval_=nFeval_>=maxFeval_;
  surpassLanda_=landa_>=maxLanda_;
  smallParamChange_=ParamChange_<ParamChangeMin_;
  smallSSChange_=(PostlogLikChange_)<PostLogLikChangeMin_;
  smallGradient_=GradNormPost_<GradientNormPostMin_;
  isNanLogPostLik_=logPostLikCurr_!=logPostLikCurr_;


  return surpassIter_||
      surpassFeval_||
      smallParamChange_||
      smallSSChange_||
      smallGradient_||
      surpassLanda_||
      isNanLogPostLik_;
}



Parameters LevenbergMarquardtMultinomial::OptimParameters()const
{
  return ParamCurr_;
}

std::size_t LevenbergMarquardtMultinomial::numEval()const
{
  return nFeval_;
}
std::size_t LevenbergMarquardtMultinomial::numIter()const
{
  return nIter_;
}
double LevenbergMarquardtMultinomial::LogLik()const
{
  return logLikCurr_;
}

double LevenbergMarquardtMultinomial::PostLogLik() const
{
  return logPostLikCurr_;
}
std::vector<double> LevenbergMarquardtMultinomial::Gradient()const
{
  return G_;
}
std::string LevenbergMarquardtMultinomial::report()const{
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


std::ostream& operator<<(std::ostream& s, LevenbergMarquardtMultinomial& LM)
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

std::vector<std::vector<double> >
ABC_Multinomial_Model::J(const Parameters &p,
                         const std::vector<std::vector<double> > &p_ex
                         , double delta) const
{


  std::vector<std::vector<double>> out(getData().numCells()
                                       ,std::vector<double>(p.size()));
  for (std::size_t i=0; i< p.size(); ++i)
    {
      unsigned ic=0;
      double dp=delta/getPrior().pStds()[i];
      Parameters pr(p);
      pr[i]=p[i]+dp;
      std::vector<std::vector<double>> p_exp_d=p_exp(pr);
      for (unsigned ii=0; ii<getData().ntot_obs().size(); ++ii)
        {
          for (unsigned j=0; j<getData().n_obs(ii).size(); ++j)
            {
              double dpdbeta=(log(p_exp_d[ii][j])-log(p_ex[ii][j]))/dp;
              out[ic][i]=dpdbeta;
              ++ic;
            }
        }
    }

  return out;


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

double ABC_Multinomial_Model::logLik(const Parameters &p) const
{
  return logLik(p_exp(p));
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

double ABC_Multinomial_Model::logPrior(const Parameters &p)const
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

std::vector<double> ABC_Multinomial_Model::PriorGradient(const Parameters &p)const
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
