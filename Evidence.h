#ifndef EVIDENCE_H
#define EVIDENCE_H

#include "Matrix.h"
#include "Distributions.h"

#include <random>

#include <cmath>



struct D_logL
{
  M_Matrix<double> G;
  M_Matrix<double> H;
};


struct mcmc_prior
{
  M_Matrix<double> param;
  double logPrior;
  D_logL D_prior;
};


struct mcmc_post: public mcmc_prior
{
  mcmc_post(mcmc_prior &&p): mcmc_prior(p), isValid(false), f(),logLik(0){}
  mcmc_post(){}
  bool isValid;
  M_Matrix<double> f;
  double logLik;
  double logbPL(double beta)const {return logPrior+logLik*beta;}
};


struct mcmc_Dpost: public mcmc_post
{
  mcmc_Dpost(mcmc_post &&p): mcmc_post(p), D_lik(){}

  D_logL D_lik;
};


template<typename Dist=MultivariateGaussian>
struct mcmc_step: public mcmc_Dpost
{
  mcmc_step(mcmc_Dpost &&p, double beta_): mcmc_Dpost(p),beta(beta_), proposed(){}
  double beta;
 double logbPL()const {return logbPL(beta);}
  Dist proposed;

};



template <class T>
// T regular type
class SamplesSeries
{
public:
  SamplesSeries(std::size_t n):
    n_(0),samples_(n){}

  bool push_back(T x)
  {
    if (n_<samples_.size())
      {
        samples_[n_]=x;
        ++n_;
        return true;
      }
    else return false;
  }

  bool full()const
  {
    return n_==samples_.size();
  }

  T sample(std::mt19937_64& mt, std::size_t i0=0)const
  {
    std::uniform_int_distribution<std::size_t> u(i0,n_-1);
    return samples_[u(mt)];
  }

  template<class F>
  // F(T)->double
  double mean(const F& f, std::size_t i0=0)
  {
    double s=f(samples_[i0]);
    for (std::size_t i=i0+1; i<n_; ++i)
      s+=f(samples_[i]);
    return s/(n_-i0);
  }

  template<class R, class F>
  SamplesSeries<R>
  getSample(const F& f, std::size_t i0=0)
  {
    SamplesSeries<R> o(samples_.size()-i0);
    for (std::size_t i=i0; i<samples_.size(); ++i)
      {
        o.push_back(f(samples_[i]));
      }
    return o;
  }

  std::size_t size()const {return n_;}


  template<class F>
  std::pair<double,double> mean_var(const F& f, std::size_t i0=0)
  {
    double m=mean(f,i0);
    double s=f(samples_[i0]);
    for (std::size_t i=i0+1; i<n_; ++i)
      s+=std::pow(f(samples_[i]),2);
    double var=s/(n_-i0)-std::pow(m,2);
    return {m,var};
  }


private:
  std::size_t n_;
  std::vector<T> samples_;
};



template <>
class SamplesSeries<double>
{
public:
  SamplesSeries(std::size_t n):
    n_(0),samples_(n),u_(0,n-1){}

  bool push_back(double && x)
  {
    if (n_<samples_.size())
      {
        samples_[n_]=x;
        ++n_;
        return true;
      }
    else return false;
  }

  bool full()const
  {
    return n_==samples_.size();
  }
  std::size_t size()const {return n_;}

  double sample(std::size_t i0,std::mt19937_64& mt)
  {
    std::uniform_int_distribution<std::size_t> u(i0,n_-1);
    return samples_[u(mt)];
  }

  template<class F>
  double mean(std::size_t i0,const F& f)
  {
    double s=f(samples_[i0]);
    for (std::size_t i=i0+1; i<n_; ++i)
      s+=f(samples_[i]);
    return s/n_;
  }

  template<class R, class F>
  SamplesSeries<R>
  getSample(std::size_t i0,const F& f)
  {
    SamplesSeries<R> o(samples_.size());
    for (std::size_t i=i0; i<samples_.size(); ++i)
      {
        o.push_back(f(samples_[i]));
      }
    return o;
  }


  std::size_t size(){return n_;}

  template<class F>
  std::pair<double,double> mean_var(std::size_t i0,const F& f)
  {
    double m=mean(i0,f);
    double s=f(samples_[i0]);
    for (std::size_t i=i0+1; i<n_; ++i)
      s+=std::pow(f(samples_[i]),2);
    double var=s/n_-std::pow(m,2);
    return {m,var};
  }


private:
  std::vector<double> samples_;
  std::size_t n_;
  std::uniform_int_distribution<std::size_t> u_;
};






template<class D, template<class> class M>
class Poisson_Likelihood
{
public:
  double logLikelihood(double landa,std::size_t k)
  {
    return k*log(landa)-landa-lgamma(k+1);
  }


  double logL(const D& data,const M_Matrix<double>& landa)
  {
    M_Matrix<std::size_t> k=data();
    double sumLogL=0;
    for (std::size_t i=0; i<k.nrows(); ++i)
      {
        for (std::size_t j=0; j<k.ncols(); ++j)
          sumLogL+=logLikelihood(landa(i,j),k(i,j));
      }
    return sumLogL;
  }

  mcmc_post get_mcmc_Post(const M<D>& model, const D& data, M_Matrix<double> param)
  {
    mcmc_prior p=model.prior(data,param);
    mcmc_post out(std::move(p));
    out.f=model.f(data,param);
    out.logLik=logL(data,out.f);
    if (std::isnan(out.logLik))
      out.isValid=false;
    else out.isValid=true;
    return out;
  }

};

template<class D, template<class> class M>
class Poisson_DLikelihood: public Poisson_Likelihood<D,M>
{
public:
  mcmc_Dpost get_mcmc_Dpost(const M<D>& model, const D& data, const M_Matrix<double>& param,double beta)
  {
    return mcmcD(model,data,param,get_mcmc_post(model,data,param,beta));
  }


  mcmc_Dpost get_mcmc_Dpost(const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post p)
  {
    mcmc_Dpost out(std::move(p));
    M_Matrix<std::size_t> k=data();
    M_Matrix<double> logLanda_0=out.f.apply([](double x){return log10(x);});
    M_Matrix<double> J=get_J(model,  data, param,logLanda_0  );
    out.D_lik.G=get_G(out.f,k,J);
    out.D_lik.H=get_H(out.f,J);
    return out;
  }

private:
  M_Matrix<double> get_G(const M_Matrix<double>& landa, const M_Matrix<std::size_t>& k, const M_Matrix<double>& J)const
  {
    M_Matrix<double> epsilon=landa-M_Matrix<double>(k);
    return epsilon*J;
  }

  M_Matrix<double> get_H(const M_Matrix<double>& landa, const M_Matrix<double>& J)const
  {
    std::size_t n=landa.size();
    std::size_t npar=J.ncols();
    M_Matrix<double> out(npar,npar,0.0);
    for (std::size_t j=0; j<npar; ++j)
      for (std::size_t j2=j; j2<npar; ++j2)
        for (std::size_t i=0; i<n; ++i)
          out(j,j2)=landa[i]*J(i,j)*J(i,j2);
    for (std::size_t j=0; j<npar; ++j)
      for (std::size_t j2=0; j2<j; ++j2)
        out(j,j2)=out(j2,j);

    return out;

  }

  M_Matrix<double>
  get_J(const M<D>& model, const D& data, const M_Matrix<double>& param,
        const M_Matrix<double>& logLanda_0  ,
        double delta=1e-4)const
  {
    M_Matrix<double> k=data();
    std::size_t n=k.size();
    std::size_t npar=param.size();
    M_Matrix<double> out(n,npar,0.0);
    for (std::size_t j=0; j<npar; ++j)
      {
        M_Matrix<double> p(param);
        p[j]+=delta;
        M_Matrix<double> logLanda_i=model.logLanda(data,p);
        for (std::size_t i=0; i<n; ++i)
          out(i,j)=(logLanda_i[i]-logLanda_0[i])/delta;

      }
    return out;


  }

};




class D_Likelihood
{
public:
  double logL()const{ return logL_;}
  M_Matrix<double> G()const{ return G_;}
  M_Matrix<double> H()const {return H_;}

  D_Likelihood(double _logL,M_Matrix<double> _G, M_Matrix<double> _H)
    :logL_(_logL),G_(_G),H_(_H){}

private:
  double logL_;
  M_Matrix<double> G_;
  M_Matrix<double> H_;


};



class Likelihood_Evaluation
{
public:
  double logLikelihood() const {return logL_;}
  double logPrior()const {return logPrior_;}
  double beta_logLikelihood()const {return beta_logL_;}
  double logProb(const M_Matrix<double>& s)const {
    return g_.logP(s);
  }
  M_Matrix<double> sample(std::mt19937_64& mt)
  {
    auto o= g_.sample(mt);
    return M_Matrix<double>(1,o.size(),o);
  }

  const M_Matrix<double>& sample()const {return parameters_;}

  const M_Matrix<double>&  G() const {return G_;}
  const M_Matrix<double>&  H() const {return H_;}
  const M_Matrix<double>&  beta_G() const {return beta_G_;}
  const M_Matrix<double>&  beta_H() const {return beta_H_;}
  bool isValid()const{return isValid_;}


  M_Matrix<double> parameters_;
  std::size_t k;
  double beta;
  bool isValid_;
  double logL_;
  double logPrior_;
  double beta_logL_;
  M_Matrix<double>  G_;
  M_Matrix<double>  H_;
  M_Matrix<double>  beta_G_;
  M_Matrix<double>  beta_H_;

  double landa_;
  M_Matrix<double> landa_beta_H_;
  M_Matrix<double> landa_beta_Hinv_;
  M_Matrix<double> d_;
  M_Matrix<double> next_parameters_;

  double exp_next_beta_logL_;


  MultivariateGaussian g_;


  void update_landa(double landa)
  {
    landa_=landa;
    for (std::size_t i=0; i<k; ++i)
      landa_beta_H_(i,i)=beta_H_(i,i)*(1+landa);
    landa_beta_Hinv_=invSafe(landa_beta_H_);
    d_=landa_beta_Hinv_*beta_G_;
    next_parameters_=parameters_+d_;
    exp_next_beta_logL_=0.5*TranspMult(d_,beta_G_)[0];
  }
};



struct ModelPredictions
{
  struct gaussian
  {
    M_Matrix<double> epsilon_;
    M_Matrix<double> variance_;
    M_Matrix<double> logVariance_;
    M_Matrix<double> weight_;

  };



  bool isValid_;
  gaussian g_;
  double logPrior_;
  double logLikelihood_;
  M_Matrix<double> priorG_;
  M_Matrix<double> priorH_;


};






template<class Model, class Data>
struct LevenbergMarquardt_gaussian_aprox
{
  const Model& model_;
  const Data& data_;
  double dp_;

  LevenbergMarquardt_gaussian_aprox(const Model& model,
                                    const Data& data,
                                    double dp=1E-5
      ):
    model_(model),data_(data),dp_(dp){}



  double get_betaLogLikelihood(const M_Matrix<double>& parameters, double beta_)const
  {
    ModelPredictions pred=model_(parameters,data_);
    return pred.logLikelihood_*beta_+pred.logPrior_;

  }



  Likelihood_Evaluation get_Hessian(const M_Matrix<double>& parameters, double beta_)const
  {
    ModelPredictions pred=model_(parameters,data_);
    Likelihood_Evaluation o;
    std::size_t n=pred.g_.epsilon_.size();

    o.logL_=-n*0.5*log(2*PI);
    for (std::size_t j=0; j<n; ++j)
      {
        o.logL_-=0.5*pred.g_.logVariance_[j];
        o.logL_-=0.5*std::pow(pred.g_.epsilon_[j],2)/pred.g_.variance_[j];

      }



    o.beta_logL_=o.logL_*beta_+pred.logPrior_;



    std::size_t k=parameters.size();
    M_Matrix<double>  J_e(n,k);
    M_Matrix<double>  J_logv(n,k);


    M_Matrix<double> p=parameters;
    for (std::size_t i=0; i<k; ++i)
      {

        double po=p[i];
        p[i]+=dp_;

        ModelPredictions pred_i=model_(p,data_);
        for (std::size_t j=0; j<n; ++j)
          {
            J_e(j,i)=(pred_i.g_.epsilon_[j]-pred.g_.epsilon_[j])/dp_;
            J_logv(j,i)=(pred_i.g_.logVariance_[j]-pred.g_.logVariance_[j])/dp_;
          }
        p[i]=po;
      }

    o.G_=M_Matrix<double>(k,1);
    o.H_=M_Matrix<double>(k,k);
    o.beta_G_=M_Matrix<double>(k,1);
    o.beta_H_=M_Matrix<double>(k,k);
    for (std::size_t i=0; i<k; ++i)
      {
        o.G_(i,0)=0;

        for (std::size_t j=0; j<n; ++j)
          {
            o.G_(i,0)+=(std::pow(pred.g_.epsilon_[j],2)/pred.g_.variance_[j]-1.0)/2.0*J_logv(j,i);
            o.G_(i,0)+=pred.g_.epsilon_[j]/pred.g_.variance_[j]*J_e(j,i);
          }
        o.beta_G_(i,0)=pred.priorG_(i,0)+o.G_(i,0)*beta_;

        for (std::size_t i2=0; i2<k; ++i2)
          {
            o.H_(i,i2)=0;
            for (std::size_t j=0; j<n; ++j)
              {
                o.H_(i,i2)+=1.0/2.0*J_logv(j,i)*J_logv(j,i2);
                o.H_(i,i2)*=1/pred.g_.variance_[j]*J_e(j,i)*J_e(j,i2);
              }
            o.beta_H_(i,i2)=pred.priorH_(i,i2)+beta_*o.H_(i,i2);
          }
      }
    return o;



  }


};



template<
    class D, template<class> class M
    , template<class,template<class> class > class D_Lik=Poisson_DLikelihood
    >
class LevenbergMarquardt_step
{
public:

  struct trust_region
  {
    double r_;

    bool operator()(double expected, double found, std::size_t k)const
    {
      if (std::pow(expected-found,2)<2*k*r_)
        return true;
      else
        return false;
    }
  };

  trust_region t_;
  double landa0_;
  double v_;
  std::size_t maxLoop_;

public:

  LevenbergMarquardt_step( double  r_trust_region,
                           double landa0=1E8,
                           double v=3,
                           std::size_t maxLoop=10)
    :t_{r_trust_region},landa0_{landa0},v_(v),maxLoop_(maxLoop)
  {
  }
  struct LM_logL:public D_logL
  {
    double logL;
    M_Matrix<double> Hinv;
    M_Matrix<double> d;
    double exp_next_logL;
  };

  static LM_logL update_landa(const mcmc_Dpost& postL,double landa,double beta)
  {
    LM_logL out;
    out.G=postL.D_prior.G+postL.D_lik.G*beta;
    out.H=postL.D_prior.H+postL.D_lik.H*beta;
    for (std::size_t i=0; i<out.H.nrows(); ++i)
      out.H(i,i)=out.H(i,i)*(1+landa);
    out.Hinv=invSafe(out.H);
    out.d=out.Hinv*out.G;
    out.exp_next_logL=0.5*TranspMult(out.d,out.G)[0];
    return out;
  }

  mcmc_step<> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, double beta)
  {
     return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta);
  }

  mcmc_step<> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,double beta)
  {
    mcmc_step<> out(std::move(p),beta);
    double landa=landa0_;
    std::size_t iloop=0;
    LM_logL LM_logL=update_landa( out.D_lik,landa);

    M_Matrix<double> next=p.param+LM_logL.d;
    mcmc_post cand=L.get_mcmc_Post(model,data,next);
    while (!t_(LM_logL.exp_next_logL,cand.logbPL(beta),p.param.size())&&iloop<maxLoop_)
      {
        landa*=v_;
        LM_logL=update_landa( out.D_lik,landa);
        next=p.param+LM_logL.d;
        cand=L.get_mcmc_Post(model,data,next);
        ++iloop;
      }
    out.proposed=MultivariateGaussian(next,LM_logL.H,LM_logL.Hinv);
    return out;
  }
};



template
<
    class D, template<class> class M
    , template<class,template<class> class >
    class D_Lik=Poisson_DLikelihood
    ,template<class, template<class> class, template<class,template<class> class > class>
    class propDist=LevenbergMarquardt_step
    >
class Metropolis_Hastings_mcmc
{
  struct test
  {
    test()=default;

    bool operator()(const mcmc_step<>& sLik
                    ,const mcmc_step<>& cLik
                    ,std::mt19937_64 &mt)
    {
      double logPcurrent=sLik.logbPL();
      double logPcandidate=cLik.logbPL();

      double logQforward=sLik.proposed.logP(sLik.param);
      double logQbackward=cLik.proposed.logP(cLik.param);

      double logA=logPcandidate-logQforward-logPcurrent+logQbackward;
      double A=std::min(1.0,exp(logA));

      std::uniform_real_distribution<double> u(0,1);
      double r=u(mt);
      bool accept_=r<A;
      return accept_;
      
    }



  };

public:

  SamplesSeries<mcmc_step<>> run
  (const propDist<D,M,D_Lik>& LMLik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<>& sDist,
   std::size_t nsamples, std::size_t nskip,std::mt19937_64& mt)const
  {
    SamplesSeries<mcmc_step<>> o(nsamples);


    if (!sDist.isValid)
      return o;
    std::size_t i=0;

    M_Matrix<double> c=sDist.proposed.sample(mt);

    mcmc_step<> cDist=LMLik.get_mcmc_step(lik,model,data,c,sDist.beta);

    test t;

    while(cDist.isValid)
      {
        if (t(sDist,cDist,mt))
          {
            sDist=std::move(cDist);
          }
        if (i%nskip==0)
          o.push_back(sDist);
        if (o.full())
          break;
        c=sDist.proposed.sample(mt);
        cDist=LMLik.get_mcmc_step(lik,model,data,c,sDist.beta);

      }
    return o;

  }

};



template<class mcmc=mcmc_step<>>
class Evidence_Evaluation
{

  std::pair<double,double> logEvidence_;
  std::vector<std::pair<double,SamplesSeries<mcmc>>> run_;
public:

  Evidence_Evaluation(std::vector<std::pair<double,SamplesSeries<mcmc>>>&& o)
  {

    double sum=0;
    double sumVar=0;
    for (int i=0; i<o.size(); ++i)
      {
        SamplesSeries<mcmc>& s=o[i].second;
        std::size_t nsamples=s.size();
        std::pair<double,double> l=s.mean_var
            ([](const mcmc& mc){return mc.logLik;},nsamples/2);
        double betan, betap;
        if (i==0)
          {
            betan=0;
          }
        else
          {
            betan=o[i].first;
          }
        if (i==o.size()-1)
          betap=1;
        else
          betap=o[i+1].first;
        double beta=o[i].first;
        double db=(beta-betan)+0.5*(betap-beta);
        sum+=db*l.first;
        sumVar+=db*l.second;
      }
    logEvidence_={sum,sqrt(sumVar)};
    run_=o;

  }


};





template<
        class D
    , template<class>
    class M
        , template<class,template<class> class >
        class D_Lik=Poisson_DLikelihood
        ,template<class, template<class> class, template<class,template<class> class > class>
        class propDist=LevenbergMarquardt_step
    ,template  < class  , template<class> class  , template<class,template<class> class >class
    ,template<class, template<class> class, template<class,template<class> class > class>
    class >
class MH=Metropolis_Hastings_mcmc>
    class Thermodynamic_Integration_mcmc
{
public:
  Evidence_Evaluation<> *
  run
  (const MH<D,M,D_Lik,propDist>& mcmc
   ,   const propDist<D,M,D_Lik>& LMLik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,const std::vector<std::pair<double, std::pair<std::size_t,std::size_t>>>& beta
      , std::mt19937_64& mt)
  {
    std::vector<std::pair<double,SamplesSeries<mcmc_step<>>>> o;

    M_Matrix<double> pinit=model.sample(mt);
    double beta0=beta[0].first;
    mcmc_step<> cDist=LMLik.get_mcmc_step(lik,model,data,pinit,beta0);

    for (std::size_t i=0; i<beta.size(); ++i)
      {

        o.push_back({beta[i].first,mcmc.run(LMLik,lik,model,data,cDist,beta[i].second.first,beta[i].second.second,beta[i].first,mt)});
       }
    return new Evidence_Evaluation<>(std::move(o));
  }
};



template<class Model,class Data,class HessianApprox=LevenbergMarquardt_gaussian_aprox<Model,Data>>
// prediction(model,data,vector)->vector<pair<double, double>> mean ,std
// prediction (data) -> vector <double>
class LevenbergMarquardt_Function
{
  /// tengo que obtener los siguiente de cada prediccion
  /// 1. el epsilon (y-y(par)) y varianza para cada prediccion.
  /// con estos construyo una mutivariate gaussian.
  ///


  struct trust_region
  {
    double r_;

    bool operator()(double expected, double found, std::size_t k)const
    {
      if (std::pow(expected-found,2)<2*k*r_)
        return true;
      else
        return false;
    }
  };

  HessianApprox c_;
  trust_region t_;
  double landa0_;
  double v_;
  std::size_t maxLoop_;

public:

  LevenbergMarquardt_Function(const Model& m, const Data& d,
                              double  r_trust_region,
                              double landa0=1E8,
                              double v=3,
                              std::size_t maxLoop=10)
    :c_{m,d},t_{r_trust_region},landa0_{landa0},v_(v),maxLoop_(maxLoop)
  {
  }

  Likelihood_Evaluation operator()(const M_Matrix<double>& parameters, double beta)const
  {
    Likelihood_Evaluation s=c_.get_Hessian(parameters,beta);
    double landa=landa0_;
    std::size_t iloop=0;
    s.update_landa(landa);
    double beta_logL=c_.get_betaLogLikelihood(s.next_parameters_,beta);
    while (!t_(s.exp_next_beta_logL_,beta_logL,parameters.size())&&iloop<maxLoop_)
      {
        landa*=v_;
        s.update_landa(landa);
        beta_logL=c_.get_betaLogLikelihood(s.next_parameters_,beta);
        ++iloop;
      }
    s.g_=MultivariateGaussian(s.next_parameters_.toVector(),s.landa_beta_H_,s.landa_beta_Hinv_);
    return s;
  }
};




template<class Model,class Data, class HessianApproximation=LevenbergMarquardt_gaussian_aprox<Model,Data>,class LikelihoodFunction=LevenbergMarquardt_Function<Model,Data,HessianApproximation>
         >
class Metropolis_Hastings
{
  LikelihoodFunction likF_;
  struct test
  {
    test()=default;

    bool operator()(const Likelihood_Evaluation& sLik
                    ,const Likelihood_Evaluation& cLik,std::mt19937_64 &mt)
    {
      logPcurrent=sLik.beta_logLikelihood();
      logPcandidate=cLik.beta_logLikelihood();

      logQforward=sLik.logProb(cLik.sample());

      logQbackward=cLik.logProb(sLik.sample());

      logA=logPcandidate-logQforward-logPcurrent+logQbackward;
      A=std::min(1.0,exp(logA));

      std::uniform_real_distribution<double> u(0,1);
      r=u(mt);
      accept_=r<A;
      return accept_;
    }


    double logPcurrent;
    double logPcandidate;

    double logQforward;

    double logQbackward;

    double logA;

    double A;
    double r;
    bool accept_;

  };

public:

  Metropolis_Hastings(const Model& m, const Data& d,double  r_trust_region)
    :likF_(m,d,r_trust_region){}

  Metropolis_Hastings(const LikelihoodFunction& l):likF_(l){}

  const LikelihoodFunction& logL()const {return likF_;}



  SamplesSeries<Likelihood_Evaluation> run
  (Likelihood_Evaluation& sDist,std::size_t nsamples, std::size_t nskip,double beta,std::mt19937_64& mt)const
  {
    SamplesSeries<Likelihood_Evaluation> o(nsamples);


    if (!sDist.isValid())
      return o;
    std::size_t i=0;

    M_Matrix<double> c=sDist.sample(mt);

    Likelihood_Evaluation cDist=likF_(c,beta);

    test t;

    while(cDist.isValid())
      {
        if (t(sDist,cDist,mt))
          {
            sDist=std::move(cDist);
          }
        if (i%nskip==0)
          o.push_back(sDist);
        if (o.full())
          break;
        c=sDist.sample(mt);
        cDist=likF_(c,beta);

      }
    return o;

  }

};



template<class Likelihood_Eval=Likelihood_Evaluation>
class Evidence_evaluation
{
  std::pair<double,double> logEvidence_;
  std::vector<std::pair<double,SamplesSeries<Likelihood_Eval>>> run_;
public:

  Evidence_evaluation(std::vector<std::pair<double,SamplesSeries<Likelihood_Eval>>>&& o)
  {

    double sum=0;
    double sumVar=0;
    for (int i=0; i<o.size(); ++i)
      {
        SamplesSeries<Likelihood_Eval>& s=o[i].second;
        std::size_t nsamples=s.size();
        std::pair<double,double> l=s.mean_var
            ([](const Likelihood_Eval& l){return l.logLikelihood();},nsamples/2);
        double betan, betap;
        if (i==0)
          {
            betan=0;
          }
        else
          {
            betan=o[i].first;
          }
        if (i==o.size()-1)
          betap=1;
        else
          betap=o[i+1].first;
        double beta=o[i].first;
        double db=(beta-betan)+0.5*(betap-beta);
        sum+=db*l.first;
        sumVar+=db*l.second;
      }
    logEvidence_={sum,sqrt(sumVar)};
    run_=o;

  }


};


template<
    class Model, class Data, class HessianApproximation
    ,class McmcAlgorithm=Metropolis_Hastings<Model,Data,HessianApproximation>
    ,class Likelihood_Eval=Likelihood_Evaluation
    >
class Thermodynamic_Integration
{
public:
  Evidence_evaluation<> *
  run(const Model& m, const Data& d,double  r_trust_region,
      const M_Matrix<double>& x
      ,const std::vector<std::pair<double, std::pair<std::size_t,std::size_t>>>& beta
      , std::mt19937_64& mt)
  {
    McmcAlgorithm mcmc(m,d,r_trust_region);
    std::vector<std::pair<double,SamplesSeries<Likelihood_Eval>>> o;
    Likelihood_Eval cDist=((mcmc.logL())(x,beta[0].first));

    for (std::size_t i=0; i<beta.size(); ++i)
      {

        o.push_back({beta[i].first,mcmc.run(cDist,beta[i].second.first,beta[i].second.second,beta[i].first,mt)});
      }
    return new Evidence_evaluation<>(std::move(o));
  }
};











#endif // EVIDENCE_H
