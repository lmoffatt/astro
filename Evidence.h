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
  mcmc_Dpost(){}
  mcmc_Dpost(mcmc_post &&p): mcmc_post(p), D_lik(){}

  D_logL D_lik;
};


template<typename Dist=MultivariateGaussian>
struct mcmc_step: public mcmc_Dpost
{
  mcmc_step(){}
  mcmc_step(mcmc_Dpost &&p, double beta_): mcmc_Dpost(p),beta(beta_), proposed(){}
  double beta;
 double logbPL()const {return mcmc_post::logbPL(beta);}
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
  static double logLikelihood(double landa,std::size_t k)
  {
    return k*log(landa)-landa-lgamma(k+1);
  }


  static double logL(const D& data,const M_Matrix<double>& landa)
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

  static mcmc_post get_mcmc_Post(const M<D>& model, const D& data, M_Matrix<double> param)
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
  static mcmc_Dpost get_mcmc_Dpost(const M<D>& model, const D& data, const M_Matrix<double>& param)
  {
    return get_mcmc_Dpost(model,data,param,Poisson_Likelihood<D,M>::get_mcmc_Post(model,data,param));
  }


  static mcmc_Dpost get_mcmc_Dpost(const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post p)
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
  static
  M_Matrix<double> get_G(const M_Matrix<double>& landa, const M_Matrix<std::size_t>& k, const M_Matrix<double>& J)
  {
    M_Matrix<double> epsilon=landa-M_Matrix<double>(k);
    return epsilon*J;
  }

  static
  M_Matrix<double> get_H(const M_Matrix<double>& landa, const M_Matrix<double>& J)
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

  static
  M_Matrix<double>
  get_J(const M<D>& model, const D& data, const M_Matrix<double>& param,
        const M_Matrix<double>& logLanda_0  ,
        double delta=1e-4)
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

  mcmc_step<> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, double beta)const
  {
     return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta);
  }

   mcmc_step<> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,double beta)const
  {
    mcmc_step<> out(std::move(p),beta);
    return update_mcmc_step(L,model,data,out,beta);
    }

  mcmc_step<>& update_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_step<>& out,double beta)const
  {
    out.beta=beta;

    double landa=landa0_;
    std::size_t iloop=0;
    LM_logL LM_logL=update_landa( out,landa,beta);

    M_Matrix<double> next=out.param+LM_logL.d;
    mcmc_post cand=L.get_mcmc_Post(model,data,next);
    while (!t_(LM_logL.exp_next_logL,cand.logbPL(beta),out.param.size())&&iloop<maxLoop_)
      {
        landa*=v_;
        LM_logL=update_landa( out,landa, beta);
        next=out.param+LM_logL.d;
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
   std::size_t nsamples,
   std::size_t nskip,
   std::mt19937_64& mt)const
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
        cDist=LMLik.update_mcmc_step(lik,model,data,cDist,beta[i].first);
        o.push_back({beta[i].first,mcmc.run(LMLik,lik,model,data,cDist,beta[i].second.first,beta[i].second.second,mt)});
       }
    return new Evidence_Evaluation<>(std::move(o));
  }
};







#endif // EVIDENCE_H
