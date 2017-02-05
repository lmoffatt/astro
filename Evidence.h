#ifndef EVIDENCE_H
#define EVIDENCE_H

#include "Matrix.h"
#include "Distributions.h"

#include <random>

#include <cmath>
#include <list>


class gaussian_lin_regr
{
public:

  double get_optimal_x()const
  {
    double xmean=SX_/SW_;
    double ymean=SY_/SW_;
    double var_x=SXX_/SW_-sqr(xmean);
    double var_y=SYY_/SW_-sqr(ymean);
    double cov=SXY_/SW_-xmean*ymean;
    double b=cov/var_x;
    double a=ymean-b*xmean;
    double out=(y_opt_-a)/b;

    std::cout<<"xmean="<<xmean<<" ymean="<<ymean<<" var_x="<<var_x;
    std::cout<<"var_y="<<var_y<<" cov="<<cov<<" b="<<b<<" a="<<a<<" out= "<<out;


    return out;
  }

  gaussian_lin_regr(double target_y,double target_y_error)
    :y_opt_(target_y), y_e_(target_y_error),SW_(0),SX_(0),SY_(0),SXX_(0),SYY_(0),SXY_(0){}


  void push_back(double x, double y, double ysd)
  {
    double w=1.0/ysd/ysd;
    //      *std::exp(-(std::pow((y-y_opt_),2)/(2.0*(sqr(ysd)+sqr(y_e_)))));
    std::cout<<" yopt= "<<y_opt_<<" y_e_="<<y_e_<<" x= "<<x<<"y ="<<y<<" ysd= "<<ysd<<" w= "<<w<<"\n";

    SW_+=w;
    SX_+=w*x;
    SY_+=w*y;
    SXX_+=w*x*x;
    SYY_+=w*y*y;
    SXY_+=w*x*y;
  }


private:
  double y_opt_;
  double y_e_;
  double SW_;
  double SX_;
  double SY_;
  double SXX_;
  double SYY_;
  double SXY_;


};





inline std::pair<double,double> logit(const std::pair<std::size_t,std::size_t>& x)
{
  std::size_t n=x.first+x.second;
  double p=(1.0+x.first)/(2.0+n);
  double s=sqrt(p*(1-p)/n);
  return logit(p,s);

}


struct D_logL
{
  M_Matrix<double> G;
  M_Matrix<double> H;


};



inline
std::ostream& operator<<(std::ostream& os,const D_logL& d)
{
  os<<"G\n"<<d.G<<"\nH\n"<<d.H;
  return os;
}



struct mcmc_prior
{
  M_Matrix<double> param;
  double logPrior=std::numeric_limits<double>::quiet_NaN();
  D_logL D_prior;
};
inline
std::ostream& operator<<(std::ostream& os,const mcmc_prior& x)
{
  os<<"param\n"<<x.param<<"\nlogPrior\n"<<x.logPrior<<"\nD_prior\n"<<x.D_prior;
  return os;
}


struct mcmc_post: public mcmc_prior
{
  mcmc_post(mcmc_prior &&p): mcmc_prior(p), isValid(false), f(),logLik(0){}
  mcmc_post(){}
  bool isValid=false;
  M_Matrix<double> f;
  double logLik=std::numeric_limits<double>::quiet_NaN();

  double logbPL(double beta)const {return logPrior+logLik*beta;}
};

inline
std::ostream& operator<<(std::ostream& os,const mcmc_post& x)
{
  const mcmc_prior& p=x;

  os<<p<<"\nisValid\n"<<x.isValid;
  if (x.isValid)
    os<<"\f\n"<<x.f<<"\nlogLik\n"<<x.logLik;
  return os;
}


struct mcmc_Dpost: public mcmc_post
{
  mcmc_Dpost(){}
  mcmc_Dpost(mcmc_post &&p): mcmc_post(p), D_lik(){}
  D_logL D_lik;
};
inline
std::ostream& operator<<(std::ostream& os,const mcmc_Dpost& x)
{
  const mcmc_post& p=x;

  os<<p;
  if (x.isValid)
    os<<"\nD_lik\n"<<x.D_lik;
  return os;
}



struct LM_MultivariateGaussian: public MultivariateGaussian
{
  LM_MultivariateGaussian(MultivariateGaussian m, double mylanda):
    MultivariateGaussian(m),landa(mylanda){}
  LM_MultivariateGaussian():landa(){}
  double landa;
};

inline
std::ostream& operator<<(std::ostream& os,const LM_MultivariateGaussian& x)
{
  const MultivariateGaussian& b=x;
  os<<b;
  os<<"\nlanda\n"<<x.landa;
  return os;
}



template<typename Dist>
struct mcmc_step: public mcmc_Dpost
{
  mcmc_step(){}
  mcmc_step(mcmc_Dpost &&p, double beta_): mcmc_Dpost(p),beta(beta_), proposed(){}
  double beta;
  double logbPL()const {return mcmc_post::logbPL(beta);}
  Dist proposed;

};


template<typename Dist>
std::ostream& operator<<(std::ostream& os,const mcmc_step<Dist>& x)
{
  const mcmc_Dpost& b=x;
  os<<b;
  os<<"\nbeta\n"<<x.beta;
  os<<"\nproposed\n"<<x.proposed;

  return os;
}


template <class T>
// T regular type
class SamplesSeries
{
public:
  SamplesSeries(std::size_t n):
    n_(0),samples_(n){}

  SamplesSeries():n_(0),samples_(){}

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

  const std::vector<T>& samples()const {return samples_;}


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

  SamplesSeries()=default;
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


template<typename T>
std::ostream& operator<<(std::ostream& os,const SamplesSeries<T>& x)
{
  os<<"\nsize\n"<<x.size();
  os<<"\nsamples\n"<<x.samples();

  return os;
}


inline double log10_guard(double x)
{
  if (x==0) return 0;
  else return log10(x);

}




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
          {
            if (std::isnan(landa(i,j)))
              return landa(i,j);
            else if (landa(i,j)!=0)
              sumLogL+=logLikelihood(landa(i,j),k(i,j));
          }
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
    return get_mcmc_Dpost(model,data,param,get_mcmc_Post(model,data,param));
  }
  static mcmc_post get_mcmc_Post(const M<D>& model, const D& data, M_Matrix<double> param)
  {
    return Poisson_Likelihood<D,M>::get_mcmc_Post(model,data,param);
  }



  static mcmc_Dpost get_mcmc_Dpost(const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post p)
  {
    mcmc_Dpost out(std::move(p));
    if (out.isValid)
      {
        M_Matrix<std::size_t> k=data();
        M_Matrix<double> logLanda_0=out.f.apply([](double x)
        {return log10_guard(x);});
        if (isnan(logLanda_0))
          out.isValid=false;
        M_Matrix<double> J=get_J(model,  data, param,logLanda_0  );
        if (J.size()==0)
          out.isValid=false;
        else
          {
            out.D_lik.G=get_G(out.f,k,J);
            out.D_lik.H=get_H(out.f,J);
          }
      }
    return out;
  }

private:
  static
  M_Matrix<double> get_G(const M_Matrix<double>& landa, const M_Matrix<std::size_t>& k, const M_Matrix<double>& J)
  {
    M_Matrix<double> out(1,J.ncols(),0.0);
    for (std::size_t j=0; j<J.ncols(); ++j)
      for (std::size_t i=0; i<landa.size(); ++i)
        out[j]+=(landa[i]-k[i])*J(i,j);
    return out;
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
          out(j,j2)+=landa[i]*J(i,j)*J(i,j2);
    for (std::size_t j=0; j<npar; ++j)
      for (std::size_t j2=0; j2<j; ++j2)
        out(j,j2)=out(j2,j);

    return out;

  }

  static
  M_Matrix<double>
  get_J(const M<D>& model, const D& data, const M_Matrix<double>& param,
        const M_Matrix<double>& logLanda_0  ,
        double delta=1e-4, double delta_div=3, double deltamin=1e-8)
  {
    M_Matrix<double> k=data();
    std::size_t n=k.size();
    std::size_t npar=param.size();
    M_Matrix<double> out(n,npar,0.0);
    for (std::size_t j=0; j<npar; ++j)
      {
        double deltarun=delta;

        M_Matrix<double> p(param);
        p[j]+=deltarun;
        M_Matrix<double> logLanda_i=model.logLanda(data,p);
        while (isnan(logLanda_i)&&deltarun>deltamin)
          {
            deltarun=deltarun/delta_div;
            p=param;
            p[j]+=deltarun;
            logLanda_i=model.logLanda(data,p);
          }
        if (isnan(logLanda_i))
          return {};

        for (std::size_t i=0; i<n; ++i)
          out(i,j)=(logLanda_i[i]-logLanda_0[i])/deltarun;

      }
    return out;


  }

};

struct trust_region
{
  static trust_region min(){return {1E-5};}

  static trust_region max(){return {1.0};}

  double r_;

  double getValue()const {return r_;}
  void setValue(double r) { r_=r;}


  bool operator()(double expected, double found, std::size_t k)const
  {
    if (!std::isnan(found)&&(std::pow(expected-found,2)<(2.0*k*r_)))
      return true;
    else
      return false;
  }

};

inline
std::ostream& operator<<(std::ostream& os, const trust_region& t)
{
  os<<t.r_;
  return os;
}


template<
    class D
    , template<class> class M
    , template<class,template<class> class > class D_Lik=Poisson_DLikelihood
    , class myPropDistribution=LM_MultivariateGaussian
    , class AP=trust_region
    >
class LevenbergMarquardt_step
{
public:

  AP t_;
  double landa0_;
  double v_;
  std::size_t maxLoop_;

  double optimal_acceptance_rate_=0.574; // Handbook of Markov Chain Montecarlo p.100
  // optimal acceptance rate for Metropolis-Adjusted Langevin Algorithm



public:

  double optimal_acceptance_rate()const {return optimal_acceptance_rate_;}

  LevenbergMarquardt_step( double landa0,
                           double v,
                           std::size_t maxLoop)
    :t_{AP::min()},landa0_{landa0},v_(v),maxLoop_(maxLoop)
  {
  }


  AP getValue()const{return t_;}



  LevenbergMarquardt_step operator()(AP new_trust_region)const
  {
    LevenbergMarquardt_step o(*this);
    o.t_=std::move(new_trust_region);
    return o;
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
    out.logL=postL.logbPL(beta);
    for (std::size_t i=0; i<out.H.nrows(); ++i)
      //    out.H(i,i)=out.H(i,i)+postL.D_lik.H(i,i)*beta*landa;
      // this alternative does not work
      out.H(i,i)=out.H(i,i)*(1+landa);
    out.Hinv=invSafe(out.H);
    if (!out.Hinv.empty())
      {
        out.d=-(out.G*out.Hinv);
        out.exp_next_logL=out.logL-0.5*multTransp(out.d,out.G)[0];
      }
    return out;
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, double beta)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta);
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,double beta)const
  {
    mcmc_step<myPropDistribution> out(std::move(p),beta);
    return update_mcmc_step(L,model,data,out,beta);
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post&& p,double beta)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),beta);
  }




  mcmc_step<myPropDistribution>& update_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_step<myPropDistribution>& out,double beta)const
  {
    out.beta=beta;


    double landa=0;
    std::size_t iloop=0;
    if (out.isValid)
      {
        LM_logL LM_logL=update_landa( out,landa,beta);
        M_Matrix<double> next;
        mcmc_post cand;
        if (!LM_logL.Hinv.empty())
          {
            next=out.param+LM_logL.d;
            cand=L.get_mcmc_Post(model,data,next);
          }
        while ((LM_logL.Hinv.empty()
                ||!t_(LM_logL.exp_next_logL,cand.logbPL(beta),out.param.size())
                )&&iloop<maxLoop_)

          {
            landa=landa0_*std::pow(v_,iloop);
            LM_logL=update_landa( out,landa, beta);
            next=out.param+LM_logL.d;
            cand=L.get_mcmc_Post(model,data,next);
            ++iloop;
          }
        if (LM_logL.Hinv.empty())
          out.isValid=false;
        else
          out.proposed=myPropDistribution(MultivariateGaussian
                                          (next,LM_logL.Hinv,LM_logL.H),landa);
      }
    return
        out;
  }


};


template<
    class D
    , template<class> class M
    , template<class,template<class> class > class D_Lik
    , class myPropD
    , class myAP
    >
std::ostream& operator<<(std::ostream& os,const LevenbergMarquardt_step<D,M,D_Lik,myPropD,myAP>& x )
{
  os<<"\ntrust_region\n"<<x.t_;
  os<<"\nlanda0\n"<<x.landa0_;
  os<<"\nlanda_mult_factor\n"<<x.v_;
  os<<"\nmaxLoop\n"<<x.maxLoop_;
  return os;
}



template
<
    class D
    , template<class> class M
    , template<class,template<class> class > class D_Lik=Poisson_DLikelihood
    , class my_PropD=LM_MultivariateGaussian
    , class AP=trust_region
    ,template<
      class
      , template<class> class
      , template<class,template<class> class > class
      , class
      , class
      >
    class propDistStep=LevenbergMarquardt_step
    >
class Metropolis_Hastings_mcmc
{
public:
  typedef my_PropD pDist;

  struct test
  {
    test()=default;

    static std::ostream& put(std::ostream& os,const mcmc_step<pDist>& sLik
                             ,const mcmc_step<pDist>& cLik)
    {
      if (cLik.isValid)
        {
          double logPcurrent=sLik.logbPL();
          double logPcandidate=cLik.logbPL();

          double logQforward=sLik.proposed.logP(cLik.param);
          double logQbackward=cLik.proposed.logP(sLik.param);
          double logForward=logPcandidate-logQforward;
          double logBackward=logPcurrent-logQbackward;

          os<<logForward<<" "<<logBackward<<" ";
        }
      return os;
    }

    static bool accept(const mcmc_step<pDist>& sLik
                       ,const mcmc_step<pDist>& cLik
                       ,std::mt19937_64 &mt)
    {
      if (!cLik.isValid)
        {

          return false;
        }
      else
        {
          double logPcurrent=sLik.logbPL();
          double logPcandidate=cLik.logbPL();

          double logQforward=sLik.proposed.logP(cLik.param);
          double logQbackward=cLik.proposed.logP(sLik.param);

          double logA=logPcandidate-logQforward-(logPcurrent-logQbackward);
          double A=std::min(1.0,exp(logA));
          std::uniform_real_distribution<double> u(0,1);
          double r=u(mt);
          bool accept_=r<A;
          return accept_;
        }
    }



  };

public:



  static std::size_t min_tryParameter(){return 5;}




  static SamplesSeries<mcmc_step<pDist>> run
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   std::size_t nsamples,
   std::size_t nskip,
   std::mt19937_64& mt)
  {
    std::size_t nsamples_0=min_tryParameter();
    std::size_t nsamplesFinal=std::max(40lu,nsamples*nskip/20);
    if (!sDist.isValid)
      return {};
    M_Matrix<double> c=sDist.proposed.sample(mt);
    mcmc_step<pDist> cDist=LM_Lik.get_mcmc_step(lik,model,data,c,sDist.beta);

    AP r_optimal=adapt_Parameter
        (LM_Lik,lik,model,data,sDist,cDist,nsamples_0,nsamplesFinal,mt);
    auto LM=LM_Lik(r_optimal);

    SamplesSeries<mcmc_step<pDist>> o(nsamples);
    while(true)
      {
        std::size_t naccepts=0;
        std::size_t nrejects=0;
        n_steps(LM,lik,model,data,sDist,cDist,nskip,naccepts,nrejects,mt);
        o.push_back(sDist);
        if (o.full())
          break;
      }
    return o;

  }


  static std::pair<std::size_t,std::size_t> try_Parameter
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,const AP& r_value
   ,mcmc_step<pDist>& sDist,
   mcmc_step<pDist>& cDist,
   std::size_t nsamples,
   std::mt19937_64& mt)
  {
    std::size_t naccepts=0;
    std::size_t nrejects=0;
    auto LM=LM_Lik(r_value);
    n_steps(LM,lik,model,data,sDist,cDist,nsamples,naccepts,nrejects,mt);
    return {naccepts,nrejects};

  }


  static void n_steps
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   mcmc_step<pDist>& cDist,
   std::size_t nsamples,
   std::size_t& naccepts,
   std::size_t& nrejects,
   std::mt19937_64& mt)
  {
    std::size_t i=0;
    AP r_value=LM.getValue();

    while(i<nsamples)
      {
        std::cout<<i<<" "<<sDist.beta<<" "<<r_value<<" ";
        test::put(std::cout,sDist,cDist);
        put(std::cout,sDist,cDist);
        if (step(LM,lik,model,data,sDist,cDist,mt))
          {
            ++naccepts;
          }
        else
          ++nrejects;
        ++i;
        std::cout<<" "<<naccepts<<" "<<nrejects<<"\n";

      }

  }


  static std::pair<double,double>
  bay_mean_sd(const std::pair<std::size_t,std::size_t> x)
  {
    std::size_t n=x.first+x.second;
    double p=(1.0+x.first)/(2+n);
    double sd=sqrt(p*(1-p)/n);
    return {p,sd};

  }








  static bool accept_Parameter
  (double opt_acc, std::pair<std::size_t, std::size_t> res)
  {
    auto m=bay_mean_sd(res);
    std::cout<<"m="<<m<<"\n";
    std::cout<<"opt_acc="<<opt_acc<<"\n";
    if ((opt_acc>m.first-2*m.second)&&(opt_acc<m.first+2*m.second))
      return true;
    else return false;
  }


  static AP adapt_Parameter
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   mcmc_step<pDist>& cDist,
   std::size_t nsamples_0,
   std::size_t nsamplesFinal,
   std::mt19937_64& mt)
  {
    AP ap0=LM_Lik.getValue();
    double opt_accpt=LM_Lik.optimal_acceptance_rate();
    double y=logit(opt_accpt);
    double yse=std::sqrt(1.0/(opt_accpt*(1-opt_accpt)*nsamplesFinal));
    gaussian_lin_regr lr(y,yse);

    std::size_t nsamples=nsamples_0;
    auto res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt);
    double x0=std::log(ap0.getValue());
    auto m0=logit(res0);
    lr.push_back(x0,m0.first,m0.second );
    ap0=AP::max();

    res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt);
    x0=std::log(ap0.getValue());
    m0=logit(res0);
    lr.push_back(x0,m0.first,m0.second );

    x0=lr.get_optimal_x();
    ap0.setValue(std::exp(x0));
    res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt);

    while (!accept_Parameter(opt_accpt,res0))
      {
        m0=logit(res0);
        lr.push_back(x0,m0.first,m0.second );
        x0=lr.get_optimal_x();
        ap0.setValue(std::exp(x0));
        res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt);
      }
    while (nsamples<nsamplesFinal)
      {
        res0+=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt);
        nsamples*=2;
        while (!accept_Parameter(opt_accpt,res0))
          {
            m0=logit(res0);
            lr.push_back(x0,m0.first,m0.second );
            x0=lr.get_optimal_x();
            ap0.setValue(std::exp(x0));
            res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt);
          }
      }
    return ap0;

  }




  static std::ostream& put
  (std::ostream& os
   ,mcmc_step<pDist>& sLik
   ,mcmc_step<pDist>& cLik)
  {
    os<<sLik.logbPL()<<" "<<cLik.logbPL()<<" ";
    os<<sLik.logLik<<" "<<cLik.logLik<<" ";
    os<<sLik.proposed.landa<<" "<<cLik.proposed.landa;
    return os;
  }



  static bool step
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist
   ,mcmc_step<pDist>& cDist
   ,std::mt19937_64& mt)
  {
    bool out;
    if (test::accept(sDist,cDist,mt))
      {
        sDist=std::move(cDist);
        out =true;
      }
    else
      out=false;
    M_Matrix<double> c=sDist.proposed.sample(mt);
    cDist=LM_Lik.get_mcmc_step(lik,model,data,c,sDist.beta);
    return out;
  }




};






template<class mcmc>
class Evidence_Evaluation
{  
  std::pair<double,double> logEvidence_;
  std::vector<std::pair<double,SamplesSeries<mcmc>>> run_;
public:

  std::pair<double,double> logEvidence()const {return logEvidence_;}
  const std::vector<std::pair<double,SamplesSeries<mcmc>>>& samples()const
  {
    return run_;
  }

  Evidence_Evaluation(std::vector<std::pair<double,SamplesSeries<mcmc>>>&& o)
  {

    double sum=0;
    double sumVar=0;
    for (std::size_t i=0; i<o.size(); ++i)
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

template<class mcmc>
std::ostream& operator<<(std::ostream& os,Evidence_Evaluation<mcmc>& x)
{
  os<<"\nlogEvidence\n"<<x.logEvidence();
  os<<"\nsamples\n"<<x.samples();
  return os;
}



template<
    class D
    , template<class>   class M
    , template<class,template<class> class >   class D_Lik=Poisson_DLikelihood
    , class my_PropD=LM_MultivariateGaussian
    , class AP=trust_region
    ,template<
      class
      , template<class> class
      , template<class,template<class> class > class
      , class
      , class
      >
    class propDist=LevenbergMarquardt_step
    ,template<
      class
      , template<class> class
      , template<class,template<class> class >class
      ,class
      ,class
      ,template<
        class
        , template<class> class
        , template<class,template<class> class > class
        , class
        , class
        >
      class
      >
    class MH=Metropolis_Hastings_mcmc>
class Thermodynamic_Integration_mcmc
{
public:
  typedef mcmc_step<my_PropD> mystep;
  typedef Evidence_Evaluation<mystep> myEvidence;

  static Evidence_Evaluation<mystep> *
  run
  (const MH<D,M,D_Lik,my_PropD,AP,propDist>& mcmc
   ,   const propDist<D,M,D_Lik,my_PropD,AP>& LMLik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,const std::vector<std::pair<double, std::pair<std::size_t,std::size_t>>>& beta
   , std::mt19937_64& mt)
  {
    std::vector<std::pair<double,SamplesSeries<mystep>>> o;


    M_Matrix<double> pinit;
    double beta0=beta[0].first;
    std::size_t ntrials=0;
    mcmc_post postL;
    while(!postL.isValid)
      {
        pinit=model.sample(mt);
        postL=lik.get_mcmc_Post(model,data,pinit);
        ++ntrials;
      }

    mystep cDist=LMLik.get_mcmc_step(lik,model,data,pinit,std::move(postL),beta0);

    for (std::size_t i=0; i<beta.size(); ++i)
      {
        cDist=LMLik.update_mcmc_step(lik,model,data,cDist,beta[i].first);
        auto s=mcmc.run(LMLik,lik,model,data,cDist,beta[i].second.first,beta[i].second.second,mt);
        o.push_back({beta[i].first,s});
      }
    return new Evidence_Evaluation<mystep>(std::move(o));
  }
};







#endif // EVIDENCE_H
