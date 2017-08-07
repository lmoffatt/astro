#ifndef EVIDENCE_H
#define EVIDENCE_H

#include "Matrix.h"
#include "Distributions.h"
#include "Optimization_BFGS.h"

#include <random>

#include <cmath>
#include <list>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <set>


template<typename T>
std::size_t ncells(const std::vector<T>& v)
{
  std::size_t count=0;
  for (auto &e:v)
    count+=e.size();
  return count;
}
template<typename T, class Predicate>
std::size_t ncells(const std::vector<T>& v,const  Predicate& p)
{
  std::size_t count=0;
  for (auto &e:v)
    if (p(e)) ++count;
  return count;
}

template<typename T, class Predicate>
std::size_t ncells(const std::vector<std::vector<T>>& v,const  Predicate& p)
{
  std::size_t count=0;
  for (auto &e:v)
    count+=ncells(e,p);
  return count;
}


class gaussian_lin_regr
{
public:

  double regression_coefficient()const
  {
    double xmean=SX_/SW_;
    double ymean=SY_/SW_;
    double var_x=SXX_/SW_-sqr(xmean);
    double cov=SXY_/SW_-xmean*ymean;
    double b=cov/var_x;
    return b;
  }


  double get_optimal_x(double target_y, std::ostream& os)const
  {
    double xmean=SX_/SW_;
    double ymean=SY_/SW_;
    double var_x=SXX_/SW_-sqr(xmean);
    double var_y=SYY_/SW_-sqr(ymean);
    double cov=SXY_/SW_-xmean*ymean;
    double b=cov/var_x;
    double a=ymean-b*xmean;
    double out=(target_y-a)/b;

    std::cout<<"get_optimal_x::\txmean\t"<<xmean<<"\tymean\t"<<ymean<<"\tvar_x\t"<<var_x;
    std::cout<<"\tvar_y\t"<<var_y<<"\tcov\t"<<cov<<"\tb\t"<<b<<"\ta\t"<<a<<"\tout\t"<<out<< "\texp(out)\t"<<std::exp(out)<<"\n";

    os<<"get_optimal_x::\txmean\t"<<xmean<<"\tymean\t"<<ymean<<"\tvar_x\t"<<var_x;
    os<<"\tvar_y\t"<<var_y<<"\tcov\t"<<cov<<"\tb\t"<<b<<"\ta\t"<<a<<"\tout\t"<<out<< "\texp(out)\t"<<std::exp(out)<<"\n";

    return std::max(min_x_,std::min(max_x,out));
  }

  gaussian_lin_regr(double min_x,double max_x)
    :min_x_(std::min(min_x,max_x)),max_x(std::max(min_x,max_x)),SW_(0),SX_(0),SY_(0),SXX_(0),SYY_(0),SXY_(0){}

  gaussian_lin_regr()
    :min_x_(),max_x(),SW_(0),SX_(0),SY_(0),SXX_(0),SYY_(0),SXY_(0){}

  void push_back(double x, double y, double ysd)
  {
    double w=1.0/ysd/ysd;
    //      *std::exp(-(std::pow((y-y_opt_),2)/(2.0*(sqr(ysd)+sqr(y_e_)))));
    std::cout<<" x= "<<x<<"y ="<<y<<" ysd= "<<ysd<<" w= "<<w<<"\n";

    SW_+=w;
    SX_+=w*x;
    SY_+=w*y;
    SXX_+=w*x*x;
    SYY_+=w*y*y;
    SXY_+=w*x*y;
  }


private:
  double min_x_;
  double max_x;
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

  double d_logLik_dBeta(double beta)const
  {
    auto Hbinv=inv(D_prior.H+beta*D_lik.H);
    auto Gb=D_prior.G+beta*D_lik.G;
    auto db=Gb*Hbinv;
    auto GdH=D_lik.G-db*D_lik.H;
    double s=xTSigmaX(GdH,Hbinv);
    //  for (std::size_t i=0; i<D_lik.H.nrows(); ++i)
    //     s+=sqr(D_lik.H(i,i));
    return s;
  }
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






template<typename Dist>
struct mcmc_step: public mcmc_Dpost
{
  mcmc_step(){}
  mcmc_step(mcmc_Dpost &&p, double beta_, std::size_t iscout_): mcmc_Dpost(p),beta(beta_),iscout(iscout_),proposed(){}
  double beta;
  std::size_t iscout;

  double logbPL()const {return mcmc_post::logbPL(beta);}
  double logbPLb(double mybeta)const {return mcmc_post::logbPL(mybeta);}
  Dist proposed;

  std::ostream& writelogLHeaderDataFrame(std::ostream& os)
  {
    os<<"beta\t";
    os<<"nscout\t";
    os<<"logPrior\t";
    os<<"logLik\t";
    os<<"logP";
    return os;
  }
  std::ostream& writelogLRowDataFrame(std::ostream& os)
  {
    os<<beta<<"\t";
    os<<iscout<<"\t";
    os<<logPrior<<"\t";
    os<<logLik<<"\t";
    os<<logbPL();
    return os;
  }


  std::ostream& writeParamRowDataFrame(std::ostream& os)
  {
    os<<"\t";
    for (std::size_t i=0; i<param.size(); ++i)
      os<<param[i]<<"\t";
    return os;
  }

  std::ostream& writeYfitRowDataFrame(std::ostream& os)
  {

    for (std::size_t i=0; i+1<f.size(); ++i)
      os<<f[i]<<"\t";
    if (f.size()>0)
      os<<f[f.size()-1];
    return os;
  }



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
    double s=0;
    for (std::size_t i=i0; i<n_; ++i)
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
    double s=0;
    for (std::size_t i=i0; i<n_; ++i)
      s+=sqr(f(samples_[i]));
    double var=s/(n_-i0)-sqr(m);
    return {m,var};
  }

  template<class F>
  std::pair<double,double> mean_std(const F& f, std::size_t i0=0)
  {
    auto o=mean_var(f,i0);
    return {o.first,std::sqrt(o.second)};
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
  std::size_t n_;
  std::vector<double> samples_;
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
    if (out.f.size()==0)
      out.isValid=false;
    else
      {
        out.logLik=logL(data,out.f);
        if (std::isnan(out.logLik))
          out.isValid=false;
        else out.isValid=true;
      }
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
        double delta=1e-5, double delta_div=10, double deltamin=1e-7)
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
        while ((isnan(logLanda_i)||(logLanda_i.empty()))&&deltarun>deltamin)
          {
            deltarun=deltarun/delta_div;
            p=param;
            p[j]+=deltarun;
            logLanda_i=model.logLanda(data,p);
          }
        if (isnan(logLanda_i)||logLanda_i.empty())
          return {};

        for (std::size_t i=0; i<n; ++i)
          out(i,j)=(logLanda_i[i]-logLanda_0[i])/deltarun;

      }
    return out;


  }

};

struct trust_region
{
  static trust_region min(){return {1E-6};}

  static trust_region max(){return {1.0};}

  static double logFactor(){return std::log(100.0);}

  double r_;

  double getValue()const {return r_;}
  void setValue(double r) { r_=r;}


  bool operator()(double expected, double found, std::size_t k)const
  {
    if (!std::isnan(found)&&(sqr(expected-found)<2.0*r_*k))
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

struct landa_report
{
  std::size_t count;
  std::size_t nanCount;
  double sumAcceptance;
  double sumSqrAcceptance;

  double sumDeltaLogL;
  double sumSqrDeltaLogL;
};

inline
std::ostream& operator<<(std::ostream& os, const landa_report& l)
{
  os<<l.count<<"\t"<<l.nanCount<<"\t"<<l.sumAcceptance<<"\t"<<l.sumSqrAcceptance<<"\t"<<l.sumDeltaLogL<<"\t"<<l.sumSqrDeltaLogL<<"\n";
  return os;
}



template <typename T>
T sample_rev_map(const std::map<double,T>& reverse_prior,std::mt19937_64& mt)
{
  std::uniform_real_distribution<> u;
  double r= u(mt);
  auto it=reverse_prior.lower_bound(r);
  return it->second;

}

template <class T>
std::map<double,T>
cumulative_reverse_map(const std::map<T,double>& normalized_map)
{
  std::map<double,T> out;
  double sump=0;
  for (auto& e:normalized_map)
    {
      sump+=e.second;
      out[sump]=e.first;
    }
  return out;
}





struct BayesIterator
{

  template <typename Data,  class T, class LikelihoodFunction_with_count>

  static std::map<T,double>
  posterior(const LikelihoodFunction_with_count& f,const Data& newData, std::map<T,double> prior,  double & Evidence, std::size_t& n)
  {
    Evidence=0;
    std::map<T,double> out(std::move(prior));
    for (auto it=out.begin(); it!=out.end(); ++it)
      {
        double lik=f(newData,it->first,n);
        it->second*=lik;
        Evidence+=it->second;
      }
    for (auto it=out.begin(); it!=out.end(); ++it)
      {
        it->second/=Evidence;

      }
    return out;
  }

  template <typename Data,  class T, class LikelihoodFunction_with_count>
  static std::pair<std::map<T,double>,std::size_t>
  multinomial_posterior(const LikelihoodFunction_with_count& f,const Data& newData, std::pair<std::map<T,double>, std::size_t> prior, double & Evidence)
  {
    std::map<T,double> out(prior.first);
    std::size_t count=prior.second;
    std::size_t n;
    std::map<T,double> sample=posterior(f,newData,prior.first,Evidence,n);

    for (auto it=out.begin(); it!=out.end(); ++it)
      {
        it->second=(it->second*count +sample[it->first])/(count+n);
      }
    count+=n;

    return {out,count};
  }

};







struct OptimalDistribution
{

  template<typename T, class GainFunction, class Data>
  static std::map<T,double> optimal(const GainFunction& g
                                    ,const Data& data
                                    , const std::map<T,double>& initDistribution,
                                    std::size_t count)
  {

    auto init=to_logit(initDistribution);
    auto g_logit=logit_to_Function<GainFunction,T,Data>(g,initDistribution,count);
    std::map<T,double> out;
    opt_max_iter res;
    // M_Matrix<double> init=init0;
    res=BFGS_optimal::opt(g_logit,data,init);
    out=logit_to_distribution(res.sample.b,initDistribution);
    //  std::cerr<<"res\n"<<res;
    //  std::cerr<<out;
    return out;

  }

  template<typename T>
  static M_Matrix<double> to_logit(const std::map<T,double>& dist)
  {
    M_Matrix<double> out(1,dist.size()-1);
    double q=1;
    auto it=dist.begin();
    for (std::size_t i=0; i<out.size(); ++i)
      {
        double p=it->second;
        out[i]=logit(p/q);
        q-=p;
        ++it;
      }
    return out;

  }

  template<typename T>
  static std::map<T,double>  logit_to_distribution(const M_Matrix<double>& logParam
                                                   ,const std::map<T,double>& dist)
  {
    std::map<T,double> out(dist);
    double q=1.0;
    auto it=out.begin();
    for (std::size_t i=0; i<logParam.size(); ++i)
      {
        double p=logistic(logParam[i]);
        it->second=p*q;
        q-=it->second;
        ++it;
      }
    it->second=q;
    return out;


  }

  template <class GainFunction, typename T,class Data>
  struct logit_to_Function
  {
    const GainFunction& g_;
    const std::map<T,double>& initDistribution;
    M_Matrix<double> bmean;
    M_Matrix<double> bstd;
    std::size_t count;

    logit_to_Function(const GainFunction& g,
                      const std::map<T,double>& initDistribution_,
                      std::size_t count_)
      :g_(g)
      ,initDistribution(initDistribution_)
      ,bmean(to_logit(initDistribution))
      ,bstd()
      ,count(count_){
      bstd=ones(bmean);
    }

    double operator()(const Data& data
                      ,const M_Matrix<double>& logitValues)const
    {
      auto lo=logit_to_distribution(logitValues,initDistribution);
      double prior=Normal(logitValues,bmean,bstd);
      return log(g_(data,lo))*count-prior;
    }

  };
};



struct Optimal_Lagrange_Distribution
{

  template<typename T, class GainFunction, class Data>
  static std::map<T,double> optimal(const GainFunction& g
                                    ,const Data& data
                                    , const std::map<T,double>& initDistribution,
                                    std::size_t count)
  {
    auto init=to_logit(initDistribution);
    auto g_logit=logit_to_landa_Function<GainFunction,T,Data>(g,initDistribution,count);
    std::map<T,double> out;
    opt_max_iter res;
    // M_Matrix<double> init=init0;
    res=BFGS_optimal::opt(g_logit,data,init);
    if (res.sample.G.empty()) return {};
    else
      {
        out=logit_to_distribution(res.sample.b,initDistribution);
        //  std::cerr<<"res\n"<<res;
        //  std::cerr<<out;
        return out;

      }
  }

  template<typename T>
  static M_Matrix<double> to_logit(const std::map<T,double>& dist)
  {
    M_Matrix<double> out(1,dist.size());
    auto it=dist.begin();
    for (std::size_t i=0; i<out.size(); ++i)
      {
        double p=it->second;
        out[i]=logit(p);
      }
    return out;

  }

  template<typename T>
  static M_Matrix<double> to_logit_landa(const std::map<T,double>& dist,double landa)
  {
    M_Matrix<double> out(1,dist.size()+1);
    auto it=dist.begin();
    for (std::size_t i=0; i<dist.size(); ++i)
      {
        double p=it->second;
        out[i]=logit(p);
      }
    out[dist.size()]=landa;
    return out;

  }


  template<typename T>
  static std::map<T,double>
  logit_to_distribution(const M_Matrix<double>& logit_landa_Param
                        ,const std::map<T,double>& dist)
  {
    std::map<T,double> out(dist);
    auto it=out.begin();
    for (std::size_t i=0; i<logit_landa_Param.size()-1; ++i)
      {
        double p=logistic(logit_landa_Param[i]);
        it->second=p;
      }
    return out;
  }

  template <class GainFunction, typename T,class Data>
  struct logit_to_landa_Function
  {
    const GainFunction& g_;
    const std::map<T,double>& initDistribution;
    M_Matrix<double> bmean;
    M_Matrix<double> bstd;
    std::size_t count;

    logit_to_landa_Function(const GainFunction& g,
                            const std::map<T,double>& initDistribution_,
                            std::size_t count_)
      :g_(g)
      ,initDistribution(initDistribution_)
      ,bmean(to_logit(initDistribution))
      ,bstd()
      ,count(count_){
      bstd=ones(bmean);
    }

    double operator()(const Data& data
                      ,const M_Matrix<double>& logit_landa_Values)const
    {

      auto lo=logit_to_distribution(logit_landa_Values,initDistribution);
      double prior=Normal(logit_landa_Values,bmean,bstd);
      double landa=logit_landa_Values[logit_landa_Values.size()];
      double sum_p=0;
      for (auto& e:lo)
        sum_p+=e.second;
      return log(g_(data,lo))*count-prior-landa*(sum_p-1);
    }

  };
};







class Landa
{
public:
  typedef std::pair<Landa,double> myParameter;


  static std::string ParName(std::size_t i)
  {
    switch (i)
      {
      case 0: return "Landa_50";
      case 1: return "Hill_Coeff";
      default: return "";
      }
  }


  static std::string ClassName(){return "Landa";}
  static trust_region min(){return {0.0};}

  static trust_region max(){return {1E9};}

  static double logFactor(){return std::log(10.0);}



  double landa_;

  double getValue()const {return landa_;}
  void setValue(double r) { landa_=r;}


  struct myAcceptProb
  {
    double operator()(Landa landa,const myParameter& param)const
    {
      double landa50=param.first.getValue();
      double h=param.second;
      return 1.0/(1.0+std::pow(landa50/(landa.getValue()+1),h));
    }
  };
  typedef myAcceptProb AcceptanceProbability;

  struct myExpectVelocity
  {
    double operator()(const Landa& landa)const
    {
      return (1.0/(1.0+landa.getValue()));
    }
  };

  typedef myExpectVelocity ExpectedVelocity;

  static std::map<myParameter, double> uniform_parameter_prior(const  std::vector<std::vector<double>>& v)
  {
    std::map<myParameter, double> out;
    double p=1.0/(v[0].size()*v[1].size());
    for (std::size_t i=0; i<v[0].size(); ++i)
      for (std::size_t j=0; j<v[1].size(); ++j)
        out[{Landa{v[0][i]},v[1][j]}]+=p;
  return out;
}


};


inline bool operator<(const Landa& one,const Landa& two)
{ return one.getValue()<two.getValue();}



template <class AP=Landa>
class Adaptive_discrete
{
public:
  AP sample(std::mt19937_64& mt)const
  {
    return sample_rev_map(rev_,mt);
  }

  void push_acceptance(AP landa,double dHd)
  {
    std::pair<std::multiset<AP>,AP> p{std::move(currentRejected_),landa};
    currentRejected_.clear();

    std::get<0>(rejAccCount_[p])++;
    std::get<1>(rejAccCount_[p])+=dHd;
    std::get<2>(rejAccCount_[p])+=dHd*dHd;

    double Evidence=0;
    parDist_=BayesIterator::multinomial_posterior
        (&likelihood,p,std::move(parDist_),Evidence);
  }


  void push_rejection(AP landa)
  {
    currentRejected_.insert(landa);

  }

  void actualize()
  {
    auto pold=this->p_;
    this->p_=OptimalDistribution::optimal(&expectedVelocity,parDist_.first,this->p_,parDist_.second);
    auto p2=Optimal_Lagrange_Distribution::optimal(&expectedVelocity,parDist_.first,pold,parDist_.second);

    double sum_p=0;
    for (auto& e:p2)
      sum_p+=e.second;

    if (false)
      {
        std::cerr<<"p logit"<<p_<<"\n";
        std::cerr<<"p lagrange"<<p2<<"\n"<<sum_p<<"\n";
      }

    this->rev_=cumulative_reverse_map(this->p_);
  }

  static double likelihood(std::pair<std::multiset<AP>,AP> data,
                           const typename AP::myParameter& par
                           , std::size_t& n)
  {
    double p=1;
    n=0;
    typename AP::myAcceptProb pA;
    for (const AP& landa:data.first)
      {
        p*=(1.0-pA(landa,par));
        n++;
      }
    p*=pA(data.second,par);
    n++;
    return p;

  }

  static double expectedVelocity(const std::map<typename AP::myParameter, double>& parDist,
                                 const std::map<AP,double>& pAp)
  {
    double sum=0;
    typename AP::ExpectedVelocity E;
    typename AP::AcceptanceProbability AcP;

    for (auto& e1:parDist)
      {
        double sum2=0;
        for (auto&e2: pAp)
          {
            sum2+=e2.second*AcP(e2.first,e1.first)*E(e2.first);
          }
        sum+=e1.second/sum2;
      }
    return sum;
  }


  Adaptive_discrete(const std::map<AP,double>& prior_landa,
                    const std::map<typename AP::myParameter, double>& prior_par):
    p_{prior_landa},parDist_{prior_par,1},
    rev_{cumulative_reverse_map(p_)},
    rejAccCount_{}{}

  Adaptive_discrete(std::map<AP,double>&& prior_landa,
                    std::map<typename AP::myParameter, double>&& prior_par):
    p_{std::move(prior_landa)},parDist_{std::move(prior_par),1},
    rev_{},
    rejAccCount_{}{
    rev_=cumulative_reverse_map(p_);
  }

  template<template<typename>class V>
  Adaptive_discrete(const V<AP>& landa,
                    const  std::vector<std::vector<double>>& par):
    Adaptive_discrete(uniform_prior(landa),AP::uniform_parameter_prior(par)){}


  Adaptive_discrete partialReset()const
  {
    return Adaptive_discrete(p_,parDist_.first);
  }

  friend
  std::ostream& operator<<(std::ostream& os, const Adaptive_discrete<AP>& me)
  {
    os<<AP::ClassName()<<" distribution\n";
    for (auto &e:me.p_)
      {
        os<<e.first<<"\t"<<e.second<<"\n";
      }
    os<<AP::ClassName()<<" parameter distribution\n";
    os<<"count\t"<<me.parDist_.second<<"\n";
    for (auto &e:me.parDist_.first)
      {
        os<<e.first<<"\t"<<e.second<<"\n";
      }

    os<<AP::ClassName()<<" reverse distribution\n";
    for (auto &e:me.rev_)
      {
        os<<e.first<<"\t"<<e.second<<"\n";
      }

    os<<AP::ClassName()<<" currently rejected\n"<<me.currentRejected_<<"\n";
    os<<AP::ClassName()<<" history of rejected accepted\n";
    for (auto &e:me.rejAccCount_)
      {
        os<<e.first<<"\t"<<e.second<<"\n";
      }
    return os;

  }

  friend
  std::ostream& operator<<(std::ostream& os, const std::vector<Adaptive_discrete<AP>>& me)
  {
    os<<AP::ClassName()<<" distribution\n";

    for (auto &e:me[0].p_)
      {
        auto a=e.first;
        os<<a<<"\t";
        for (std::size_t i=0; i<me.size(); ++i)
          {
            auto it=me[i].p_.find(a);
            os<<*it<<"\t";
          }
        os<<"\n";
      }
    os<<AP::ClassName()<<" parameter distribution\n";
    os<<"count\t"<<me[0].parDist_.second<<"\n";
    for (auto &e:me[0].parDist_.first)
      {
        auto a=e.first;
        os<<a<<"\t";
        for (std::size_t i=0; i<me.size(); ++i)
          {
            auto it=me[i].parDist_.first.find(a);
            os<<*it<<"\t";
          }
        os<<"\n";

      }

    os<<AP::ClassName()<<" reverse distribution\n";
    std::vector<typename std::map<double,AP>::const_iterator> its(me.size());
    for (std::size_t i=0; i<me.size(); ++i)
      its[i]=me[i].rev_.begin();

    while (its[0]!=me[0].rev_.end())
      {
        for (std::size_t i=0; i<me.size(); ++i)
          {
            os<<i<<" "<<its[i]->first<<"\t"<<its[i]->second<<"\t";
            ++its[i];
          }
        os<<"\n";
      }
    os<<AP::ClassName()<<" history of rejected accepted\n";
    std::vector<std::map<AP,std::tuple<double,double,std::size_t>>> hist(me.size());
    for (std::size_t i=0; i<me.size(); ++i) hist[i]=me[i].history();
    for (auto &e:hist[0])
      {
        auto a=e.first;
        os<<a<<"\t";
        for (std::size_t i=0; i<hist.size(); ++i)
          os<<i<<": "<<hist[i][a]<<"\t";
        os<<"\n";
      }
    return os;

  }

  std::map<AP,std::tuple<double,double,std::size_t>>
  history()const
  {
    std::map<AP,std::tuple<std::size_t,std::size_t,double>> o;
    for (auto &e:rejAccCount_)
      {
        for (auto & rej:e.first.first)
          std::get<1>(o[rej])+=std::get<0>(e.second);
        std::get<0>(o[e.first.second])+=std::get<0>(e.second);
        std::get<2>(o[e.first.second])+=std::get<1>(e.second);
      }
    std::map<AP,std::tuple<double,double,std::size_t>> out;
    for (auto &e: o)
      {
        std::size_t n=std::get<0>(e.second)+std::get<1>(e.second);
        double pAcc=std::get<0>(e.second)/n;
        double movMedio=std::get<2>(e.second)/std::get<0>(e.second);
        out[e.first]={pAcc,movMedio,n};
      }
    return out;

  }

private:
  std::map<AP,double> p_;

  std::pair<std::map<typename AP::myParameter, double>,std::size_t> parDist_;
  std::map<double,AP> rev_;

  std::multiset<AP> currentRejected_;
  std::map<std::pair<std::multiset<AP>,AP>,std::tuple<std::size_t,double,double>> rejAccCount_;

  template<template<typename>class V>
  static std::map<AP,double> uniform_prior(const V<AP>& v)
  {
    std::map<AP,double> out;
    double p=1.0/v.size();
    for (auto it=v.begin(); it!=v.end(); ++it)
      out[*it]+=p;
    return out;
  }



};





inline
std::ostream& operator<<(std::ostream& os, const Landa& t)
{
  os<<t.getValue();
  return os;
}

inline
std::istream& operator>>(std::istream& is, Landa& t)
{
  double landa;
  is>>landa;
  t.setValue(landa);
  return is;
}



class Beta

{
public:
  static std::string ClassName(){return "Beta";}


  Beta(std::size_t n, double med_beta, std::size_t n2=4, double min_beta=1e-7):asc_beta_(n+n2)
  {
    double f=std::pow(med_beta,1.0/(n-1));
    double f2=std::pow(med_beta/min_beta,1.0/n2);

    for (std::size_t i=0; i<n; ++i)
      asc_beta_[i]=(std::pow(f,i)+(n-i-1)*med_beta)/(1.0+(n-i-1)*med_beta);
    for (std::size_t i=0; i<n2; ++i)
      asc_beta_[n+i]=med_beta/(std::pow(f2,i+1));


  }

  Beta():asc_beta_(){}

  Beta(std::vector<double>&& beta):asc_beta_(beta){}


  std::vector<double>const & getValue()const {return asc_beta_;}

  std::vector<double>& getValue() {return asc_beta_;}

  std::size_t size()const {return asc_beta_.size();}

  static
  std::vector<double> resize(std::size_t newsize, const std::vector<double>& beta)
  {
    std::vector<double> newBeta(newsize);
    double f=1.0*beta.size()/newsize;
    for (std::size_t i=0; i<newsize; ++i)
      newBeta[i]=interpolate(i*f,beta);
    return newBeta;
  }

private:
  std::vector<double> asc_beta_;

  static double interpolate(double n, const std::vector<double>& beta)
  {
    std::size_t i=std::floor(n);
    double f=n-i;
    if (i+1<beta.size())
      return beta[i]+f*beta[i+1];
    else return beta[i];


  }
};

inline
std::ostream& operator<<(std::ostream& os, const Beta& t)
{
  os<<t.getValue();
  return os;
}

inline
std::istream& operator>>(std::istream& is, Beta& t)
{
  std::vector<double> betas;
  is>>betas;
  t.getValue()=betas;
  return is;
}





// implements
//Dynamic temperature selection for parallel tempering in Markov chain Monte Carlo simulations W. D. Vousden, W. M. Farr and I. Mandel
//doi:10.1093/mnras/stv2422
//https://arxiv.org/abs/1501.05823

class Adaptive_Beta
{
public:

  void push_acceptance(std::size_t i)
  {
    ++accepts_[i].first;
  }

  void push_rejection(std::size_t i)
  {
    ++accepts_[i].second;
  }





  void actualize(std::size_t n)
  {
    nsamples+=n;
    auto S=get_S(beta_);
    auto A=get_A(accepts_);
    double k=get_k(nu_,t0_,nsamples);
    for (std::size_t i=0; i< S.size(); ++i)
      {
        S[i]+=k*(A[i]-A[i+1]);
      }
    beta_=to_Beta(beta_,S);

    clear_accepts(accepts_);
  }


  friend
  std::ostream& operator<<(std::ostream& os, const Adaptive_Beta& me)
  {
    auto A=get_A(me.accepts_);
    os<<Beta::ClassName()<<" acceptance\n";
    for (std::size_t i=0; i<A.size(); ++i)
      {
        os<<me.beta_.getValue()[i]<<"\t"<<me.beta_.getValue()[i+1]<<"\t"<<A[i]<<"\t"<<me.accepts_[i]<<"\n";
      }

    return os;

  }

  static double get_k(double nu, std::size_t nsamples50, std::size_t nsamples)
  {
    return 1.0/nu*nsamples50/(nsamples50+nsamples);
  }

  Adaptive_Beta(std::size_t N,double beta_min,double nu, std::size_t nsamples50)

    :nsamples(0),
      beta_min_(beta_min),  // minimal value of beta
      N_(N), // number of beta intervals
      nu_{nu},  // reciprocal of the initial amplitude of adjustments
      t0_{nsamples50},  //lag parameter
      beta_{N,beta_min},
      accepts_{N-1,{0,0}}{}


  std::size_t size()const {return N_;}

  Beta const& getBeta()const {return beta_;}


private:
  std::size_t nsamples;
  double beta_min_;  // minimal value of beta
  std::size_t N_; // number of beta intervals
  double nu_;  // reciprocal of the initial amplitude of adjustments
  std::size_t t0_;  //lag parameter

  Beta beta_;
  std::vector<std::pair<std::size_t,std::size_t>> accepts_;



  static std::vector<double> get_S(const Beta& b)
  {
    auto n=b.getValue().size();
    std::vector<double> out(n-2);
    for (std::size_t i=0; i< out.size(); ++i)
      out[i]=std::log(1.0/b.getValue()[i+1]-1.0/b.getValue()[i]);
    return out;
  }
  static Beta& to_Beta(Beta& b, std::vector<double> S)
  {
    for (std::size_t i=0; i< S.size(); ++i)
      b.getValue()[i+1]=1.0/(1.0/b.getValue()[i]+std::exp(S[i]));
    return b;
  }

  static std::vector<double> get_A(const std::vector<std::pair<std::size_t,std::size_t>>& accept)
  {
    std::vector<double> o(accept.size(),0);
    for (std::size_t i=0; i<o.size(); ++i)
      if ((accept[i].first+accept[i].second)>0)
        o[i]=1.0*accept[i].first/(accept[i].first+accept[i].second);
    return o;
  }

  static void clear_accepts(std::vector<std::pair<std::size_t,std::size_t>>& accep)
  {
    for (auto&e: accep)
      {
        e.first=0;
        e.second=0;
      }
  }
};




// implements
// usa

class Luciano_Adaptive_Beta
{
public:


  void push_back(double beta, double logLik)
  {
    std::get<0>(data_[beta])+=logLik;
    std::get<1>(data_[beta])+=logLik*logLik;
    ++std::get<2>(data_[beta]);
  }
  void push_acceptance(double betaless, double betamore)
  {
    ++accepts_[{betaless,betamore}].first;
  }

  void push_rejection(double betaless, double betamore)
  {
    ++accepts_[{betaless,betamore}].second;
  }

  friend
  std::ostream& operator<<(std::ostream& os, const Luciano_Adaptive_Beta& me)
  {
    os<<"Beta\n"<<me.beta_;
    os<<"History of likelihoods\n";

    for (auto &e: me.data_)
      {
        double x=e.first;
        double y=std::get<0>(e.second)/std::get<2>(e.second);
        double sd=std::sqrt(std::get<1>(e.second)/std::get<2>(e.second)-y*y);
        std::size_t n=std::get<2>(e.second);
        os<<x<<"\t"<<y<<"\t"<<sd<<"\t"<<n<<"\n";
      }
    return os;

  }

  void actualize()
  {
    std::vector<double> betanew;
    double be_0=1;
    betanew.push_back(be_0);
    double be_1=beta_.getValue()[1];
    for (std::size_t i=1; i< beta_.size(); ++i)
      {
        double be0=beta_.getValue()[i-1];
        double be1=beta_.getValue()[i];
        double be2;
        if (i+1<beta_.size())
          be2=beta_.getValue()[i+1];
        else
          be2=0;
        const auto it0=data_.lower_bound(be0);
        const auto it2=data_.upper_bound(be2);
        double dlogdb=d_logLik_d_beta(it0,it2);
        double db=1.0/std::sqrt(dlogdb);
        be_1=be_0-db*factor_;
        double be12=(be1+be2)/2;
        while (be_1>be12)
          {
            betanew.push_back(be_1);
            be_0=be_1;
            be_1=be_0-db*factor_;
          }
      }
    if (beta_.size()>betanew.size())
      betanew=Beta::resize(beta_.size(),betanew);
    else if (betanew.size()>Nmax_)
      betanew=Beta::resize(Nmax_,betanew);

    beta_=Beta(std::move(betanew));
  }


  void reset()
  {
    data_.clear();
  }

  Luciano_Adaptive_Beta(std::size_t Ninitial, std::size_t Nmax, double beta_min,double factor=1): Nmax_(Nmax),factor_(factor),beta_{Ninitial,beta_min}, data_()
  {}

  void init(const mcmc_Dpost s, double factor=1)

  {
    beta_=getBeta(s,factor);
  }


  std::size_t size()const {return beta_.size();}

  Beta const& getBeta()const {return beta_;}


private:
  std::size_t Nmax_;
  double factor_;

  Beta beta_;

  std::map<double, std::tuple<double, double, std::size_t>, std::greater<double> >  data_;
  std::map<std::pair<double,double>,std::pair<std::size_t,std::size_t>> accepts_;


  template<class It>
  static
  double d_logLik_d_beta(const It& begin,const It& end)
  {
    gaussian_lin_regr lr;
    for (It it=begin; it!=end; ++it)
      {
        double x=it->first;
        double y=std::get<0>(it->second)/std::get<2>(it->second);
        double sd=std::sqrt((std::get<1>(it->second)/std::get<2>(it->second)-y*y)/std::get<2>(it->second));
        lr.push_back(x,y,sd);
      }
    return lr.regression_coefficient();
  }




  static
  Beta getBeta(const mcmc_Dpost s, double factor)
  {
    std::vector<double> b;

    double be=1;
    while (be>0)
      {
        b.push_back(be);
        double db=1.0/std::sqrt(s.d_logLik_dBeta(be));
        be-=db*factor;
      }
    return Beta(std::move(b));

  }


};



struct Likelihood_Record_steps{
  struct Particle
  {
    std::vector<std::size_t> ibeta;
    double logLikInit;
    double logLikNext;
  };
  std::size_t size()const {return beta.size();}



  M_Matrix<double> beta;
  std::vector<Particle> particle;
  std::map<double,std::size_t> iParticle_init;
  std::map<double,std::size_t> iParticle_end;
  std::map<double,std::size_t> Beta_to_i;

  template<typename mcmc>
  Likelihood_Record_steps(const Beta& b,const std::vector<mcmc>& sDist):
    beta(1,b.getValue().size(),b.getValue()),
    particle(sDist.size())
  {
    for (std::size_t i=0; i<sDist.size(); ++i)
      {
        particle[i].ibeta.push_back(i);
        particle[i].logLikInit=sDist[i].logLik;
        iParticle_init[beta[i]]=i;
        iParticle_end[beta[i]]=i;
        Beta_to_i[beta[i]]=i;
      }
  }
  template<typename mcmc>
  void end_record(const Beta& b,const std::vector<mcmc>& sDist)
  {
    for (std::size_t i=0; i<sDist.size(); ++i)
      {
        particle[iParticle_end[b.getValue()[i]]].logLikNext=sDist[i].logLik;
      }

  }
  template<typename mcmc>
  void new_record(const Beta& b,const std::vector<mcmc>& sDist)
  {
    for (std::size_t i=0; i<sDist.size(); ++i)
      {
        Likelihood_Record_steps::Particle& p=particle[iParticle_end[b.getValue()[i]]];
        p.logLikInit=sDist[i].logLik;
        std::size_t j=p.ibeta.back();
        p.ibeta.clear();
        p.ibeta.push_back(j);
      }
    iParticle_init=iParticle_end;
  }

  void beta_jump(double betaless, double betamore)
  {
    auto iless=iParticle_end[betaless];
    auto imore=iParticle_end[betamore];
    auto ibetaless=Beta_to_i[betaless];
    auto ibetamore=Beta_to_i[betamore];
    particle[iless].ibeta.back()=ibetamore;
    particle[imore].ibeta.back()=ibetaless;
    iParticle_end[betaless]=imore;
    iParticle_end[betamore]=iless;
  }
  void push_step()
  {
    for (auto &p:particle)
      {
        std::size_t i=p.ibeta.back();
        p.ibeta.push_back(i);
      }
  }




};

struct Likelihood_Record{
  struct Particle
  {
    bool isValid=false;
    double logLikInit;
    double logLikNext;
  };

  std::size_t size()const {return particle.size();}
  std::vector<std::vector<Particle>> particle;
  void newParticles()
  {
    std::vector<Particle> p(desc_beta.size());
    particle.push_back(p);
  }

  std::vector<double> desc_beta;

};

template<class logLikelihood>
struct Gradient_Finite_Difference
{
  static std::string ClassName(){return "Gradient_Finite_Difference";}
  Gradient_Finite_Difference(const logLikelihood& l, double dx=1e-6):logL(l),dx_(dx){}
  const logLikelihood& logL;
  double dx_;
  M_Matrix<double> operator()(const M_Matrix<double>& x)const
  {
    M_Matrix<double> o(1,x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      {
        auto xp=x;
        auto xn=x;
        xp[i]+=dx_;
        xn[i]-=dx_;
        double Lp=logL(xp);
        double Ln=logL(xn);
        o[i]=(Lp-Ln)/dx_/2.0;

      }
    return o;
  }


};


template<class function>
struct Jacobian_Finite_Difference
{
  static std::string ClassName(){return "Jacobian_Finite_Difference";}
  Jacobian_Finite_Difference(const function& l, double dx=1e-6):fun(l),dx_(dx){}
  const function& fun;
  double dx_;
  M_Matrix<double> operator()(const M_Matrix<double>& x)const
  {
    auto y=fun(x);
    M_Matrix<double> o(y.size(),x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      {
        auto xp=x;
        auto xn=x;
        xp[i]+=dx_;
        xn[i]-=dx_;
        M_Matrix<double> Lp(y.size(),1,fun(xp));
        M_Matrix<double> Ln(y.size(),1,fun(xn));
        o(":",i,(Lp-Ln)*(0.5/dx_));

      }
    return o;
  }

};

template<class logLikelihood>
struct Hessian_Finite_Difference
{
  static std::string ClassName(){return "Hessian_Finite_Difference";}
  Hessian_Finite_Difference(const logLikelihood& l, double dx=1e-4):logL(l),dx_(dx){}

  const logLikelihood& logL;
  double dx_;
  M_Matrix<double> operator()(const M_Matrix<double>& x)const
  {
    double L0=logL(x);
    std::size_t n=x.size();
    M_Matrix<double> o(n,n);

    for (std::size_t i=0; i<n; ++i)
      {
        auto xp=x;
        auto xn=x;

        xp[i]+=dx_;
        xn[i]-=dx_;

        double Lp=logL(xp);
        double Ln=logL(xn);
        o(i,i)=(Lp-2*L0+Ln)/dx_/dx_;
        for (std::size_t j=i+1; j<n; ++j)
          {
            auto xpp=xp;
            auto xpn=xp;
            auto xnp=xn;
            auto xnn=xn;

            xpp[j]+=dx_;
            xnn[j]-=dx_;

            xnp[j]+=dx_;
            xpn[j]-=dx_;


            double Lpp=logL(xpp);
            double Lnn=logL(xnn);
            double Lpn=logL(xpn);
            double Lnp=logL(xnp);


            o(i,j)=(Lpp-Lnp-Lpn+Lnn)/(dx_*dx_*4.0);
            o(j,i)=o(i,j);

          }

      }
    return o;
  }
};

template<class logLikelihood>
Gradient_Finite_Difference<logLikelihood> make_Gradient_Finite_Difference(const logLikelihood& l, double dx=1e-7)
{
  return Gradient_Finite_Difference<logLikelihood>(l,dx);
}

template<class function>
Jacobian_Finite_Difference<function> make_Jacobian_Finite_Difference(const function& l, double dx=1e-8)
{
  return Jacobian_Finite_Difference<function>(l,dx);
}




template<class logLikelihood>
Hessian_Finite_Difference<logLikelihood> make_Hessian_Finite_Difference(const logLikelihood& l, double dx=1e-5)
{
  return Hessian_Finite_Difference<logLikelihood>(l,dx);
}




template<class F0, class F1>
bool test_against(const F0& f0, const F1& f1, const M_Matrix<double>& x, double resolution,double tolerance)
{
  std::cout<<"test  "<<F0::ClassName();
  std::cout<<"against  "<<F1::ClassName()<<"\n";
  std::cout<<"tested on x="<<x<<"\n";
  auto Gf=f0(x);
  auto Gd=f1(x);

  bool out=true;
  std::set<std::size_t> s;
  for (std::size_t i=0; i<Gf.size(); ++i)
    {
      std::cout<<i<<"\t"<<Gf[i]<<"\t"<<Gd[i]<<"\t";
      if (Gf[i]==0)
        {
          if (std::abs(Gd[i])>resolution)
            {
              std::cout<<Gd[i]<<"\trejected\n";
              out=false;
              s.insert(i);
            }
          else
            std::cout<<"accepted\n";

        }
      else
        {
          double r=(Gf[i]-Gd[i])/Gf[i];
          std::cout<<r;
          if (std::abs(r)>tolerance)
            {
              std::cout<<"\trejected\n";
              s.insert(i);
              out=false;
            }
          else
            std::cout<<"\taccepted\n";

        }

    }
  for (auto& i:s)
    {
      std::cout<<i<<"\t"<<Gf[i]<<"\t"<<Gd[i]<<"\t";
      if (Gf[i]==0)
        {
          std::cout<<Gd[i]<<"\trejected\n";
        }
      else
        {
          double r=(Gf[i]-Gd[i])/Gf[i];
          std::cout<<r;
          std::cout<<"\trejected\n";
        }

    }
  return out;
}




struct Master_Tempering_Likelihood
{
  struct Prior{
    std::vector<double> desc_beta_for_dL;
    std::vector<double> log_desc_beta_for_dL;

    double mL0;
    double sL0;

    double mmlogdL;
    double smlogdL;

    double mlogslogdL;
    double slogslogdL;

    double mloglandalogdL;
    double sloglandalogdL;

    double mlogepsilon2logdL;
    double slogepsilon2logdL;


    double mmlogdelta;
    double smlogdelta;

    double mlogslogdelta;
    double slogslogdelta;
    Prior(std::vector<double> _beta_for_dL,
          double _mL0,
          double _sL0,
          double _mmlogdL,
          double _smlogdL,
          double _mlogslogdL,
          double _slogslogdL,

          double _mloglandalogdL,
          double _sloglandalogdL,
          double _mlogepsilon2logdL,
          double _slogepsilon2logdL,
          double _mmlogdelta,
          double _smlogdelta,
          double _mlogslogdelta,
          double _slogslogdelta):
      desc_beta_for_dL(_beta_for_dL),log_desc_beta_for_dL(_beta_for_dL.size()),
      mL0(_mL0),sL0(_sL0),mmlogdL(_mmlogdL),smlogdL(_smlogdL),
      mlogslogdL(_mlogslogdL),slogslogdL(_slogslogdL),
      mloglandalogdL(_mloglandalogdL),sloglandalogdL(_sloglandalogdL),
      mlogepsilon2logdL(_mlogepsilon2logdL),slogepsilon2logdL(_slogepsilon2logdL),
      mmlogdelta(_mmlogdelta),smlogdelta(_smlogdelta), mlogslogdelta(_mlogslogdelta),slogslogdelta(_slogslogdelta){
      for(std::size_t i=0; i<log_desc_beta_for_dL.size(); ++i)
        log_desc_beta_for_dL[i]=std::log(desc_beta_for_dL[i]);

    }



  };

  class Lfit
  {
  public:

    Lfit(std::vector<double>betas):b(betas.rbegin(),betas.rend())
    {

      for (std::size_t i=0; i<b.size(); ++i)
        this->betas_[b[i]]=i;
    }

    double dL(double beta,const M_Matrix<double>& dLs)const
    {
      auto it0=betas_.find(beta);
      if (it0!=betas_.end())
        return dLs[it0->second];
      else
        {

          auto it=betas_.lower_bound(beta);
          auto it2=betas_.upper_bound(beta);
          double b1=it2->first;
          double b0=it->first;
          double dL1=dLs[it2->second];
          double dL0=dLs[it->second];

          double out=(beta-b0)/(b1-b0)*dL0+(-beta+b1)/(b1-b0)*dL1;
          return out;
        }

    }

    M_Matrix<double> ddL_ddL(double beta,const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> out(1,dLs.size(),0.0);
      auto it0=betas_.find(beta);
      if (it0!=betas_.end())
        {
          out(0,it0->second)=1;
          return out;
        }
      else
        {
          auto it=betas_.lower_bound(beta);
          auto it2=betas_.upper_bound(beta);
          double b1=it2->first;
          double b0=it->first;


          out(0,it->second)=(beta-b0)/(b1-b0);
          out(0,it2->second)=(-beta+b1)/(b1-b0);
          return out;
        }

    }




    M_Matrix<double> dL(const std::vector<double>& beta,const M_Matrix<double>& dLs)const
    {
      if (beta==b)
        return dLs;
      else {
          M_Matrix<double> out(1,beta.size());
          for(std::size_t i=0; i<beta.size(); ++i)
            out[i]=dL(beta[i],dLs);
          return out;
        }
    }

    M_Matrix<double> ddL_ddL(const std::vector<double>& beta,const M_Matrix<double>& dLs)const
    {
      if (beta==b)
        return eye<double>(dLs.size());
      else {
          M_Matrix<double> out(beta.size(),beta.size());
          for(std::size_t i=0; i<beta.size(); ++i)
            out(i,":",ddL_ddL(beta[i],dLs));
          return out;
        }
    }


    double L(double beta,const M_Matrix<double>& dLs,const M_Matrix<double>& Ls)const
    {
      auto it0=betas_.find(beta);
      if (it0!=betas_.end())
        return Ls[it0->second];
      else
        {

          auto it=betas_.lower_bound(beta);
          auto it2=betas_.upper_bound(beta);
          double b1=it2->first;
          double b0=it->first;
          double dL1=dLs[it2->second];
          double dL0=dLs[it->second];
          double out=Ls[it->first];
          out+=0.5*(beta-b0+(b1-beta)*(beta-b0)/(b1-b0))*dL0
              +0.5*sqr(beta-b0)/(b1-b0)*dL1;
          return out;
        }

    }




    M_Matrix<double> dL_ddL(double beta,const M_Matrix<double>& dLs_ddL)const
    {
      auto it0=betas_.find(beta);
      if (it0!=betas_.end())
        return dLs_ddL(it0->second,":");
      else
        {

          auto it=betas_.lower_bound(beta);
          auto it2=betas_.upper_bound(beta);
          double b1=it2->first;
          double b0=it->first;
          auto out=dLs_ddL(it->second,":");
          out(0,it->second)+=0.5*(beta-b0+(b1-beta)*(beta-b0)/(b1-b0));
          out(0,it2->second)+=0.5*sqr(beta-b0)/(b1-b0);

          return out;
        }

    }


    M_Matrix<double> L(const std::vector<double>& beta,double L0,const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> Ls=Lis(L0,dLs);
      if (beta==b)
        return Ls;
      else
        {

          M_Matrix<double> out(1,beta.size());
          for(std::size_t i=0; i<beta.size(); ++i)
            out[i]=L(beta[i],dLs,Ls);
          return out;
        }
    }


    M_Matrix<double> dL_ddL(const std::vector<double>& beta,const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> dLs_ddL=dLis_ddLj(dLs);
      if (beta==b)
        return dLs_ddL;
      else
        {

          M_Matrix<double> out(beta.size(),beta.size());
          for(std::size_t i=0; i<beta.size(); ++i)
            out(i,":",dL_ddL(beta[i],dLs_ddL));
          return out;
        }
    }



    M_Matrix<double> Lis(double L0,const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> out(1,dLs.size());
      double r=L0;
      double br=0;
      double dLr=dLs[0];
      for (std::size_t i=0; i<dLs.size(); ++i )
        {
          double db=b[i]-br;
          double dLv=(dLs[i]+dLr)*0.5;
          r+=db*(dLv);
          out[i]=r;
          br=b[i];
          dLr=dLs[i];

        }
      return out;
    }


    M_Matrix<double> dLis_ddLj(const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> out(dLs.size(),dLs.size(),0.0);

      double br=0;
      std::size_t ir=0;
      for (std::size_t i=0; i<dLs.size(); ++i )
        {
          double db=(b[i]-br);
          for (std::size_t j=i; j<dLs.size(); ++j)
            {
              out(j,i)+=db*0.5;
              out(j,ir)+=db*0.5;
            }
          br=b[i];
          ir=i;
        }
      return out;

    }



    M_Matrix<double> dLis_dL0dLj(const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> out(dLs.size(),1+dLs.size(),0.0);

      double br=0;
      std::size_t ir=0;
      for (std::size_t j=0; j<dLs.size(); ++j)
        {
          out(j,0)=1;
        }
      for (std::size_t i=0; i<dLs.size(); ++i )
        {
          double db=(b[i]-br);
          for (std::size_t j=i; j<dLs.size(); ++j)
            {
              out(j,1+i)+=db*0.5;
              out(j,1+ir)+=db*0.5;
            }
          br=b[i];
          ir=i;
        }


      return out;
    }



    M_Matrix<double> cov_L0dLj(double /*L0*/,const M_Matrix<double>& dLs,
                               const M_Matrix<double>& Hinv)const
    {
      std::size_t n=dLs.size();
      M_Matrix<double> cov(1+n, 1+n);
      cov(0,0)=Hinv(0,0);
      for (std::size_t i=0; i<n; ++i)
        {
          cov(0,1+i)=Hinv(0,7+i)*dLs[i];
          cov(1+i,0)=cov(0,1+i);
          for (std::size_t j=0; j<n; ++j)
            {
              cov(1+i,1+j)=Hinv(7+i,7+j)*dLs[i]*dLs[j];
            }
        }
      return cov;

    }

    M_Matrix<double> cov_Lis(const M_Matrix<double>& dLs,
                             const M_Matrix<double>& covL0dLj)const
    {
      auto b=dLis_dL0dLj(dLs);
      return b*multTransp(covL0dLj,b);
    }


    double E(double L0,const M_Matrix<double>& dLs)const
    {
      auto n=b.size()-1;
      double o=L0;
      o+=(
            +0.5*(b[n]-b[0])*(b[0])
          +1.0/3.0*sqr(b[0])
          +0.5*(b[n]-b[1+1])*(b[1+1]-b[1])
          +0.5*(b[n]-b[1])*(b[1]-b[1-1])
          +1.0/6.0*sqr(b[1]-b[1-1])
          +1.0/3.0*sqr(b[1+1]-b[1])
          )*dLs[0];
      for (std::size_t i=1; i<n; ++i)
        {
          o+=
              (+0.5*(b[n]-b[i+1])*(b[i+1]-b[i])
              +0.5*(b[n]-b[i])*(b[i]-b[i-1])
              +1.0/6.0*sqr(b[i]-b[i-1])
              +1.0/3.0*sqr(b[i+1]-b[i])
              )*dLs[i];
        }
      o+=(1.0/6.0*sqr(b[n]-b[n-1])
          )*dLs[n];
      auto Li=Lis(L0,dLs);
      double out;
      double rL=L0;
      double br=0;
      double dLr=dLs[0];
      for (std::size_t i=0; i<dLs.size(); ++i )
        {
          double db=b[i]-br;
          double dLv=(dLs[i]+2*dLr)/6.0;
          out+=rL*db+dLv*sqr(db);
          br=b[i];
          dLr=dLs[i];
          rL=Li[i];

        }
      // return out;

      return o;
    }


    M_Matrix<double>  dE_dL0dLj(const M_Matrix<double>& dLs)const
    {
      M_Matrix<double> out(1, 1+dLs.size());
      auto n=b.size()-1;
      out(0,0)=b[n];
      out(0,1)= (
            +0.5*(b[n]-b[0])*(b[0])
          +1.0/3.0*sqr(b[0])
          +0.5*(b[n]-b[1+1])*(b[1+1]-b[1])
          +0.5*(b[n]-b[1])*(b[1]-b[1-1])
          +1.0/6.0*sqr(b[1]-b[1-1])
          +1.0/3.0*sqr(b[1+1]-b[1])
          );
      for (std::size_t i=1; i<n; ++i)
        {
          out(0,i+1)=
              (+0.5*(b[n]-b[i+1])*(b[i+1]-b[i])
              +0.5*(b[n]-b[i])*(b[i]-b[i-1])
              +1.0/6.0*sqr(b[i]-b[i-1])
              +1.0/3.0*sqr(b[i+1]-b[i])
              );
        }
      out(0,n+1)=(1.0/6.0*sqr(b[n]-b[n-1])
          );
      return out;
    }


    double s2_E(const M_Matrix<double>& dLs, const M_Matrix<double>& covL0dLjs)const
    {
      auto b=dE_dL0dLj(dLs);
      return xTSigmaX(b,covL0dLjs);

    }


    std::vector<double> BetaOpt(const M_Matrix<double> dLs)const
    {
      std::vector<double> out;
      double b=1;
      while (b>0)
        {
          out.push_back(b);
          double d=dL(b,dLs);
          b-=sqrt(1.0/d);
        }
      return out;
    }


    std::vector<double> beta()const { return b;}

  private:
    std::map<double,std::size_t> betas_;
    std::vector<double> b; //ascending
  };

  struct Parameters{
    double L0;  // likelihood at beta=0;
    double mlogdL;  //
    double logvlogdL;
    double logrlandalogdL;
    double logepsilon2logL;

    double mlogdelta; // mean log of delta
    double logvlogdelta;

    M_Matrix<double> dL_j;
    // std::vector<std::vector<double>> delta_ij;
    M_Matrix<double> log_dL_j;
    // std::vector<std::vector<double>> log_delta_ij;

    // auxiliares
    double delta;
    double d0;
    double d1;
    double d4;

    double vlogdelta;
    Lfit lfit_;
    std::size_t size()const
    {
      return 7+dL_j.size();
    }


    double E()const
    {
      double e= lfit_.E(L0,dL_j);

      auto de=lfit_.dE_dL0dLj(dL_j);
      auto b=lfit_.beta();
      double br=0;
      double sum=L0;
      // double Lr=L0;
      auto Li=Ls();
      for (std::size_t i=0; i<Li.size(); ++i)
        {
          sum+=(Li[i])*(b[i]-br)/2.0;
          //Lr=Li[i];
          br=b[i];
        }
      return e;

    }

    M_Matrix<double> Ls()const
    {
      return lfit_.Lis(L0,dL_j);
    }



    static M_Matrix<double> G(
        double dL_dL0,
        double dL_dmlogdL,
        double dL_dslogdL,
        double dL_dlandalogdL,
        double dL_depsilonlogdL,
        double dL_dmlogdelta,
        double dL_dslogdelta,
        M_Matrix<double> dL_dlogdLj)
    {
      std::vector<double> o;
      o.push_back(dL_dL0);
      o.push_back(dL_dmlogdL);
      o.push_back(dL_dslogdL);
      o.push_back(dL_dlandalogdL);
      o.push_back(dL_depsilonlogdL);


      o.push_back(dL_dmlogdelta);
      o.push_back(dL_dslogdelta);
      std::vector<double> v=dL_dlogdLj.toVector();
      o.insert(o.end(),v.begin(),v.end());
      return M_Matrix<double>(1,o.size(),o);

    }


    static M_Matrix<double> H(
        double d2L_dL02,double d2L_dL0_dmdelta,double d2L_dL0_dsdelta,M_Matrix<double> d2L_dL0_dlogdL,

        double d2L_dmlogdL2,M_Matrix<double> d2L_dmlogdL2_ddL,
        double dLdlogsimga2,double dLdlogsigma2_dloglambda, double d2L_logsigma2_dlogepsilon2,
        double dL2dloglambda2,double d2L_logsigma2_logepsilon2,
        double d2L_dlogepsilon,
        double d2L_dmlogdelta2, double d2L_dmdelta_dsdelta, M_Matrix<double> d2L_dmdelta_ddL,
        double d2L_dslogdelta2, M_Matrix<double> d2L_dsdelta_ddL,
        M_Matrix<double> d2L_dlogdLj1_dlogdLj2)
    {
      std::size_t nL=d2L_dL0_dlogdL.size();

      std::size_t n=7+nL;
      M_Matrix<double> o(n,n,0.0);
      o(0,0)=d2L_dL02;
      o(0,5)=d2L_dL0_dmdelta;
      o(5,0)=d2L_dL0_dmdelta;
      o(0,6)=d2L_dL0_dsdelta;
      o(6,0)=d2L_dL0_dsdelta;
      for (std::size_t i=0; i<nL; ++i)
        {
          o(0,7+i)=d2L_dL0_dlogdL[i];
          o(7+i,0)=d2L_dL0_dlogdL[i];
        }
      o(1,1)= d2L_dmlogdL2;
      for (std::size_t i=0; i<nL; ++i)
        {
          o(1,7+i)=d2L_dmlogdL2_ddL[i];
          o(7+i,1)=d2L_dmlogdL2_ddL[i];
        }


      o(2,2)=dLdlogsimga2;
      o(2,3)=dLdlogsigma2_dloglambda;
      o(3,2)=o(2,3);
      o(2,4)=d2L_logsigma2_dlogepsilon2;
      o(4,2)=o(2,4);
      o(3,3)= dL2dloglambda2;
      o(3,4)=d2L_logsigma2_logepsilon2;
      o(4,3)=o(3,4);
      o(4,4)=d2L_dlogepsilon;
      o(5,5)=d2L_dmlogdelta2;
      o(5,6)=d2L_dmdelta_dsdelta;
      o(6,5)=d2L_dmdelta_dsdelta;
      for (std::size_t i=0; i<nL; ++i)
        {
          o(5,7+i)=d2L_dmdelta_ddL[i];
          o(7+i,5)=d2L_dmdelta_ddL[i];
        }

      o(6,6)=d2L_dslogdelta2;
      for (std::size_t i=0; i<nL; ++i)
        {
          o(6,7+i)=d2L_dsdelta_ddL[i];
          o(7+i,6)=d2L_dsdelta_ddL[i];
        }


      for (std::size_t i=0; i<nL; ++i)
        {
          for (std::size_t j=0; j<nL; ++j)
            {
              o(7+i,7+j)=d2L_dlogdLj1_dlogdLj2(i,j);
              o(7+j,7+i)=d2L_dlogdLj1_dlogdLj2(i,j);
            }

        }


      return o;
    }


    static M_Matrix<double> X(
        double L0,
        double mlogdL,  //
        double logvlogdL,
        double loglandalogdL,
        double logepsilon2logL,
        double mlogdelta,
        double logvlogdelta,
        M_Matrix<double> logdLj)
    {
      std::vector<double> o;
      o.push_back(L0);
      o.push_back(mlogdL);
      o.push_back(logvlogdL);
      o.push_back(loglandalogdL);
      o.push_back(logepsilon2logL);
      o.push_back(mlogdelta);
      o.push_back(logvlogdelta);
      std::vector<double> v=logdLj.toVector();
      o.insert(o.end(),v.begin(),v.end());
      return M_Matrix<double>(1,o.size(),o);

    }

    static M_Matrix<double> X(const std::vector<Likelihood_Record>& s)
    {
      Parameters p(s);
      return p();
    }

    static M_Matrix<double> X(const std::vector<Likelihood_Record>& s,
                              const Prior& pr,std::mt19937_64& mt)
    {
      Parameters p(s,pr,mt);
      return p();
    }


    struct getX
    {
      static
      double getL0(const Likelihood_Record& r)
      {
        std::size_t i=0;

        double sum=r.particle[0][i].logLikInit;
        for (i=1; i< r.particle[0].size(); ++i)
          {
            if (sum>r.particle[0][i].logLikInit)
              sum=r.particle[0][i].logLikInit;
          }
        return sum;
      }
      static
      M_Matrix<double> L_to_dL(const std::vector<double>& beta,double L0,const M_Matrix<double>& L)
      {
        std::size_t n=beta.size();

        M_Matrix<double> Lo=sort(L, std::greater<double>());
        if (Lo[n-1]<L0)
          std::swap(Lo[n-1],L0);
        M_Matrix<double> dL(1,beta.size(),0.0);
        double Lr=L0;
        double br=0;
        for (std::size_t i=0; i<beta.size(); ++i)
          {
            double dbeta=beta[n-i-1]-br;
            double dLr=Lo[n-i-1]-Lr;
            dL[i]=(dLr)/(dbeta);
            Lr=Lo[n-i-1]; br=beta[n-i-1];
          }
        return dL;
      }

      static
      std::pair<std::size_t, std::size_t> getIndex(const std::vector<Likelihood_Record>& r, std::size_t index)
      {
        std::size_t istep=0;
        std::size_t ipar=0;
        while(istep+r[ipar].size()<index)
          {
            istep+=r[ipar].size();
            ++ipar;
          }

        std::size_t m=index-istep;
        return {ipar,m};

      }
      static std::size_t getNumSteps(const std::vector<Likelihood_Record>& r)
      {
        std::size_t numsteps=0;
        for (auto& e:r)
          numsteps+=e.size();
        return numsteps;

      }

      static double getAcceptanceRatio(const std::vector<Likelihood_Record>& r)
      {
        std::size_t nA=0;
        std::size_t n=0;

        for (const Likelihood_Record& e:r)
          for (auto& pv:e.particle)
            for (auto& p:pv)
              {
                ++n;
                if (p.isValid) ++nA;
              }
        double ratio=1.0*nA/n;
        return ratio;

      }



      static
      M_Matrix<double> getLeq(const std::vector<Likelihood_Record>& r,  std::size_t nsteps,double delta, double L0)
      {
        double d=1-std::exp(-delta*nsteps);
        std::size_t n=r.back().particle.back().size();
        M_Matrix<double> Lf (1,n);
        for (std::size_t i=0;i< n; ++i)
          {
            double L1=r.back().particle.back()[i].logLikInit;
            double dL=L0-L1;
            double LL=L0-dL/d;

            Lf[i]=LL;
          }
        return Lf;
      }

      static
      std::vector<std::vector<double>> get_logdelta_ij(
          const std::vector<Likelihood_Record>& r,
          double L0,
          const M_Matrix<double>& dL)
      {
        Lfit Lf(r.back().desc_beta);

        std::vector<std::vector<double>> out;
        for (std::size_t i=0; i< r.size(); ++i)
          {
            auto Ln=Lf.L(r[i].desc_beta,L0,dL);
            for (std::size_t ii=0; ii<r[i].particle.size(); ++ii)
              {
                std::vector<double> deltas_j;

                for (std::size_t j=0; j<r[i].desc_beta.size(); ++j)

                  {
                    if (r[i].particle[ii][j].isValid)
                      {
                        double DL=+r[i].particle[ii][j].logLikNext
                            -r[i].particle[ii][j].logLikInit;
                        double DE=Ln[j]-r[i].particle[ii][j].logLikNext;
                        double deldd=std::abs(DL/DE);
                        deltas_j.push_back(std::log(deldd));
                      }
                  }
                out.push_back(deltas_j);
              }
          }
        return out;
      }

      static
      double get_logdelta(const std::vector<Likelihood_Record>& r)
      {
        const double mindelta=0.05;
        std::size_t numsteps=getNumSteps(r);
        double pACC=getAcceptanceRatio(r);
        double L0=0;
        double L1=0;
        double L2=0;
        for (std::size_t i=0; i<r[0].particle[0].size(); ++i)
          L0+=r[0].particle[0][i].logLikInit;
        L0/=r[0].particle[0].size();
        for (std::size_t i=0; i<r.back().particle.back().size(); ++i)
          L2+=r.back().particle.back()[i].logLikInit;
        L2/=r.back().particle.back().size();
        auto in=getIndex(r,numsteps/2);
        for (std::size_t i=0; i<r[in.first].particle[in.second].size(); ++i)
          L1+=r[in.first].particle[in.second][i].logLikInit;
        L1/=r[in.first].particle[in.second].size();
        double dL0=L0-L1;
        double dL1=L0-L2;
        double del=dL1/dL0-1;
        double delta;
        if (del<1)
          delta=-log(del)/numsteps*2*pACC;
        else
          delta=mindelta/numsteps*2*pACC;
        return log(delta);
      }


    };

    M_Matrix<double> operator()()const

    {
      return X(L0,mlogdL,logvlogdL,logrlandalogdL,logepsilon2logL,mlogdelta,logvlogdelta,log_dL_j);
    }

    Parameters(const std::vector<Likelihood_Record>& s):
      lfit_(s.back().desc_beta)
    {
      L0=getX::getL0(s[0]);
      mlogdelta=getX::get_logdelta(s);
      logvlogdelta=std::log(std::abs(mlogdelta))*0.5;
      std::size_t nsteps=getX::getNumSteps(s);
      auto Lq=getX::getLeq(s,nsteps,exp(mlogdelta),L0);
      dL_j=getX::L_to_dL(s.back().desc_beta,L0,Lq);
      this->logvlogdL=1;
      this->logrlandalogdL=-1;
      logepsilon2logL=-5;
      log_dL_j=dL_j.apply([](double x){return std::log(x);});
      this->mlogdL=mean(log_dL_j);
      delta=std::exp(mlogdelta);
      d0=1.0/(1.0+delta);
      d1=delta/(1.0+delta);
      d4=sqr(delta/sqr(1+delta));
      vlogdelta=std::exp(logvlogdelta);


    }

    Parameters(const std::vector<Likelihood_Record>& s,const Prior& p, std::mt19937_64& mt):
      lfit_(s.back().desc_beta)
    {
      using normal=
      std::normal_distribution<double>;
      L0=normal(p.mL0,p.sL0)(mt);
      mlogdelta=normal(p.mmlogdelta,p.smlogdelta)(mt);
      logvlogdelta=normal(p.mlogslogdelta, p.slogslogdelta)(mt);
      std::size_t nsteps=getX::getNumSteps(s);
      auto Lq=getX::getLeq(s,nsteps,exp(mlogdelta),L0);
      dL_j=getX::L_to_dL(s.back().desc_beta,L0,Lq);
      this->logvlogdL=normal(p.mlogslogdL, p.slogslogdL)(mt);
      this->logrlandalogdL=normal(p.mloglandalogdL, p.sloglandalogdL)(mt);
      logepsilon2logL=normal(p.mlogepsilon2logdL, p.slogepsilon2logdL)(mt);
      log_dL_j=dL_j.apply([](double x){return std::log(x);});
      this->mlogdL=normal(p.mmlogdL,p.smlogdL)(mt);
      delta=std::exp(mlogdelta);
      d0=1.0/(1.0+delta);
      d1=delta/(1.0+delta);
      d4=sqr(delta/sqr(1+delta));
      vlogdelta=std::exp(logvlogdelta);


    }

    Parameters(const std::vector<Likelihood_Record>& s,const Prior& p, std::mt19937_64& mt, bool ):
      lfit_(s.back().desc_beta)
    {
      using normal=
      std::normal_distribution<double>;
      L0=normal(p.mL0,p.sL0)(mt);
      mlogdelta=normal(p.mmlogdelta,p.smlogdelta)(mt);
      logvlogdelta=normal(p.mlogslogdelta, p.slogslogdelta)(mt);
      std::size_t nsteps=getX::getNumSteps(s);
      auto Lq=getX::getLeq(s,nsteps,exp(mlogdelta),L0);
      dL_j=getX::L_to_dL(s.back().desc_beta,L0,Lq);
      logvlogdL=normal(p.mlogslogdL, p.slogslogdL)(mt);
      logrlandalogdL=normal(p.mloglandalogdL, p.sloglandalogdL)(mt);
      logepsilon2logL=normal(p.mlogepsilon2logdL, p.slogepsilon2logdL)(mt);
      log_dL_j=dL_j.apply([](double x){return std::log(x);});
      this->mlogdL=normal(p.mmlogdL,p.smlogdL)(mt);


    }


    Parameters(const M_Matrix<double>& param,
               const std::vector<Likelihood_Record>&/* s*/,
               std::vector<double> betafordL)
      :
        dL_j{1,betafordL.size()}, log_dL_j(1,betafordL.size()),lfit_(betafordL)

    {
      std::size_t i=0;
      L0=param[i];
      ++i;
      mlogdL=param[i];
      ++i;
      logvlogdL=param[i];
      ++i;
      logrlandalogdL=param[i];
      ++i;
      logepsilon2logL=param[i];
      ++i;

      mlogdelta=param[i];
      ++i;
      logvlogdelta=param[i];
      ++i;
      for (std::size_t j=0; j<dL_j.size();     ++j)
        {
          log_dL_j[j]=param[i];
          dL_j[j]=exp(param[i]);
          ++i;
        }
      delta=std::exp(mlogdelta);
      d0=1.0/(1.0+delta);
      d1=delta/(1+delta);
      d4=sqr(delta/sqr(1+delta));
      vlogdelta=std::exp(logvlogdelta);


    }




  };


  struct Parameters_SE: public Parameters
  {

    Parameters_SE(const M_Matrix<double>& param,
                  const std::vector<Likelihood_Record>& s,
                  std::vector<double> betafordL,
                  const M_Matrix<double> _Hinv):
      Parameters(param,s,betafordL),
      CovHinv(-_Hinv),
      covL0dLj(lfit_.cov_L0dLj(L0,dL_j,-_Hinv)){}



    double L0_er()const {
      return std::exp(std::sqrt(CovHinv(0,0)))/std::abs(L0);
    }
    double L0_se()const {
      return std::sqrt(CovHinv(0,0));
    }

    double mlogdL_se()const
    {
      return std::sqrt(CovHinv(1,1));

    }
    double logvdL_se()const
    {
      return std::sqrt(CovHinv(2,2));

    }
    double loglandadL_se()const
    {
      return std::sqrt(CovHinv(3,3));

    }
    double logepsilondL_se()const
    {
      return std::sqrt(CovHinv(4,4));

    }

    double mlogdelta_se() const {
      return std::sqrt(CovHinv(5,5));
    }
    double mlogdelta_re() const {
      return std::sqrt(CovHinv(5,5))/std::abs(mlogdelta);
    }



    double logvlogdelta_se()const {
      return std::sqrt(CovHinv(6,6));
    }

    double vlogdelta_re()const {
      return std::exp(std::sqrt(CovHinv(6,6)))-1.0;
    }


    double vlogdelta_se()const {
      return std::sqrt(CovHinv(6,6))*std::exp(logvlogdelta);
    }



    M_Matrix<double> logdL_j_se()const {
      M_Matrix<double> o(1,log_dL_j.size());
      for (std::size_t i=0; i<o.size();++i)
        o[i]=std::sqrt(CovHinv(7+i,7+i));
      return o;
    }

    M_Matrix<double> dL_j_se()const {
      M_Matrix<double> o(1,log_dL_j.size());
      for (std::size_t i=0; i<o.size();++i)
        o[i]=std::sqrt(CovHinv(7+i,7+i)*sqr(dL_j[i]));
      return o;
    }



    M_Matrix<double> logdL_j_er()const {
      M_Matrix<double> o(1,log_dL_j.size());
      for (std::size_t i=0; i<o.size();++i)
        o[i]=std::exp(CovHinv(7+i,7+i)/2)-1.0;
      return o;

    }




    M_Matrix<double> dL_j_er()const {
      M_Matrix<double> o(1,dL_j.size());
      for (std::size_t i=0; i<o.size();++i)
        o[i]=exp(std::sqrt(CovHinv(7+i,7+i)))-1.0;
      return o;
    }

    M_Matrix<double> L_j_se()const
    {
      auto cov=lfit_.cov_Lis(dL_j,covL0dLj);
      M_Matrix<double> out(1,dL_j.size());
      for (std::size_t i=0; i< out.size(); ++i)
        out[i]=std::sqrt(cov(i,i));
      return out;
    }

    double E_se()const
    {
      return std::sqrt(lfit_.s2_E(dL_j,covL0dLj));
    }



    M_Matrix<double> CovHinv;
    M_Matrix<double> covL0dLj;

  };






  struct logL{
    const Beta& beta;
    const std::vector<Likelihood_Record>& s;
    const Prior& p;


    logL(const Beta& beta_,
         const std::vector<Likelihood_Record>& s_,
         const Prior& p_):beta(beta_),s(s_),p(p_){}


    double operator ()(const M_Matrix<double> parameters)const
    {
      Parameters par(parameters,s,p.desc_beta_for_dL);
      return loglik(par,s,p);

    }


    static double loglik(const Parameters& par,
                         const std::vector<Likelihood_Record>& v,
                         const Prior& p)
    {
      double sumlogL_L=0;

      for (std::size_t ii=0; ii<v.size(); ++ii)
        {
          const Likelihood_Record &s=v[ii];

          Lfit lf(p.desc_beta_for_dL);
          M_Matrix<double> Leq=lf.L(s.desc_beta,par.L0,par.dL_j);
          M_Matrix<double> dLi=lf.dL(s.desc_beta,par.dL_j);


          for (std::size_t i=0; i<s.particle.size(); ++i)
            {

              for (std::size_t j=0; j<s.particle[i].size(); ++j)
                {

                  const Likelihood_Record::Particle& p=s.particle[i][j];
                  if (p.isValid)
                    {
                      double L=p.logLikNext;
                      double Linit=p.logLikInit;
                      double Lfinal=Leq[j];
                      double dL=dLi[j];
                      double mL=par.d0*Linit+par.d1*Lfinal;
                      double varianceL=dL+sqr(Lfinal-Linit)*par.d4*par.vlogdelta;
                      double logL_L=Normal(L,mL,varianceL,true);
                      sumlogL_L+=logL_L;
                    }
                }

            }

        }

      double logL_L0=Normal(par.L0,p.mL0,p.sL0);
      double logL_dL=univariate_gaussian_process_distribution::logL
          (par.mlogdL,std::exp(par.logvlogdL),std::exp(par.logrlandalogdL),
           std::exp(par.logepsilon2logL),
           par.log_dL_j,p.log_desc_beta_for_dL);
      double logL_mlogdL=Normal(par.mlogdL,p.mmlogdL,p.smlogdL);
      double logL_slogdL=Normal(par.logvlogdL,p.mlogslogdL,p.slogslogdL);
      double logL_landalogdL=Normal
          (par.logrlandalogdL,p.mloglandalogdL,p.sloglandalogdL);
      double logL_elogdL=Normal(par.logepsilon2logL,p.mlogepsilon2logdL,p.slogepsilon2logdL);

      double logL_mdelta=Normal(par.mlogdelta,p.mmlogdelta,p.smlogdelta);
      double logL_sdelta=Normal(par.logvlogdelta, p.mlogslogdelta,p.slogslogdelta);

      double out= logL_L0+sumlogL_L+logL_dL+logL_mdelta+logL_sdelta+logL_mlogdL+
          logL_slogdL+logL_landalogdL+logL_elogdL;
      if (std::isfinite(out)&&out>1e8)
        return out;
      else
        return out;
    }
  };


  struct G{
    static std::string ClassName(){return "Master_Tempering_Gradient";}

    const logL& logL_;
    const std::vector<Likelihood_Record>& v;
    const Prior& p ;

    G(const logL& L): logL_(L),v(L.s), p(L.p){}


    bool test(const M_Matrix<double>& x, double dx=1e-7, double tol=1e-4)
    {

      auto Gfd=make_Gradient_Finite_Difference(logL_,dx);
      return test_against(*this,Gfd,x,dx,tol);

    }


    M_Matrix<double> operator()(const M_Matrix<double>& x)const
    {
      Parameters par(x,v,p.desc_beta_for_dL);
      double dL_dL0=0;
      double dL_dmlogdelta=0;
      double dL_dlogvlogdelta=0;
      M_Matrix<double> dL_ddL(1, par.dL_j.size(),0.0);


      for (std::size_t ii=0; ii<v.size(); ++ii)
        {
          const Likelihood_Record &s=v[ii];

          Lfit lf(p.desc_beta_for_dL);
          M_Matrix<double> Leq=lf.L(s.desc_beta,par.L0,par.dL_j);

          M_Matrix<double> dLi=lf.dL(s.desc_beta,par.dL_j);

          M_Matrix<double> dLeq_ddL=lf.dL_ddL(s.desc_beta,par.dL_j);

          M_Matrix<double> ddLi_ddL=lf.ddL_ddL(s.desc_beta,par.dL_j);


          for (std::size_t i=0; i<s.particle.size(); ++i)
            {

              for (std::size_t j=0; j<s.particle[i].size(); ++j)
                {

                  const Likelihood_Record::Particle& part=s.particle[i][j];
                  if (part.isValid)
                    {
                      double L=part.logLikNext;
                      double Li=part.logLikInit;
                      double Lf=Leq[j];
                      double dL=dLi[j];
                      double mL=par.d0*Li+par.d1*Lf;
                      double vL=dL+sqr(Lf-Li)*par.d4*par.vlogdelta;
                      double dL_dm=(L-mL)/vL;
                      double dL_dv=0.5*(sqr(L-mL)/vL-1);
                      dL_dL0+= dL_dm*par.d1
                          +dL_dv/vL*(2.0*(Lf-Li)*par.d4*par.vlogdelta);


                      for (std::size_t jj=0; jj<par.dL_j.size(); ++jj)
                        {
                          dL_ddL(0,jj)+= dL_dm*par.d1*dLeq_ddL(j,jj)
                              +dL_dv/vL*(ddLi_ddL(j,jj)+dLeq_ddL(j,jj)*
                                         2.0*(Lf-Li)*par.d4*par.vlogdelta);
                        }
                      dL_dmlogdelta+=dL_dm*(Lf-Li)*par.d1/(1+par.delta)
                          +dL_dv/vL*sqr(Lf-Li)*par.d4*par.vlogdelta*2*
                          (1-par.delta)/(1+par.delta);

                      dL_dlogvlogdelta+=dL_dv/vL*(sqr(Lf-Li)*par.d4*par.vlogdelta);


                    }
                }

            }
        }
      for (std::size_t jj=0; jj<par.dL_j.size(); ++jj)
        {
          dL_ddL(0,jj)*=par.dL_j[jj];
        }

      dL_dL0+=-(par.L0-p.mL0)/sqr(p.sL0);
      dL_dmlogdelta-=(par.mlogdelta-p.mmlogdelta)/sqr(p.smlogdelta);
      dL_dlogvlogdelta-=(par.logvlogdelta-p.mlogslogdelta)/sqr(p.slogslogdelta);
      M_Matrix<double> m=par.log_dL_j;
      m=par.mlogdL;

      auto cov=univariate_gaussian_process_distribution::cov
          (p.log_desc_beta_for_dL,std::exp(par.logvlogdL),
           std::exp(par.logrlandalogdL),std::exp(par.logepsilon2logL));
      auto covinv=invSafe(cov);

      auto Ggp=univariate_gaussian_process_distribution::G
          (par.log_dL_j,m,p.log_desc_beta_for_dL,cov,covinv,
           std::exp(par.logvlogdL),std::exp(par.logrlandalogdL), std::exp(par.logepsilon2logL));




      double dL_dmlogdL=0;
      std::size_t i;
      for ( i=0; i<par.log_dL_j.size(); ++i)
        {
          dL_ddL(0,i)-=Ggp[i];
          dL_dmlogdL+=Ggp[i];
        }
      dL_dmlogdL-=(par.mlogdL-p.mmlogdL)/sqr(p.smlogdL);
      double dL_dslogdL=Ggp[i]-(par.logvlogdL-p.mlogslogdL)/sqr(p.slogslogdL);
      ++i;
      double dL_dlandalogdL=Ggp[i]-(par.logrlandalogdL-p.mloglandalogdL)/sqr(p.sloglandalogdL);
      ++i;
      double dL_depsilonlogdL=Ggp[i]-(par.logepsilon2logL-p.mlogepsilon2logdL)/sqr(p.slogepsilon2logdL);
      ++i;

      return Parameters::G(dL_dL0,dL_dmlogdL,dL_dslogdL,dL_dlandalogdL,dL_depsilonlogdL,dL_dmlogdelta,dL_dlogvlogdelta,dL_ddL);
    }
  };


  struct H{
    static std::string ClassName(){return "Master_Tempering_Hessian";}

    const logL& logL_;

    const std::vector<Likelihood_Record>& v;
    const Prior& p;

    H(const logL& L): logL_(L),v(L.s), p(L.p){}

    bool test(const M_Matrix<double>& /*x*/,double/* res=1e-7*/,double/* tol=1e-2*/)
    {
      /*
      auto Hfd=GN_H(logL_);
      return test_against(*this,Hfd,x,res,tol);
*/
      return false;
    }




    M_Matrix<double> operator()(const M_Matrix<double>& x) const
    {
      Parameters par(x,v,p.desc_beta_for_dL);
      std::size_t nL=par.dL_j.size();


      double d2L_dL02=0;
      double d2L_dL0_dmdelta=0;
      double d2L_dL0_dsdelta=0;
      M_Matrix<double> d2L_dL0_ddL(1,nL,0.0);

      double d2L_dmlogdL2=0;


      double d2L_dmdelta2=0;
      double d2L_dmdelta_dsdelta=0;
      M_Matrix<double> d2L_dmdelta_ddL(1,nL, 0.0);

      double d2L_dsdelta2=0;
      M_Matrix<double> d2L_dsdelta_ddL(1,nL, 0.0);

      M_Matrix<double> d2L_ddL_ddL(nL, nL,0.0);

      for (std::size_t ii=0; ii<v.size(); ++ii)
        {
          const Likelihood_Record &s=v[ii];

          Lfit lf(p.desc_beta_for_dL);
          M_Matrix<double> Leq=lf.L(s.desc_beta,par.L0,par.dL_j);

          M_Matrix<double> dLi=lf.dL(s.desc_beta,par.dL_j);

          M_Matrix<double> dLeq_ddL=lf.dL_ddL(s.desc_beta,par.dL_j);

          M_Matrix<double> ddLi_ddL=lf.ddL_ddL(s.desc_beta,par.dL_j);


          for (std::size_t i=0; i<s.particle.size(); ++i)
            {

              for (std::size_t j=0; j<s.particle[i].size(); ++j)
                {

                  const Likelihood_Record::Particle& part=s.particle[i][j];
                  if (part.isValid)
                    {
                      double Li=part.logLikInit;
                      double Lf=Leq[j];
                      double dL=dLi[j];
                      double vL=dL+sqr(Lf-Li)*par.d4*par.vlogdelta;

                      for (std::size_t jj=0; jj<par.dL_j.size(); ++jj)
                        {
                          for (std::size_t jjj=0; jjj<par.dL_j.size(); ++jjj)
                            {
                              d2L_ddL_ddL(jj,jjj)+=
                                  (-1.0/vL*sqr(par.d1)*dLeq_ddL(j,jj)*dLeq_ddL(j,jjj)
                                   -0.5/vL*(ddLi_ddL(j,jj)+dLeq_ddL(j,jj)*
                                            2.0*(Lf-Li)*par.d4*par.vlogdelta)
                                   *1.0/vL*(ddLi_ddL(j,jjj)+dLeq_ddL(j,jjj)*
                                            2.0*(Lf-Li)*par.d4*par.vlogdelta))*
                                  par.dL_j[jj]*par.dL_j[jjj]
                                  ;

                            }
                          d2L_dL0_ddL(0,jj)+=
                              (-1.0/vL*sqr(par.d1)*dLeq_ddL(j,jj)
                               -0.5/vL*(ddLi_ddL(j,jj)+dLeq_ddL(j,jj)*
                                        2.0*(Lf-Li)*par.d4*par.vlogdelta)
                               *1.0/vL*(2.0*(Lf-Li)*par.d4*par.vlogdelta))*
                              par.dL_j[jj];
                          ;
                          d2L_dmdelta_ddL(0,jj)+=
                              (-1.0/vL*(Lf-Li)*par.d1/(1+par.delta)
                               *par.d1*dLeq_ddL(j,jj)
                               -0.5/vL*sqr(Lf-Li)*par.d4*par.vlogdelta
                               *2*(1-par.delta)/(1+par.delta)
                               *1.0/vL*(ddLi_ddL(j,jj)+dLeq_ddL(j,jj)
                                        *2.0*(Lf-Li)*par.d4*par.vlogdelta))*
                              par.dL_j[jj];

                          ;
                          d2L_dsdelta_ddL(0,jj)+=
                              (-0.5/vL*sqr(Lf-Li)*par.d4*par.vlogdelta
                               *1.0/vL*(ddLi_ddL(j,jj)+dLeq_ddL(j,jj)
                                        *2.0*(Lf-Li)*par.d4*par.vlogdelta))*
                              par.dL_j[jj];

                        }
                      d2L_dL02+=
                          -1.0/vL*sqr(par.d1)
                          -0.5*sqr(1.0/vL*(2.0*(Lf-Li)*par.d4*par.vlogdelta))
                          ;

                      d2L_dL0_dmdelta+=
                          -1.0/vL*(Lf-Li)*par.d1/(1+par.delta)
                          *par.d1
                          -0.5/vL*sqr(Lf-Li)*par.d4*par.vlogdelta
                          *2*(1-par.delta)/(1+par.delta)
                          *1.0/vL*(2.0*(Lf-Li)*par.d4*par.vlogdelta)

                          ;
                      d2L_dL0_dsdelta+=
                          -0.5/vL*sqr(Lf-Li)*par.d4*par.vlogdelta
                          *1.0/vL*(2.0*(Lf-Li)*par.d4*par.vlogdelta);


                      d2L_dmdelta2+=
                          -1.0/vL*sqr((Lf-Li)*par.d1/(1+par.delta))
                          -0.5*sqr(1.0/vL*sqr(Lf-Li)*par.d4*par.vlogdelta
                                   *2*(1-par.delta)/(1+par.delta))

                          ;
                      d2L_dsdelta2+=
                          -0.5*sqr(1.0/vL*sqr(Lf-Li)*par.d4*par.vlogdelta);

                      d2L_dmdelta_dsdelta+=
                          -0.5*sqr(1.0/vL*sqr(Lf-Li)*par.d4*par.vlogdelta)
                          *2*(1-par.delta)/(1+par.delta)

                          ;


                    }
                }

            }
        }
      d2L_dmdelta2+= -1.0/sqr(p.smlogdelta);
      d2L_dsdelta2+= -1.0/sqr(p.slogslogdelta);


      auto cov=univariate_gaussian_process_distribution::cov
          (p.log_desc_beta_for_dL,std::exp(par.logvlogdL),
           std::exp(par.logrlandalogdL), std::exp(par.logepsilon2logL));
      auto covinv=invSafe(cov);

      auto Hgp=univariate_gaussian_process_distribution::H
          (p.log_desc_beta_for_dL,cov,covinv,
           std::exp(par.logvlogdL),std::exp(par.logrlandalogdL), std::exp(par.logepsilon2logL));

      std::size_t i;
      M_Matrix<double> d2L_dmlogdL2_ddL(1,nL,0.0);
      for ( i=0; i<par.dL_j.size(); ++i)
        for (std::size_t j=0; j<par.dL_j.size(); ++j)
          {
            d2L_ddL_ddL(i,j)+=Hgp(i,j);
            d2L_dmlogdL2+=Hgp(i,j);
            d2L_dmlogdL2_ddL(0,i)-=Hgp(i,j);
          }

      d2L_dmlogdL2-=1.0/sqr(p.smlogdL);
      double d2L_dlogsimga2=Hgp(i,i);
      d2L_dlogsimga2-=1.0/sqr(p.slogslogdL);


      double dLdlogsigma2_dlogrlambda=Hgp(i,i+1);
      double d2L_dlogsigma2_dlogepsilon=Hgp(i,i+2);
      ++i;

      double d2L_dloglambda2=Hgp(i,i);
      d2L_dloglambda2-=1.0/sqr(p.sloglandalogdL);
      double  d2L_dloglambda_dlogepsilon=Hgp(i,i+1);

      d2L_dL02+=-1.0/sqr(p.sL0);
      ++i;
      double d2L_dlogepsilon=Hgp(i,i);
      d2L_dlogepsilon-=1.0/sqr(p.slogepsilon2logdL);



      return Parameters::H(d2L_dL02,d2L_dL0_dmdelta,d2L_dL0_dsdelta,d2L_dL0_ddL,
                           d2L_dmlogdL2,d2L_dmlogdL2_ddL,
                           d2L_dlogsimga2,dLdlogsigma2_dlogrlambda,d2L_dlogsigma2_dlogepsilon,

                           d2L_dloglambda2,d2L_dloglambda_dlogepsilon,
                           d2L_dlogepsilon,
                           d2L_dmdelta2, d2L_dmdelta_dsdelta,d2L_dmdelta_ddL,
                           d2L_dsdelta2,d2L_dsdelta_ddL,
                           d2L_ddL_ddL);
    }
  };
  /*
  struct GN_H{
    struct f{
      const std::vector<Likelihood_Record>& v;
      const Prior& p;


      f(const std::vector<Likelihood_Record>& s_,
        const Prior& p_):v(s_),p(p_){}


      std::vector<double> operator ()(const M_Matrix<double> parameters)const
      {
        Parameters par(parameters,v,p.desc_beta_for_dL);


        std::vector<double> out;
        std::size_t isample=0;

        for (std::size_t ii=0; ii<v.size(); ++ii)
          {
            const Likelihood_Record &s=v[ii];

            Lfit lf(p.desc_beta_for_dL);
            M_Matrix<double> Leq=lf.L(s.desc_beta,par.L0,par.dL_j);
            M_Matrix<double> dLi=lf.dL(s.desc_beta,par.dL_j);


            for (std::size_t i=0; i<s.particle.size(); ++i)
              {
                std::size_t idelta=0;

                for (std::size_t j=0; j<s.particle[i].size(); ++j)
                  {

                    const Likelihood_Record::Particle& p=s.particle[i][j];
                    if (p.isValid)
                      {
                        //double L=p.logLikNext;
                        double Linit=p.logLikInit;
                        double Lfinal=Leq[j];
                        double dL=dLi[j];
                        // double logL_L=Normal(L,mL,varianceL,true);
                     //   out.push_back(mL);
                     //   out.push_back(std::log(varianceL));
                        // out+=logL_L;
                      }
                  }

              }

          }

        //double logL_L0=Normal(par.L0,p.mL0,p.sL0);
        out.push_back(par.L0);
        out.push_back(std::log(p.sL0));
        //double logL_dL=univariate_gaussian_process_distribution::logL
        //    (par.mlogdL,std::exp(par.logvlogdL),std::exp(par.loglandalogdL),
        //     par.log_dL_j,p.log_desc_beta_for_dL);
        for (std::size_t i5=0; i5<par.log_dL_j.size(); ++i5)
          for (std::size_t j5=0; j5<par.log_dL_j.size(); ++j5)
            {
              if (i5==j5)
                out.push_back(par.log_dL_j[i5]-par.mlogdL);
              else
                out.push_back(par.log_dL_j[i5]-par.mlogdL+par.log_dL_j[j5]);
            }
        out.push_back(par.logvlogdL);
        //out.push_back(par.logvlogdL+par.loglandalogdL);
        out.push_back(par.loglandalogdL);

        // double logL_mlogdL=Normal(par.mlogdL,p.mmlogdL,p.smlogdL);
        out.push_back(par.mlogdL);
        out.push_back(std::log(sqr(p.smlogdL)));

        //double logL_slogdL=Normal(par.logvlogdL,p.mlogslogdL,p.slogslogdL);
        out.push_back(par.logvlogdL);
        out.push_back(std::log(sqr(p.slogslogdL)));


        //double logL_landalogdL=Normal
        //     (par.loglandalogdL,p.mloglandalogdL,p.sloglandalogdL);
        out.push_back(par.loglandalogdL);
        out.push_back(std::log(sqr(p.sloglandalogdL)));


        //double logL_mdelta=Normal(par.mlogdelta,p.mmlogdelta,p.smlogdelta);
        out.push_back(par.mlogdelta);
        out.push_back(std::log(sqr(p.smlogdelta)));


        //  double logL_sdelta=Normal(par.logvlogdelta, p.mlogslogdelta,p.slogslogdelta);
        out.push_back(par.logvlogdelta);
        out.push_back(std::log(sqr(p.slogslogdelta)));


        // return logL_L0+out+logL_dL+logL_mdelta+logL_sdelta+logL_mlogdL+
        //     logL_slogdL+logL_landalogdL;

        return out;
      }




    };

    struct d2{
      const std::vector<Likelihood_Record>& v;
      const Prior& p;


      d2(const std::vector<Likelihood_Record>& s_,
         const Prior& p_):v(s_),p(p_){}


      std::vector<double> operator ()(const M_Matrix<double> parameters)const
      {
        Parameters par(parameters,v,p.desc_beta_for_dL);


        std::vector<double> out;
        std::size_t isample=0;

        for (std::size_t ii=0; ii<v.size(); ++ii)
          {
            const Likelihood_Record &s=v[ii];

            Lfit lf(p.desc_beta_for_dL);
            M_Matrix<double> Leq=lf.L(s.desc_beta,par.L0,par.dL_j);
            M_Matrix<double> dLi=lf.dL(s.desc_beta,par.dL_j);


            for (std::size_t i=0; i<s.particle.size(); ++i)
              {
                std::size_t idelta=0;

                for (std::size_t j=0; j<s.particle[i].size(); ++j)
                  {

                    const Likelihood_Record::Particle& p=s.particle[i][j];
                    if (p.isValid)
                      {
                        //double L=p.logLikNext;
                        //double Linit=p.logLikInit;
                        //double Lfinal=Leq[j];
                        double dL=dLi[j];
                        double delta=par.delta_ij[isample][idelta];
                        // double mL=Linit/(1+delta)+delta/(1+delta)*Lfinal;
                        double varianceL=delta/(delta+1)*dL;
                        // double logL_L=Normal(L,mL,varianceL,true);
                        out.push_back(-1.0/varianceL);
                        out.push_back(-0.5);
                        // out+=logL_L;
                        ++idelta;
                      }
                  }
                //double logL_deltas=
                //  Normal(par.log_delta_ij[isample],par.mlogdelta,exp(par.logvlogdelta),true);
                for (std::size_t i4=0; i4<par.log_delta_ij[isample].size(); ++i4)
                  {
                    out.push_back(-1.0/exp(par.logvlogdelta));
                    out.push_back(-0.5);
                  }
                //out+=logL_deltas;
                ++isample;

              }

          }

        //double logL_L0=Normal(par.L0,p.mL0,p.sL0);
        out.push_back(-1.0/sqr(p.sL0));
        out.push_back(-0.5);
        auto rcov=univariate_gaussian_process_distribution::rcov
            (p.log_desc_beta_for_dL,
             std::exp(par.loglandalogdL));
        auto rcovinv=invSafe(rcov);

        auto Hgp=univariate_gaussian_process_distribution::H
            (p.log_desc_beta_for_dL,rcov,rcovinv,
             std::exp(par.logvlogdL),std::exp(par.loglandalogdL));
        //        double logL_dL=univariate_gaussian_process_distribution::logL
        //            (par.mlogdL,std::exp(par.logvlogdL),std::exp(par.loglandalogdL),
        //             par.log_dL_j,p.log_desc_beta_for_dL);
        auto rSinv=invSafe(rcov*std::exp(par.logvlogdL));

        for (std::size_t i5=0; i5<par.log_dL_j.size(); ++i5)
          for (std::size_t j5=0; j5<par.log_dL_j.size(); ++j5)
            {
              out.push_back(-rSinv(i5,j5)*0.5);
            }
        out.push_back(Hgp(par.log_dL_j.size(),par.log_dL_j.size()));
        // out.push_back(Hgp(par.log_dL_j.size(),par.log_dL_j.size()+1));
        out.push_back(Hgp(par.log_dL_j.size()+1,par.log_dL_j.size()+1));

        // double logL_mlogdL=Normal(par.mlogdL,p.mmlogdL,p.smlogdL);
        out.push_back(-1.0/sqr(p.smlogdL));
        out.push_back(-0.5);

        //double logL_slogdL=Normal(par.logvlogdL,p.mlogslogdL,p.slogslogdL);
        out.push_back(-1.0/sqr(p.slogslogdL));
        out.push_back(-0.5);


        //double logL_landalogdL=Normal
        //     (par.loglandalogdL,p.mloglandalogdL,p.sloglandalogdL);
        out.push_back(-1.0/sqr(p.sloglandalogdL));
        out.push_back(-0.5);


        //double logL_mdelta=Normal(par.mlogdelta,p.mmlogdelta,p.smlogdelta);
        out.push_back(-1.0/sqr(p.smlogdelta));
        out.push_back(-0.5);


        //  double logL_sdelta=Normal(par.logvlogdelta, p.mlogslogdelta,p.slogslogdelta);
        out.push_back(-1.0/sqr(p.slogslogdelta));
        out.push_back(-0.5);


        // return logL_L0+out+logL_dL+logL_mdelta+logL_sdelta+logL_mlogdL+
        //     logL_slogdL+logL_landalogdL;

        return out;
      }
    };


    static std::string ClassName(){return "Master_Tempering_Gauss_Newton_Hessian";}

    const logL& logL_;

    const std::vector<Likelihood_Record>& v;
    const Prior& p;

    GN_H(const logL& L): logL_(L),v(L.s), p(L.p){}



    M_Matrix<double> operator()(const M_Matrix<double>& x) const
    {
      f j(v,p);
      d2 d(v,p);
      auto J=make_Jacobian_Finite_Difference(j)(x);
      auto D=d(x);
      return JTd2J(J,D);
    }

  };

*/

  struct Hinv
  {
    std::size_t diag_size;

    Hinv(std::size_t diagSize): diag_size(diagSize){}
    M_Matrix<double> operator()(const M_Matrix<double>& x)const
    {
      if (diag_size>0)
        {
          std::size_t n=x.nrows();
          std::size_t nA=n-diag_size;
          std::size_t nC=diag_size;
          M_Matrix<double> out(n,n);
          M_Matrix<double> A(nA,nA);
          M_Matrix<double> B(nA,nC);
          M_Matrix<double> Dinv(1,nC);
          for (std::size_t i=0; i<nA; ++i)
            for (std::size_t j=0; j<nA; ++j)
              A(i,j)=x(i,j);
          for (std::size_t i=0; i<nA; ++i)
            for (std::size_t j=0; j<nC; ++j)
              B(i,j)=x(i,nA+j);
          for (std::size_t j=0; j<nC; ++j)
            Dinv[j]=1.0/x(nA+j,nA+j);

          auto C=Transpose(B);
          auto Ainv=invSafe(A-xdiagXT(B,Dinv));
          if (Ainv.empty())
            return Ainv;
          else
            {
              auto Binv=MultDiag(-Ainv*B, Dinv);

              auto Di=diag(Dinv)-DiagMult(Dinv,C)*Binv;

              for (std::size_t i=0; i<nA; ++i)
                for (std::size_t j=0; j<nA; ++j)
                  out(i,j)=Ainv(i,j);
              for (std::size_t i=0; i<nA; ++i)
                for (std::size_t j=0; j<nC; ++j)
                  {
                    out(i,nA+j)=Binv(i,j);
                    out(nA+j,i)=Binv(i,j);
                  }

              for (std::size_t i=0; i<nC; ++i)
                for (std::size_t j=0; j<nC; ++j)
                  out(nA+i,nA+j)=Di(i,j);

              return out;
            }
        }
      else
        return invSafe(x);
    }


    M_Matrix<double> operator()(const M_Matrix<double>& H,double landa)const
    {
      M_Matrix<double> Hlanda=H+diag(H)*landa;
      return (*this)(Hlanda);
    }

  };

};



inline
std::ostream& operator<<(std::ostream& os, const Master_Tempering_Likelihood::Parameters_SE par)
{
  os<<"Master_Tempering_Likelihood::Parameters\n";
  os<<" L0="<<par.L0<<"\t"<<par.L0_se()<<"\n";
  os<<" mlogdL="<<par.mlogdL<<"\t"<<par.mlogdL_se()<<"\n";
  os<<" logvlogdL="<<par.logvlogdL<<"\t"<<par.logvdL_se()<<"\n";
  os<<" loglandalogdL="<<par.logrlandalogdL<<"\t"<<par.loglandadL_se()<<"\n";
  os<<" logepsilonlogdL="<<par.logepsilon2logL<<"\t"<<par.logepsilondL_se()<<"\n";
  os<<" mlogdelta="<<par.mlogdelta<<"\t"<<par.mlogdelta_se()<<"\n";
  os<<" logvlogdelta="<<par.logvlogdelta<<"\t"<<par.logvlogdelta_se()<<"\n";
  auto beta=par.lfit_.beta();
  os<<" deltaLikelihoods\n";
  for (std::size_t i=0; i<par.dL_j.size(); ++i)
    {
      os<<beta[i]<<"\t"<<par.dL_j[i]<<"\t"<<par.dL_j_se()[i]<<"\n";

    }
  auto Lse=par.L_j_se();
  auto L=par.Ls();
  os<<"logLiks\n";
  for (std::size_t i=0; i<par.dL_j.size(); ++i)
    {
      os<<beta[i]<<"\t"<<L[i]<<"\t"<<Lse[i]<<"\n";

    }
  os<<"Evidence="<<par.E()<<"\t"<<par.E_se()<<"\n";


  return os;

}





class Master_Adaptive_Beta_New
{
public:




  void push_step()
  {
    data_.back().newParticles();
  }

  template<class mcmc>
  void new_acceptance(std::size_t i, const mcmc& sDist,const mcmc& cDist)
  {
    this->data_.back().particle.back()[i].isValid=true;
    this->data_.back().particle.back()[i].logLikInit=sDist.logLik;
    this->data_.back().particle.back()[i].logLikNext=cDist.logLik;
  }
  template<class mcmc>
  void new_rjection(std::size_t i, const mcmc& sDist)
  {
    this->data_.back().particle.back()[i].isValid=false;
    this->data_.back().particle.back()[i].logLikInit=sDist.logLik;
  }


  void push_acceptance(double betaless, double betamore)
  {
    ++accepts_[{betaless,betamore}].first;

  }

  void push_rejection(double betaless, double betamore)
  {
    ++accepts_[{betaless,betamore}].second;
  }

  friend
  std::ostream& operator<<(std::ostream& os, const Master_Adaptive_Beta_New& me)
  {
    os<<"Beta\n"<<me.desc_beta_;
    os<<"History of likelihoods\n";

    return os;

  }



  void reset()
  {
    data_.clear();
  }

  Master_Adaptive_Beta_New(std::size_t Ninitial,  double beta_min, std::size_t N_2, double beta_infimo): desc_beta_{Ninitial,beta_min, N_2, beta_infimo}, data_{}
  {
    Likelihood_Record r;
    r.desc_beta=desc_beta_.getValue();
    data_.push_back(r);
  }




  std::size_t size()const {return desc_beta_.size();}

  Beta const& getBeta()const {return desc_beta_;}


private:

  Beta desc_beta_;
  std::map<std::pair<double,double>,std::pair<std::size_t,std::size_t>> accepts_;

  std::vector<Likelihood_Record>  data_;

};


struct LM_MultivariateGaussian: public MultivariateGaussian
{
  LM_MultivariateGaussian(MultivariateGaussian m
                          , Landa mylanda
                          ,double my_exp_delta_next_logL):
    MultivariateGaussian(m),
    landa(mylanda),
    exp_delta_next_logL(my_exp_delta_next_logL)
  {}
  LM_MultivariateGaussian():landa(),exp_delta_next_logL(){}
  Landa landa;
  double exp_delta_next_logL;
};

inline
std::ostream& operator<<(std::ostream& os,const LM_MultivariateGaussian& x)
{
  const MultivariateGaussian& b=x;
  os<<b;
  os<<"\nlanda\t"<<x.landa;
  os<<"\texp_delta_next_logL\t"<<x.exp_delta_next_logL<<"\n";
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

  double landa0_;
  double v_;
  std::size_t maxLoop_;

  //  double optimal_acceptance_rate_=0.574; // Handbook of Markov Chain Montecarlo p.100
  // optimal acceptance rate for Metropolis-Adjusted Langevin Algorithm

  double optimal_acceptance_rate_=0.234; // Handbook of Markov Chain Montecarlo p.96
  // optimal acceptance rate for Metropolis Algotithm


public:

  double optimal_acceptance_rate()const {return optimal_acceptance_rate_;}

  LevenbergMarquardt_step( double landa0,
                           double v,
                           std::size_t maxLoop)
    :landa0_{landa0},v_(v),maxLoop_(maxLoop)
  {
  }






  struct LM_logL:public D_logL
  {
    double logL;
    M_Matrix<double> Hinv;
    M_Matrix<double> d;
    double exp_delta_next_logL;
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
        out.exp_delta_next_logL=-0.5*multTransp(out.d,out.G)[0];
      }
    return out;
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, double beta)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta);
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, double beta,double landa)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta,landa);
  }


  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,double beta)const
  {
    mcmc_step<myPropDistribution> out(std::move(p),beta);
    return update_mcmc_step(L,model,data,out,beta);
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,double beta,double landa)const
  {
    mcmc_step<myPropDistribution> out(std::move(p),beta);
    return update_mcmc_step(L,model,data,out,beta,landa);
  }


  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post&& p,double beta)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),beta);
  }

  mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post&& p,double beta,double landa)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),beta,landa);
  }



  mcmc_step<myPropDistribution>& update_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_step<myPropDistribution>& out,double beta, AP t_)const
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
        while
            ((LM_logL.Hinv.empty()
              ||!t_(LM_logL.exp_delta_next_logL
                    ,cand.logbPL(beta)-out.logbPL()
                    ,out.param.size())
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


  mcmc_step<myPropDistribution>& update_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_step<myPropDistribution>& out,double beta, double landa)const
  {
    out.beta=beta;
    if (out.isValid)
      {
        LM_logL LM_logL=update_landa( out,landa,beta);
        if (!LM_logL.Hinv.empty())
          {
            M_Matrix<double> next=out.param+LM_logL.d;
            out.proposed=myPropDistribution(MultivariateGaussian
                                            (next,LM_logL.Hinv,LM_logL.H),landa);
          }
        else
          out.isValid=false;
      }
    return
        out;
  }




};



template<
    class D
    , template<class> class M
    , template<class,template<class> class > class D_Lik
    , class myPropDistribution>
class LevenbergMarquardt_step<D,M,D_Lik,myPropDistribution,Landa>
{
public:


  //  double optimal_acceptance_rate_=0.574; // Handbook of Markov Chain Montecarlo p.100
  // optimal acceptance rate for Metropolis-Adjusted Langevin Algorithm

  //  double optimal_acceptance_rate_=0.234; // Handbook of Markov Chain Montecarlo p.96
  // optimal acceptance rate for Metropolis Algotithm


public:


  struct LM_logL:public D_logL
  {
    double logL;
    M_Matrix<double> Hinv;
    M_Matrix<double> d;
    double exp_delta_next_logL;
  };

  static LM_logL update_landa(const mcmc_Dpost& postL,const Landa& landa,double beta)
  {
    LM_logL out;
    out.G=postL.D_prior.G+postL.D_lik.G*beta;
    out.H=postL.D_prior.H+postL.D_lik.H*beta;
    out.logL=postL.logbPL(beta);
    for (std::size_t i=0; i<out.H.nrows(); ++i)
      //    out.H(i,i)=out.H(i,i)+postL.D_lik.H(i,i)*beta*landa;
      // this alternative does not work
      out.H(i,i)=out.H(i,i)*(1+landa.getValue());
    out.Hinv=invSafe(out.H);
    if (!out.Hinv.empty())
      {
        out.d=-(out.G*out.Hinv);
        out.exp_delta_next_logL=-0.5*multTransp(out.d,out.G)[0];
      }
    return out;
  }

  static mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param,const Landa& landa, double beta, std::size_t iscout)
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),landa,beta, iscout);
  }


  static mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,const Landa& landa,double beta, std::size_t iscout)
  {
    mcmc_step<myPropDistribution> out(std::move(p),beta,iscout);
    return update_mcmc_step(L,model,data,out,landa,beta);
  }


  static mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post&& p,const  Landa& landa,double beta, std::size_t iscout)
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),landa ,beta, iscout);
  }





  static  mcmc_step<myPropDistribution>& update_mcmc_step(const D_Lik<D,M> , const M<D>& , const D& , mcmc_step<myPropDistribution>& out,const Landa& landa,double beta)
  {
    out.beta=beta;
    if (out.isValid)
      {
        LM_logL LM_logL=update_landa( out,landa,beta);
        if (!LM_logL.Hinv.empty())
          {
            M_Matrix<double> next=out.param+LM_logL.d;
            out.proposed=myPropDistribution(MultivariateGaussian
                                            (next,LM_logL.Hinv,LM_logL.H),
                                            landa,
                                            LM_logL.exp_delta_next_logL);
            if ((out.proposed.logDetCov()>0)|| std::isnan( out.proposed.logDetCov()))
              {
                out.isValid=false;
                std::cerr<<"\ninvalid logDetCov\n";
              }
          }
        else
          out.isValid=false;
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
  os<<"\nlanda0\n"<<x.landa0_;
  os<<"\nlanda_mult_factor\n"<<x.v_;
  os<<"\nmaxLoop\n"<<x.maxLoop_;
  return os;
}

template<class mcmc>
class Tempered_Evidence_Evaluation
{
  SamplesSeries<std::pair<Beta,std::vector<mcmc>>> run_;
  std::pair<double,double> logEvidence_;
public:

  std::pair<double,double> logEvidence()const {return logEvidence_;}
  const SamplesSeries<std::pair<Beta,std::vector<mcmc>>>& samples()const
  {
    return run_;
  }

  Tempered_Evidence_Evaluation(SamplesSeries<std::pair<Beta,std::vector<mcmc>>> o)
    :run_(std::move(o))
  {
    logEvidence_=run_.mean_var(&Evidence_pair,run_.size()/2);


  }

  static double Evidence_pair(const std::pair<Beta,std::vector<mcmc>>& sample)
  {
    return Evidence(sample.first,sample.second);
  }
  static double Evidence(const Beta mybeta,const std::vector<mcmc>& sample)
  {
    auto n=mybeta.getValue().size();
    double beta0=mybeta.getValue()[0];
    double sum=0;
    double sumdb=0;
    double logLik0=sample[0].logLik;
    for (std::size_t i=1; i<sample.size(); ++i)
      {
        double beta=mybeta.getValue()[i];
        double db=beta0-beta;
        double logLik=sample[i].logLik;
        sum+=db*(logLik0+logLik)/2;
        sumdb+=db;
        if ((i==n-1)&&(beta>0))
          {
            double db0=beta;
            sum+=(db0*(1+db0/db/2))*logLik-sqr(db0)/db/2*logLik0;
            sumdb+=db0;
          }
        beta0=beta;
        logLik0=logLik;
      }
    return sum;

  }


};



template
<
    class D
    , template<class> class M
    , template<class,template<class> class > class D_Lik//=Poisson_DLikelihood
    , class my_PropD//=LM_MultivariateGaussian
    , class AP//=Landa
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
                       ,std::mt19937_64 &mt,
                       double& dHd,
                       double& logQforward,
                       double& logQbackward )
    {
      if (!cLik.isValid)
        {
          return false;
        }
      else
        {
          double logPcurrent=sLik.logbPL();
          double logPcandidate=cLik.logbPL();

          logQforward=sLik.proposed.logP(cLik.param);
          logQbackward=cLik.proposed.logP(sLik.param);

          double logA=logPcandidate-logQforward-(logPcurrent-logQbackward);
          double A=std::min(1.0,exp(logA));
          std::uniform_real_distribution<double> u(0,1);
          auto d=cLik.param-sLik.param;
          auto H=sLik.beta*sLik.D_lik.H+sLik.D_prior.H;
          dHd=0.5*xTSigmaX(d,H);

          double r=u(mt);
          bool accept_=r<A;
          if (accept_)
            {
              return true;
            }
          else
            {
              return false;
            }
        }
    }



    static bool accept_swap(const mcmc_step<pDist>& sLik
                            ,const mcmc_step<pDist>& cLik
                            ,std::mt19937_64 &mt, double& s_dHd, double& c_dHd)
    {
      if (!cLik.isValid)
        {
          return false;
        }
      else
        {
          double logPcurrent=sLik.logbPL();
          double logPcandidate=cLik.logbPL();

          double logPSwapcurrent=sLik.logbPLb(cLik.beta);
          double logPSwapcandidate=cLik.logbPLb(sLik.beta);


          double logA=logPSwapcandidate+logPSwapcurrent-(logPcurrent+logPcandidate);
          double A=std::min(1.0,exp(logA));
          std::uniform_real_distribution<double> u(0,1);
          auto d=cLik.param-sLik.param;
          auto s_H=sLik.beta*sLik.D_lik.H+sLik.D_prior.H+sLik.beta*cLik.D_lik.H+cLik.D_prior.H;
          auto c_H=cLik.beta*sLik.D_lik.H+sLik.D_prior.H+cLik.beta*cLik.D_lik.H+cLik.D_prior.H;
          s_dHd=0.25*xTSigmaX(d,s_H);
          c_dHd=0.25*xTSigmaX(d,c_H);

          double r=u(mt);
          bool accept_=r<A;
          if (accept_)
            {
              return true;
            }
          else
            {
              return false;
            }
        }
    }



  };

public:



  static std::size_t min_tryParameter(){return 20;}




  static SamplesSeries<mcmc_step<pDist>>
  run_new
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   Adaptive_discrete<AP>& landa,
   const std::tuple<double,std::size_t,std::size_t>& betas,
   std::mt19937_64& mt
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::size_t nsamples=std::get<1>(betas);
    std::size_t nskip=std::get<2>(betas);

    if (sDist.isValid)
      {
        SamplesSeries<mcmc_step<pDist>> o(nsamples);
        landa.actualize();
        mcmc_step<pDist> cDist;
        n_steps_try(LM_Lik,lik,model,data,sDist,cDist,landa,mt,nskip,os,startTime,timeOpt);
        landa.actualize();
        //std::cerr<<landa;
        //os<<landa;
        landa=landa.partialReset();

        while (!o.full())
          {
            n_steps_try(LM_Lik,lik,model,data,sDist,cDist,landa,mt,nskip,os,startTime,timeOpt);
            o.push_back(sDist);
            landa.actualize();
            std::cerr<<landa;
            os<<landa;
          }
        std::pair<double,double> l=o.mean_std
            ([](const mcmc_step<pDist>& mc){return mc.logLik;},nsamples/2);

        std::cout<<"\nmcmc:: beta=\t"<<sDist.beta<<"\tlogLik\t"<<l<<"\tnsamples\t"<<nsamples/2<<"\n";
        os<<"\nmcmc:: beta=\t"<<sDist.beta<<"\tlogLik\t"<<l<<"\n";

        return o;
      }
    else
      return {};
  }



  static SamplesSeries<mcmc_step<pDist>> run
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   std::size_t nsamples,
   std::size_t nskip,
   std::mt19937_64& mt
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::size_t nsamples_0=min_tryParameter();
    std::size_t nsamplesFinal=std::max(40lu,nsamples*nskip/20);
    if (!sDist.isValid)
      return {};
    M_Matrix<double> c=sDist.proposed.sample(mt);
    mcmc_step<pDist> cDist=LM_Lik.get_mcmc_step(lik,model,data,c,sDist.beta);

    AP r_optimal=adapt_Parameter
        (LM_Lik,lik,model,data,sDist,cDist,nsamples_0,nsamplesFinal,mt,os,startTime,timeOpt);
    auto LM=LM_Lik(r_optimal);

    std::size_t naccepts=0;
    std::size_t nrejects=0;

    SamplesSeries<mcmc_step<pDist>> o(nsamples);
    while(true)
      {
        os<<"niter::\t"<<o.size()<<"\t";
        std::cout<<"niter::\t"<<o.size()<<"\t";

        n_steps(LM,lik,model,data,sDist,cDist,nskip,naccepts,nrejects,mt,os,startTime,timeOpt);
        o.push_back(sDist);
        if (o.full())
          break;
      }
    std::pair<double,double> l=o.mean_std
        ([](const mcmc_step<pDist>& mc){return mc.logLik;},nsamples/2);

    std::cout<<"\nmcmc:: beta=\t"<<sDist.beta<<"\tlogLik\t"<<l<<"\tnsamples\t"<<nsamples/2<<"\n";
    os<<"\nmcmc:: beta=\t"<<sDist.beta<<"\tlogLik\t"<<l<<"\n";
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
   std::mt19937_64& mt
   ,std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::size_t naccepts=0;
    std::size_t nrejects=0;
    auto LM=LM_Lik(r_value);
    sDist=LM.get_mcmc_step(lik,model,data,sDist.param,sDist.beta);
    M_Matrix<double> c=sDist.proposed.sample(mt);
    cDist=LM.get_mcmc_step(lik,model,data,c,sDist.beta);

    n_steps(LM,lik,model,data,sDist,cDist,nsamples,naccepts,nrejects,mt,os,startTime,timeOpt);
    return {naccepts,nrejects};

  }


  static void
  try_Parameters
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   mcmc_step<pDist>& cDist
   ,Adaptive_discrete<AP>& landa
   ,std::size_t nsamples_per_par,
   std::mt19937_64& mt
   ,std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt
   )
  {
    landa.actualize();
    n_steps_try(LM_Lik,lik,model,data,sDist,cDist,landa,mt,nsamples_per_par,os,startTime,timeOpt);
    landa.actualize();
    // std::cerr<<landa;
    // os<<landa;
    landa=landa.partialReset();

    for (std::size_t i=0; i< 5; ++i)
      {
        n_steps_try(LM_Lik,lik,model,data,sDist,cDist,landa,mt,nsamples_per_par,os,startTime,timeOpt);
        landa.actualize();
        //  std::cerr<<landa;
        //  os<<landa;
      }
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
   std::mt19937_64& mt
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::size_t i=0;
    AP r_value=LM.getValue();

    while(i<nsamples)
      {

        auto tnow=std::chrono::steady_clock::now();
        auto d=tnow-startTime;
        double t0=timeOpt;
        timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
        auto timeIter_=60*(timeOpt-t0);


        std::cout<<i<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<sDist.beta<<"\t"<<r_value<<"\t";
        test::put(std::cout,sDist,cDist);
        put(std::cout,sDist,cDist);

        os<<"n_steps::"<<i<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<sDist.beta<<" "<<r_value<<" ";
        test::put(os,sDist,cDist);
        put(os,sDist,cDist);



        if (step(LM,lik,model,data,sDist,cDist,mt))
          {
            ++naccepts;
          }
        else
          ++nrejects;
        ++i;
        std::cout<<" "<<naccepts<<" "<<nrejects<<"\n";
        os<<" "<<naccepts<<" "<<nrejects<<"\n";
        os.flush();

      }

  }


  static void n_steps_try
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   mcmc_step<pDist>& cDist,
   Adaptive_discrete<AP>& pars,
   std::mt19937_64& mt
   , std::size_t nsteps
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {

    for( std::size_t i=0;i<nsteps; ++i)
      {

        AP landa=pars.sample(mt);
        double dHd;
        double logQforward;
        double logQbackward;
        if (step(LM_Lik,lik,model,data,sDist,cDist,landa,
                 dHd,logQforward, logQbackward,mt,i,os,startTime,timeOpt))
          {
            pars.push_acceptance(landa,dHd);
          }
        else
          pars.push_rejection(landa);
      }


  }


  static void n_steps_template_try
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,std::vector<mcmc_step<pDist>>& sDist,
   std::vector<Adaptive_discrete<AP>>& pars,
   std::vector<std::mt19937_64>& mt
   ,double pJump,
   std::size_t nsteps
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {

    for( std::size_t i=0;i<nsteps; ++i)
      {
        std::vector<double> dHd;
        tempered_step(LM_Lik,lik,model,data,sDist,pars,dHd,pJump,mt,i,os,startTime,timeOpt);
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
   std::mt19937_64& mt
   ,std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    AP ap0=AP::max();
    double opt_accpt=LM_Lik.optimal_acceptance_rate();
    double x0=std::log(ap0.getValue());
    double xmax=x0;
    double dx=AP::logFactor();
    gaussian_lin_regr lr(std::log(AP::min().getValue()),std::log(AP::max().getValue()));

    double y=logit(opt_accpt);

    std::size_t nsamples=nsamples_0;
    auto res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);
    auto m0=logit(res0);
    lr.push_back(x0,m0.first,m0.second );
    x0-=dx;
    ap0.setValue(std::exp(x0));
    res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);
    m0=logit(res0);
    lr.push_back(x0,m0.first,m0.second );

    while (res0.first<res0.second)
      {
        x0-=dx;
        ap0.setValue(std::exp(x0));
        res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);
        m0=logit(res0);
        lr.push_back(x0,m0.first,m0.second );
      }
    x0=lr.get_optimal_x(y,os);
    ap0.setValue(std::exp(x0));
    res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);

    double x00=xmax;
    while (!accept_Parameter(opt_accpt,res0)&&x0!=x00)
      {
        m0=logit(res0);
        lr.push_back(x0,m0.first,m0.second );
        x00=x0;
        x0=lr.get_optimal_x(y,os);

        ap0.setValue(std::exp(x0));
        res0=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);
      }
    while (nsamples<nsamplesFinal)
      {
        res0+=try_Parameter(LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);
        nsamples*=2;

        x00=xmax+1;
        while (!accept_Parameter(opt_accpt,res0)&&x0!=x00)
          {
            m0=logit(res0);
            lr.push_back(x0,m0.first,m0.second );
            x00=x0;
            x0=lr.get_optimal_x(y,os);
            if (x0!=x00)
              {
                ap0.setValue(std::exp(x0));
                res0=try_Parameter
                    (LM_Lik,lik,model,data,ap0,sDist,cDist,nsamples,mt,os,startTime,timeOpt);
              }
          }
      }
    return ap0;

  }


  static void
  adapt_Parameter_new
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist,
   Adaptive_discrete<AP>& aps,
   std::size_t nsamples_0,
   std::mt19937_64& mt
   ,std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {

    std::size_t nsamples=nsamples_0;
    mcmc_step<pDist> cDist;
    try_Parameters(LM_Lik,lik,model,data,sDist,cDist,aps,nsamples,mt,os,startTime,timeOpt);





  }





  static std::ostream& put
  (std::ostream& os
   ,mcmc_step<pDist>& sLik
   ,mcmc_step<pDist>& cLik,
   const double& dHd,
   const double& logQforward,const double& logQbackward)
  {
    os<<" dHd "<<dHd<<" "<<" logQf "<<logQforward<<" logQb "<<logQbackward;
    os<<sLik.logbPL()<<" "<<cLik.logbPL()<<" ";
    os<<sLik.logPrior<<" "<<cLik.logPrior<<" ";
    os<<sLik.logLik<<" "<<cLik.logLik<<" ";
    os<<sLik.proposed.landa<<" "<<cLik.proposed.landa<<" ";
    //os<<sLik.proposed.exp_delta_next_logL/sLik.param.size()<<" ";
    //os<<cLik.proposed.exp_delta_next_logL/cLik.param.size()<<" ";
    return os;
  }

  static std::ostream& put
  (std::ostream& os
   ,const std::vector<mcmc_step<pDist>>& sLik
   ,const std::vector<bool>& accepts,
   const std::vector<double>& dHd,
   const std::vector<double>& logQforward,
   const std::vector<double>& logQbackward
   )
  {
    for (std::size_t  i=0; i<sLik.size(); ++i)
      {
        os<<sLik[i].iscout<<"\t";
        os<<sLik[i].beta<<"\t";
        if(accepts[i]) os<<"acc\t";
        else
          os<<"rej\t";
        os<<sLik[i].proposed.landa<<"\t";
        os<<" dHd "<<dHd[i]<<"\t ";
        os<<" logQf "<<logQforward[i]<<"\t ";
        os<<" logQb "<<logQbackward[i]<<"\t ";
        os<<sLik[i].logbPL()<<"\t";
        os<<sLik[i].logPrior<<"\t";
        os<<sLik[i].logLik<<"\n";
      }
    return os;
  }



  static bool step
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,mcmc_step<pDist>& sDist
   ,mcmc_step<pDist>& cDist
   ,const AP& landa
   ,double& dHd
   ,double& logQforward
   ,double& logQbackward
   ,std::mt19937_64& mt
   ,std::size_t i
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    bool out;
    LM_Lik.update_mcmc_step(lik,model,data,sDist,landa,sDist.beta);
    M_Matrix<double> c=sDist.proposed.sample(mt);
    cDist=LM_Lik.get_mcmc_step(lik,model,data,c,landa,sDist.beta,sDist.iscout);
    if (test::accept(sDist,cDist,mt,dHd,logQforward,logQbackward))
      {
        out =true;
      }
    else
      {
        out=false;
      }
    auto tnow=std::chrono::steady_clock::now();
    auto d=tnow-startTime;
    double t0=timeOpt;
    timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
    auto timeIter_=60*(timeOpt-t0);
    std::cout<<i<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<sDist.beta;
    //test::put(std::cout,sDist,cDist);
    put(std::cout,sDist,cDist,dHd,logQforward,logQbackward);

    os<<"n_steps::"<<i<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<sDist.beta;
    test::put(os,sDist,cDist);
    put(os,sDist,cDist,dHd,logQforward,logQbackward);
    if (out)
      {
        sDist=std::move(cDist);
        std::cout<<"Acc\n";
        os<<"Acc\n";
      }
    else
      {
        std::cout<<"Rej\n";
        os<<"Rej\n";

      }



    return out;
  }


  static void tempered_step
  (const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   ,std::vector<mcmc_step<pDist>>& sDist
   ,std::vector<Adaptive_discrete<AP>>& pars
   ,Master_Adaptive_Beta_New& aBeta
   ,std::vector<double>& dHd
   ,std::vector<double>& logQforward
   ,std::vector<double>& logQbackward
   , double p_Tjump
   ,std::vector<std::mt19937_64>& mt
   ,std::size_t isamples
   ,std::size_t isubSamples
   , bool does_stdout
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::vector<bool> out(sDist.size());
    aBeta.push_step();

    //#pragma omp parallel for
    for(std::size_t i=0; i<sDist.size(); ++i)
      {
        AP landa;


        mcmc_step<pDist> cDist;
        //      std::size_t icount=0;
        while(!cDist.isValid)
          {
            landa=pars[i].sample(mt[i]);
            LM_Lik.update_mcmc_step
                (lik,model,data,sDist[i],landa,aBeta.getBeta().getValue()[i]);
            M_Matrix<double> c=sDist[i].proposed.sample(mt[i]);
            cDist=LM_Lik.get_mcmc_step
                (lik,model,data,c,landa,aBeta.getBeta().getValue()[i],sDist[i].iscout);
            if (!cDist.isValid)
              pars[i].push_rejection(landa);

          }

        if (test::accept(sDist[i],cDist,mt[i],dHd[i],logQforward[i],logQbackward[i]))
          {
            aBeta.new_acceptance(i,sDist[i],cDist);
            pars[i].push_acceptance(landa,dHd[i]);
            sDist[i]=std::move(cDist);
            out[i] =true;
          }
        else
          {
            aBeta.new_rjection(i,sDist[i]);
            out[i]=false;
            pars[i].push_rejection(landa);

          }

      }


    auto tnow=std::chrono::steady_clock::now();
    auto d=tnow-startTime;
    double t0=timeOpt;
    timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
    auto timeIter_=60*(timeOpt-t0);
    double evidence=Tempered_Evidence_Evaluation<mcmc_step<pDist>>::Evidence(aBeta.getBeta(),sDist);
    std::cout<<"isample::"<<isamples<<"\t"<<"isubSample::"<<isubSamples<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<"Evidence\t"<<evidence<<"\n";

    os<<"isample::"<<isamples<<"\t"<<"isubSample::"<<isubSamples<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<"Evidence\t"<<evidence<<"\n";

    put(os,sDist,out,dHd,logQforward,logQbackward);
    if (does_stdout)
    put(std::cout,sDist,out,dHd,logQforward,logQbackward);


    for(std::size_t i=1; i<sDist.size(); ++i)
      {
        std::uniform_real_distribution<> u;
        double r=u(mt[i]);
        double s_dHd, c_dHd;
        if (r<p_Tjump)
          {
            if (test::accept_swap(sDist[i-1],
                                  sDist[i],
                                  mt[i], s_dHd, c_dHd))
              {
                std::swap(sDist[i-1],sDist[i]);
                std::swap(sDist[i-1].beta,sDist[i].beta);
                os<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                if (does_stdout)
                std::cout<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                aBeta.push_acceptance(sDist[i-1].beta,sDist[i].beta);
              }
            else
              aBeta.push_rejection(sDist[i-1].beta,sDist[i].beta);
            ;
          }
      }

  }






};



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
std::ostream& operator<<(std::ostream& os
                         ,const Metropolis_Hastings_mcmc<D,M,D_Lik,my_PropD,AP,propDistStep>& mcmc)
{
  return os;
}





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
    std::pair<double,double>  l=o[0].second.mean_var
        ([](const mcmc& mc){return mc.logLik;},o[0].second.size()/2);

    double beta=o[0].first;
    double sumdb=0;

    for (std::size_t i=1; i<o.size(); ++i)
      {
        double beta0=beta;
        auto l0=l;
        beta=o[i].first;
        SamplesSeries<mcmc>& s=o[i].second;
        std::size_t nsamples=s.size();
        l=s.mean_var
            ([](const mcmc& mc){return mc.logLik;},nsamples/2);
        double db=beta-beta0;
        sum+=db*(l.first+l0.first)/2;
        sumVar+=db*(l.second+l0.second)/2;
        sumdb+=db;
        if ((i==1)&&(beta0>0))
          {
            double db0=beta0;
            sum+=(db0*(1+db0/db/2))*l0.first-sqr(db0)/db/2*l.first;
            sumVar+=(db0*(1+db0/db/2))*l0.second+sqr(db0)/db/2*l.second;
            sumdb+=db0;
          }
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

template<class mcmc>
std::ostream& operator<<(std::ostream& os,Tempered_Evidence_Evaluation<mcmc>& x)
{
  os<<"\nlogEvidence\n"<<x.logEvidence();
  os<<"\nsamples\n"<<x.samples();
  return os;
}


template<
    class D
    , template<class>   class M
    , template<class,template<class> class >   class D_Lik//=Poisson_DLikelihood
    , class my_PropD//=LM_MultivariateGaussian
    , class AP//=trust_region
    ,template<
      class
      , template<class> class
      , template<class,template<class> class > class
      , class
      , class
      >
    class propDist//=LevenbergMarquardt_step
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
   ,const propDist<D,M,D_Lik,my_PropD,AP>& LMLik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   , const Adaptive_discrete<AP>& landa_Dist0
   ,const std::vector<std::tuple<double,std::size_t,std::size_t>>& beta
   ,std::mt19937_64& mt
   ,std::ofstream& os
   ,const std::chrono::steady_clock::time_point& startTime
   ,double& timeOpt)
  {

    std::vector<std::pair<double,SamplesSeries<mystep>>> o;
    Adaptive_discrete<AP> landa_Dist(landa_Dist0);

    M_Matrix<double> pinit;
    double beta0=std::get<0>(beta[0]);
    std::size_t ntrials=0;
    mcmc_post postL;
    while(!postL.isValid)
      {
        pinit=model.sample(mt);
        postL=lik.get_mcmc_Post(model,data,pinit);
        ++ntrials;
      }

    AP landa=landa_Dist.sample(mt);
    mystep cDist=LMLik.get_mcmc_step(lik,model,data,pinit,std::move(postL),landa,beta0,0);

    for (std::size_t i=0; i<beta.size(); ++i)
      {
        double betaval=std::get<0>(beta[i]);
        LMLik.update_mcmc_step(lik,model,data,cDist,landa,betaval);
        Adaptive_discrete<AP> landa_Dist_run=landa_Dist.partialReset();

        auto s=mcmc.run_new
            (LMLik,lik,model,data,cDist,landa_Dist_run,beta[i],mt,os,startTime,timeOpt);
        o.push_back({betaval,s});
      }
    auto out=new Evidence_Evaluation<mystep>(std::move(o));
    std::cout<<"LogEvidence= "<<out->logEvidence()<<"\n";


    return out;
  }
};


template<class Landa>
class Tempering_Dynamics_Model
{
public:
  struct history
  {
    std::map<std::size_t, std::vector<std::pair<Landa,double>>> h;
  };




};






template<
    class D
    , template<class>   class M
    , template<class,template<class> class >   class D_Lik//=Poisson_DLikelihood
    , class my_PropD//=LM_MultivariateGaussian
    , class AP//=trust_region
    ,template<
      class
      , template<class> class
      , template<class,template<class> class > class
      , class
      , class
      >
    class propDist//=LevenbergMarquardt_step
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
class Template_Tempering_mcmc
{
public:
  typedef mcmc_step<my_PropD> mystep;

  static void
  run
  (const MH<D,M,D_Lik,my_PropD,AP,propDist>& mcmc
   ,const propDist<D,M,D_Lik,my_PropD,AP>& LMLik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   , const Adaptive_discrete<AP>& landa_Dist0
   ,const Master_Adaptive_Beta_New& beta0
   , double maxTime
   , std::size_t nsamples
   ,std::size_t nskip
   ,double pTjump
   ,std::mt19937_64::result_type seed
   ,std::string EviName
   ,std::string EviNameLog0
   ,const std::chrono::steady_clock::time_point& startTime
   ,double& timeOpt
   ,std::size_t maxSimFileSize
   ,bool does_stdout)
  {
    std::mt19937_64 mt;
    mt.seed(seed);
    std::size_t i_sim=0;

    std::string f_par_name0=EviName+"_par.txt";
    std::string f_logL_name0=EviName+"_logL.txt";
    std::string f_fit_name0=EviName+"_fit.txt";
    std::string f_sim_name0=EviName+"_sim.txt";

    std::string f_par_name=f_par_name0    +"."+std::to_string(i_sim);
    std::string f_logL_name=f_logL_name0    +"."+std::to_string(i_sim);
    std::string f_fit_name=f_fit_name0    +"."+std::to_string(i_sim);
    std::string f_sim_name=f_sim_name0    +"."+std::to_string(i_sim);
    std::string f_log_name=EviNameLog0    +"."+std::to_string(i_sim);

    std::ofstream os;
       std::ofstream f_par;
    std::ofstream f_logL;
    std::ofstream f_sim;
    std::ofstream f_fit;

    os.open(f_log_name.c_str(), std::ofstream::out | std::ofstream::app);
    f_par.open(f_par_name.c_str(), std::ofstream::out | std::ofstream::app);
    f_logL.open(f_logL_name.c_str(), std::ofstream::out | std::ofstream::app);
    f_sim.open(f_sim_name.c_str(), std::ofstream::out | std::ofstream::app);
    f_fit.open(f_fit_name.c_str(), std::ofstream::out | std::ofstream::app);



    std::stringstream ss;

    ss<<"model\t";
    ss<<"seed\t";
    ss<<"time\t";
    ss<<"nsteps\t";
    ss<<"nsample\t";
    ss<<"Evidence\t";

    auto s=mystep{};
    s.writelogLHeaderDataFrame(ss);


    f_logL<<ss.str()<<std::endl;
    f_par<<ss.str()<<"\t";
    f_sim<<ss.str()<<"\t";
    f_fit<<ss.str()<<"\t";



    auto param=model.getPrior();
    param.writeHeaderDataFrame(f_par);
    f_par<<std::endl;

    auto sim=model.getSimulation(param);
    sim.writeHeaderDataFrame(f_sim);
    f_sim<<std::endl;

    auto& cl=model.getLikelihood();
    cl.writeYfitHeaderDataFrame(f_fit);
    f_fit<<std::endl;

    Master_Adaptive_Beta_New beta(beta0);
    auto n=beta.size();
    std::vector<mystep> sDists(n);

    std::vector<Adaptive_discrete<AP>> pars(n,landa_Dist0);
    std::vector<double> dHd(beta.size());
    std::vector<double> logQforward(beta.size());
    std::vector<double> logQbackward(beta.size());
    std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
    std::vector<std::mt19937_64> mts(n);
    for (std::size_t i=0; i<n; ++i)
      mts[i].seed(useed(mt));
    std::size_t o=0;
    for (std::size_t i=0; i<beta.size(); ++i)
      {
        M_Matrix<double> pinit;
        double beta0=beta.getBeta().getValue()[i];
        std::size_t ntrials=0;
        mcmc_post postL;
        while(!postL.isValid)
          {
            pinit=model.sample(mt);
            postL=lik.get_mcmc_Post(model,data,pinit);
            ++ntrials;
          }
        AP landa=pars[i].sample(mts[i]);
        sDists[i]=LMLik.get_mcmc_step(lik,model,data,pinit,std::move(postL),landa,beta0,i);
        LMLik.update_mcmc_step(lik,model,data,sDists[i],landa,beta0);
      }


    while (o<nsamples&&timeOpt<maxTime*60)
      {
        //        for (std::size_t i=0; i<n;++i)
        //          pars[i].actualize();



        for( std::size_t i=0;i<nskip; ++i)
          {
            mcmc.tempered_step
                (LMLik,lik,model,data,sDists,pars,beta,dHd,logQforward,logQbackward,pTjump,mts,o,i,does_stdout,os,startTime,timeOpt);
          }

        o++;

        auto tnow=std::chrono::steady_clock::now();
        auto d=tnow-startTime;
        double t0=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
        double evidence=Tempered_Evidence_Evaluation<mystep>::Evidence(beta.getBeta(),sDists);

        std::stringstream ss1;
        ss1<<model.id()<<"\t";
        ss1<<seed<<"\t";
        ss1<<t0<<"\t";
        ss1<<o*sDists.size()<<"\t";
        ss1<<o<<"\t";
        ss1<<evidence<<"\t";

        for (std::size_t i=0; i<sDists.size(); ++i)
          {

            auto param=model.getParameter(sDists[i].param);
            auto sim=model.getSimulation(sDists[i].param);
            (sDists[i].writelogLRowDataFrame(f_logL<<ss1.str()))<<"\n";

            (sDists[i].writelogLRowDataFrame(f_par<<ss1.str()))<<"\t";
            param.writeRowDataFrame(f_par)<<"\n";

            (sDists[i].writelogLRowDataFrame(f_sim<<ss1.str()))<<"\t";
            sim.writeRowDataFrame(f_sim)<<"\n";
            (sDists[i].writelogLRowDataFrame(f_fit<<ss1.str()))<<"\t";
            sDists[i].writeYfitRowDataFrame(f_fit)<<"\n";
            f_logL.flush();
            f_par.flush();
            f_fit.flush();

            f_sim.flush();
            os.flush();
            std::size_t f_sim_pos=f_sim.tellp()+f_logL.tellp()+f_fit.tellp()+f_par.tellp()+os.tellp();
            if (f_sim_pos>maxSimFileSize)
              {
                 f_sim.close();
                 f_fit.close();
                 f_logL.close();
                 f_par.close();
                 os.close();

                 std::cerr<<f_sim_name<<"is completed !!\n";
                 std::cerr<<f_fit_name<<"is completed !!\n";
                 std::cerr<<f_logL_name<<"is completed !!\n";
                 std::cerr<<f_log_name<<"is completed !!\n";
                 std::cerr<<f_par_name<<"is completed !!\n";

                 ++i_sim;

                 f_sim_name=f_sim_name0    +"."+std::to_string(i_sim);
                 f_sim.open(f_sim_name.c_str(), std::ofstream::out | std::ofstream::app);

                 f_log_name=EviNameLog0+"."+std::to_string(i_sim);
                 os.open(f_log_name.c_str(), std::ofstream::out | std::ofstream::app);


                 f_par_name=f_par_name0    +"."+std::to_string(i_sim);
                 f_par.open(f_par_name.c_str(), std::ofstream::out | std::ofstream::app);

                 f_fit_name=f_fit_name0    +"."+std::to_string(i_sim);
                 f_fit.open(f_fit_name.c_str(), std::ofstream::out | std::ofstream::app);

                 f_logL_name=f_logL_name0    +"."+std::to_string(i_sim);
                 f_logL.open(f_logL_name.c_str(), std::ofstream::out | std::ofstream::app);

                 std::cerr<<f_sim_name<<"is opened !!\n";
                 std::cerr<<f_par_name<<"is opened !!\n";
                 std::cerr<<f_logL_name<<"is opened !!\n";
                 std::cerr<<f_fit_name<<"is opened !!\n";
                 std::cerr<<f_log_name<<"is opened !!\n";
              }


          }




        for (std::size_t i=0; i<n;++i)
          {
            pars[i].actualize();
          }
        // std::cerr<<pars;
        os<<pars;
        // std::cout<<beta;
        os<<beta;
        if (n!=beta.size())
          {
            for (std::size_t i=n; i<beta.size(); ++i)
              {
                M_Matrix<double> pinit;
                double beta0=beta.getBeta().getValue()[i];
                std::size_t ntrials=0;
                mcmc_post postL;
                while(!postL.isValid)
                  {
                    pinit=model.sample(mt);
                    postL=lik.get_mcmc_Post(model,data,pinit);
                    ++ntrials;
                  }
                pars.push_back(pars[n-1]);
                mts.push_back(std::mt19937_64{});
                mts[i].seed(useed(mt));

                AP landa=pars[i].sample(mts[i]);
                sDists.push_back
                    (LMLik.get_mcmc_step(lik,model,data,pinit,std::move(postL),landa,beta0,i));
                LMLik.update_mcmc_step(lik,model,data,sDists[i],landa,beta0);
              }

          }

        n=beta.size();
      }


  }




};







#endif // EVIDENCE_H
