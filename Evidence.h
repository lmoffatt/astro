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

class gaussian_lin_regr
{
public:
  
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
  mcmc_step(mcmc_Dpost &&p, double beta_): mcmc_Dpost(p),beta(beta_), proposed(){}
  double beta;
  double logbPL()const {return mcmc_post::logbPL(beta);}
  double logbPLb(double mybeta)const {return mcmc_post::logbPL(mybeta);}
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
    this->p_=OptimalDistribution::optimal(&expectedVelocity,parDist_.first,this->p_,parDist_.second);
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


  Beta(std::size_t n, double min_beta):beta_(n)
  {
    beta_[0]=1;
    double f=std::pow(min_beta,1.0/(n-1));
    for (std::size_t i=1; i<n; ++i)
      beta_[i]=beta_[i-1]*f;

  }

  Beta():beta_(){}

  std::vector<double>const & getValue()const {return beta_;}

  std::vector<double>& getValue() {return beta_;}

private:
  std::vector<double> beta_;

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
  double nu_;  // reciprocal of the initial amplitude of adjustments
  std::size_t t0_;  //lag parameter
  double beta_min_;  // minimal value of beta
  std::size_t N_; // number of beta intervals

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
  
  static mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param,const Landa& landa, double beta)
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),landa,beta);
  }
  
  
  static mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_Dpost&& p,const Landa& landa,double beta)
  {
    mcmc_step<myPropDistribution> out(std::move(p),beta);
    return update_mcmc_step(L,model,data,out,landa,beta);
  }
  
  
  static mcmc_step<myPropDistribution> get_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, const M_Matrix<double>& param, mcmc_post&& p,const  Landa& landa,double beta)
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),landa ,beta);
  }
  
  
  
  
  
  static  mcmc_step<myPropDistribution>& update_mcmc_step(const D_Lik<D,M> L, const M<D>& model, const D& data, mcmc_step<myPropDistribution>& out,const Landa& landa,double beta)
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
  std::pair<double,double> logEvidence_;
  SamplesSeries<std::pair<Beta,std::vector<mcmc>>> run_;
public:

  std::pair<double,double> logEvidence()const {return logEvidence_;}
  const SamplesSeries<std::pair<Beta,std::vector<mcmc>>>& samples()const
  {
    return run_;
  }

  Tempered_Evidence_Evaluation(SamplesSeries<std::pair<Beta,std::vector<mcmc>>> o)
    :run_(std::move(o)),logEvidence_{o.mean_var
                                     (&Evidence_pair,o.size()/2)}
  {


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
                       ,std::mt19937_64 &mt, double& dHd)
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
        std::cerr<<landa;
        os<<landa;
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
    std::cerr<<landa;
    os<<landa;
    landa=landa.partialReset();
    
    for (std::size_t i=0; i< 5; ++i)
      {
        n_steps_try(LM_Lik,lik,model,data,sDist,cDist,landa,mt,nsamples_per_par,os,startTime,timeOpt);
        landa.actualize();
        std::cerr<<landa;
        os<<landa;
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
        if (step(LM_Lik,lik,model,data,sDist,cDist,landa,dHd,mt,i,os,startTime,timeOpt))
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
   double dHd)
  {
    os<<" dHd "<<dHd<<" ";
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
   const std::vector<double>& dHd)
  {
    for (std::size_t  i=0; i<sLik.size(); ++i)
      {
        os<<sLik[i].beta<<"\t";
        if(accepts[i]) os<<"acc\t";
        else
          os<<"rej\t";
        os<<sLik[i].proposed.landa<<"\t";
        os<<" dHd "<<dHd[i]<<"\t ";
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
   ,std::mt19937_64& mt
   ,std::size_t i
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    bool out;
    LM_Lik.update_mcmc_step(lik,model,data,sDist,landa,sDist.beta);
    M_Matrix<double> c=sDist.proposed.sample(mt);
    cDist=LM_Lik.get_mcmc_step(lik,model,data,c,landa,sDist.beta);
    if (test::accept(sDist,cDist,mt,dHd))
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
    put(std::cout,sDist,cDist,dHd);
    
    os<<"n_steps::"<<i<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<sDist.beta;
    test::put(os,sDist,cDist);
    put(os,sDist,cDist,dHd);
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
   ,Adaptive_Beta& aBeta
   ,std::vector<double>& dHd
   , double p_Tjump
   ,std::vector<std::mt19937_64>& mt
   ,std::size_t isamples
   ,std::size_t isubSamples
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::vector<bool> out(sDist.size());
#pragma omp parallel for
    for(std::size_t i=0; i<sDist.size(); ++i)
      {
        AP landa=pars[i].sample(mt[i]);
        LM_Lik.update_mcmc_step(lik,model,data,sDist[i],landa,aBeta.getBeta().getValue()[i]);
        M_Matrix<double> c=sDist[i].proposed.sample(mt[i]);
        mcmc_step<pDist> cDist=LM_Lik.get_mcmc_step(lik,model,data,c,landa,aBeta.getBeta().getValue()[i]);
        if (test::accept(sDist[i],cDist,mt[i],dHd[i]))
          {
            pars[i].push_acceptance(landa,dHd[i]);
            sDist[i]=std::move(cDist);
            out[i] =true;
          }
        else
          {
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
    std::cout<<"isample::"<<isamples<<"\t"<<"isubSample::"<<isubSamples<<"\t"<<timeOpt<<"\t"<<timeIter_<<"Evidence\t"<<evidence<<"\n";
    
    os<<"isample::"<<isamples<<"\t"<<"isubSample::"<<isubSamples<<"\t"<<timeOpt<<"\t"<<timeIter_<<"Evidence\t"<<evidence<<"\n";
    
    put(os,sDist,out,dHd);
    put(std::cout,sDist,out,dHd);
    
    
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
                std::cout<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                aBeta.push_acceptance(i-1);
              }
            else
              aBeta.push_rejection(i-1);
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
    mystep cDist=LMLik.get_mcmc_step(lik,model,data,pinit,std::move(postL),landa,beta0);

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
  typedef Tempered_Evidence_Evaluation<mystep> myEvidence;

  static myEvidence *
  run
  (const MH<D,M,D_Lik,my_PropD,AP,propDist>& mcmc
   ,const propDist<D,M,D_Lik,my_PropD,AP>& LMLik
   ,const D_Lik<D,M>& lik
   ,const M<D>& model
   ,const D& data
   , const Adaptive_discrete<AP>& landa_Dist0
   ,const Adaptive_Beta& beta0
   , std::size_t nsamples
   ,std::size_t nskip
   ,double pTjump
   ,std::mt19937_64& mt
   ,std::ofstream& os
   ,const std::chrono::steady_clock::time_point& startTime
   ,double& timeOpt)
  {
    Adaptive_Beta beta(beta0);
    auto n=beta.size();
    std::vector<mystep> sDists(n);
    std::vector<Adaptive_discrete<AP>> pars(n,landa_Dist0);
    std::vector<double> dHd(beta.size());
    std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
    std::vector<std::mt19937_64> mts(n);
    for (std::size_t i=0; i<n; ++i)
      mts[i].seed(useed(mt));
    SamplesSeries<typename std::pair<Beta,std::vector<mystep>>> o(nsamples);
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
        sDists[i]=LMLik.get_mcmc_step(lik,model,data,pinit,std::move(postL),landa,beta0);
        LMLik.update_mcmc_step(lik,model,data,sDists[i],landa,beta0);
      }

    while (!o.full())
      {
        for (std::size_t i=0; i<n;++i)
          pars[i].actualize();

        for( std::size_t i=0;i<nskip; ++i)
          {
            mcmc.tempered_step
                (LMLik,lik,model,data,sDists,pars,beta,dHd,pTjump,mts,o.size(),i,os,startTime,timeOpt);
          }

        o.push_back({beta.getBeta(),sDists});
        for (std::size_t i=0; i<n;++i)
          {
            pars[i].actualize();
            std::cerr<<pars[i];
            os<<pars[i];
          }
        std::cout<<beta;
        os<<beta;
        beta.actualize(nskip);


      }
    auto out=new myEvidence(o);
    std::cout<<"LogEvidence= "<<out->logEvidence()<<"\n";


    return out;
  }




};







#endif // EVIDENCE_H
