#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "Matrix.h"

#include<random>

#ifndef PI____
#define PI____
const double PI  =3.141592653589793238463;
#endif


template<class P>
struct complement_prob
{
  complement_prob(const P& p):p_{p}{}
  template<typename... Ts>
  double operator()(Ts... xs)const
  {
    return 1.0-p_(xs...);
  }
private:
  const P& p_;
};

template<class P>
complement_prob<P>
Complement_prob(const P& p){return complement_prob<P>(p);}


template<class P>
struct log_of
{
  log_of(const P& p):p_{p}{}
  template<typename... Ts>
  double operator()(Ts... xs)const
  {
    return std::log(p_(xs...));
  }
private:
  const P& p_;
};


template<class P>
log_of<P> Log_of(const P& p){return log_of<P>(p);}

template<class P>
struct exp_of
{
  exp_of(const P& p):p_{p}{}
  template<typename... Ts>
  double operator()(Ts... xs)const
  {
    return std::exp(p_(xs...));
  }
private:
  const P& p_;
};


template<class E>
class MultivariateGaussian
{
public:

  MultivariateGaussian(const M_Matrix<E> &mean,
                       const M_Matrix<E>& cov):
    mean_(mean),
    cov_(cov),
    covinv_(inv(cov).first),
    cho_cov_(chol(cov,"lower")),
    logDetCov_(logDiagProduct(cho_cov_))
  {}

  MultivariateGaussian(const M_Matrix<E>&mean
                       , const M_Matrix<E> &cov
                       , const M_Matrix<E> &covInv):
    mean_(mean),
    cov_(cov),
    covinv_(covInv),
    cho_cov_(chol(cov,"lower").first),
    logDetCov_(logDiagProduct(cho_cov_))
  {}






  MultivariateGaussian()=default;

  double logP(const M_Matrix<E> &x) const
  {
    if (mean_.size()>0)

      return -0.5*size()*log(PI)-logDetCov()-chi2(x);
    else return std::numeric_limits<double>::quiet_NaN();
  }
  double P(const M_Matrix<E>& x)const
  {
    return exp(logP(x));
  }

  void autoTest(std::mt19937_64& mt,std::size_t n)const
  {
    std::cerr<<"chi test n="<<size()<<" chis\n";
    double chisum=0;
    double chisqr=0;
    for (std::size_t i=0; i<n; ++i)
      {
        auto s=sample(mt);
        auto chi=chi2(s);
        //  std::cerr<<chi<<" ";
        chisum+=chi;
        chisqr+=chi*chi;
      }
    chisum/=n;
    chisqr-=n*chisum*chisum;
    chisqr/=(n-1);
    std::cerr<<"\n chimean="<<chisum<<" chisqr="<<chisqr;
  }

  double chi2(const M_Matrix<E> &x) const
  {
    if (!mean_.empty())
      return 0.5*xTSigmaX(x-mean_,covinv_);
    else return std::numeric_limits<double>::quiet_NaN();
  }

  double operator()(const M_Matrix<E>& x)const
  {
    return logP(x);
  }
  M_Matrix<E> sample(std::mt19937_64& mt)const
  {
    M_Matrix<E> r;
    std::normal_distribution<> normal;
    if (this->size()>0)
      {
        auto z=Rand(mean_,normal,mt);
        r=mean_+multTransp(z,cho_cov_);
      }
    return r;
  }
  M_Matrix<E> operator()(std::mt19937_64& mt)const
  {
    return sample(mt);
  }
  friend
  std::istream& operator>>(std::istream& is, MultivariateGaussian& x)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>x.mean_;
    std::getline(is,line);
    std::getline(is,line);
    is>>x.cov_;
    return is;
  }


  const M_Matrix<E>& Mean()const
  {
    return mean_;
  }
  M_Matrix<E> Cov()const
  {
    return inv(covinv_).first;
  }

  const M_Matrix<E>& CovInv()const
  {
    return covinv_;
  }




  const M_Matrix<E>& Chol()const{return cho_cov_;}

  double logDetCov()const {return logDetCov_;}

  virtual M_Matrix<E> logPGradient(const M_Matrix<E>& x)const
  {
    return M_Matrix<E>(1u,x.size(),(x-Mean())*Cov());
  }
  virtual M_Matrix<E> logPHessian(const M_Matrix<E>& )const
  {
    return Cov();
  }



  std::size_t size()const
  {
    return mean_.size();
  }



  MultivariateGaussian(const MultivariateGaussian& other)=default;
  MultivariateGaussian& operator=(const MultivariateGaussian& other)=default;


  ~MultivariateGaussian(){};

private:
  M_Matrix<E> mean_;
  M_Matrix<E> cov_;
  M_Matrix<E> covinv_;
  M_Matrix<E> cho_cov_;
  double logDetCov_;

  // Distribution interface
public:
  virtual std::__cxx11::string myClass() const
  {
    return "MultivariateGaussian";
  }
};



template<class E>
std::ostream& operator<<(std::ostream& os,const MultivariateGaussian<E>& x)
{
  os<<"\nMean\n"<<x.Mean();
  os<<"\nCov\n"<<x.Cov();
  return os;
}
class Beta_Distribution
{
public:
  Beta_Distribution(const std::pair<double,double>& v):
    a_(v){}

  Beta_Distribution():a_{0.5,0.5}{}

  std::pair<double,double>& Parameters(){return a_;}

  std::pair<double,double>const & Parameters()const {return a_;}

  double count()const {return a_.first+a_.second;}

  static Beta_Distribution UniformPrior()
  {
    return Beta_Distribution({1.0,1.0});
  }
  static Beta_Distribution UnInformativePrior()
  {
    return Beta_Distribution({0.5,0.5});
  }


  double p()const {return a_.first/(a_.first+a_.second);}

  void push_accept()
  {
    ++a_.first;
  }
  void push_reject()
  {
    ++a_.second;
  }


  double operator()(std::mt19937_64& mt)
  {
    std::gamma_distribution<double> ga(a_.first,2.0);
    std::gamma_distribution<double> gb(a_.second,2.0);

    double a=ga(mt);
    double b=gb(mt);
    return a/(a+b);
  }

  friend
  std::istream& operator>>(std::istream& is, Beta_Distribution& me)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>me.a_;
    return is;
  }

  friend
  std::ostream& operator<<(std::ostream& os,const Beta_Distribution& me)
  {
    os<<"\nBeta alfa beta\n";
    os<<me.a_;
    return os;
  }

private:
  std::pair<double,double> a_;

};



class Dirichlet_Distribution
{
public:
  Dirichlet_Distribution(const std::vector<double>& v):
    a_(v){}

  std::vector<double> operator()(std::mt19937_64& mt)
  {
    std::vector<double> out(a_.size());
    double sum=0;
    for (std::size_t i=0; i<a_.size(); ++i)
      {
        std::gamma_distribution<double> g(a_[i]);
        out[i]=g(mt);
        sum+=out[i];
      }
    for (auto& o:out)
      o/=sum;
    return out;
  }

  std::istream& operator>>(std::istream& is)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>a_;
    return is;
  }

  std::ostream& operator<<(std::ostream& os)const
  {
    os<<"\nDirichlet alfas\n";
    os<<a_;
    return os;
  }

private:
  std::vector<double> a_;

};

template<typename T>
std::map<T,double> operator*(std::map<T,double> in, double x)
{
  std::map<T,double> out(in);
  for (auto& e:out) e.second*=x;
  return out;
}

template<typename T>
class Dirichlet_map
{
  static Dirichlet_map setPrior(const std::map<T,double> a, double p)
  {
    std::map<T,double> o(a);
    for (auto& e:o) e.second=p;
    return Dirichlet_map(o);
  }

public:
  Dirichlet_map(const std::map<T,double> a):a_(a){}

  static Dirichlet_map UniformPrior(const std::map<T,double> a)
  {
    return setPrior(a,1.0);
  }

  static Dirichlet_map UninformativePrior(const std::map<T,double> a)
  {
    return setPrior(a,0.5);
  }




  Dirichlet_map()=default;
  std::size_t size()const {return a_.size();}
  double count()const {
    double sum=0;
    for (auto e:a_) sum+=e.second;
    return sum;
  }
  std::map<T,double> operator()(std::mt19937_64& mt)
  {
    std::map<T,double> out;
    double sum=0;

    for (auto it=a_.begin(); it!=a_.end(); ++it)
      {
        std::gamma_distribution<double> g(it->second);
        out[it->first]=g(mt);
        sum+=out[it->first];
      }
    for (auto& o:out)
      o.second/=sum;
    return out;
  }

  Dirichlet_map& operator+=(const Dirichlet_map& other)
  {
    for (auto& e:a_)
      {
        auto it=other.a_.find(e.first);
        if (it!=other.a_.end())
          e.second+=it->second;
      }
    return *this;
  }

  Dirichlet_map operator+(const Dirichlet_map& other)const
  {
    Dirichlet_map out(a_);
    out+=other;
    return out;
  }

  std::map<T,double> p()const
  {
    std::map<T,double> out(a_);
    double sum=count();
    for (auto& e:out) e.second/=sum;
    return out;
  }

  friend
  std::istream& operator>>(std::istream& is, Dirichlet_map& me)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>me.a_;
    return is;
  }

  friend
  std::ostream& operator<<(std::ostream& os,const Dirichlet_map& me)
  {
    os<<"\nDirichlet alfas\n";
    os<<me.a_;
    return os;
  }

private:
  std::map<T,double> a_;
};


template<typename T>
class Beta_map
{

public:
  Beta_map(const std::map<T,Beta_Distribution> a):a_(a){}

  void reduce(double nmax)
  {
    double f=nmax/count();
    if (f<1.0)
      {
        for (auto& e:a_)
          {
            e.second.Parameters().first*=f;
            e.second.Parameters().second*=f;

          }
      }
  }

  static Beta_map UniformPrior(const std::map<T,double> a)
  {
    std::map<T,Beta_Distribution> o;
    for (auto& e:a)
      o[e.first]=Beta_Distribution::UniformPrior();
    return Beta_map(o);
  }


  static Beta_map UnInformativePrior(const std::map<T,double> a)
  {
    std::map<T,Beta_Distribution> o;
    for (auto& e:a)
      o[e.first]=Beta_Distribution::UnInformativePrior();
    return Beta_map(o);
  }




  Beta_map()=default;

  std::size_t size()const {return a_.size();}

  double count()const {
    double sum=0;
    for (auto& e:a_) sum+=e.second.count();
    return sum;
  }


  std::map<T,double> operator()(std::mt19937_64& mt)
  {
    std::map<T,double> out;
    double sum=0;

    for (auto it=a_.begin(); it!=a_.end(); ++it)
      {
        std::gamma_distribution<double> g(it->second);
        out[it->first]=g(mt);
        sum+=out[it->first];
      }
    for (auto& o:out)
      o.second/=sum;
    return out;
  }

  Beta_map& operator+=(const Beta_map& other)
  {
    for (auto& e:a_)
      {
        auto it=other.a_.find(e.first);
        if (it!=other.a_.end())
          e.second.Parameters()+=it->second.Parameters();
      }
    return *this;
  }

  Beta_map operator+(const Beta_map& other)const
  {
    Beta_map out(a_);
    out+=other;
    return out;
  }

  std::map<T,double> p()const
  {
    std::map<T,double> out;
    for (auto& e:a_)
      out[e.first]=e.second.p();
    return out;
  }

  Beta_Distribution& operator[](const T& x)
  {
    return a_[x];
  }

  Beta_Distribution operator[](const T& x)const
  {
    auto it=a_.find(x);
    if (it!=a_.end())
      return it.second;
    else
      return {};
  }


  friend
  std::istream& operator>>(std::istream& is, Beta_map& me)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>me.a_;
    return is;
  }

  friend
  std::ostream& operator<<(std::ostream& os,const Beta_map& me)
  {
    os<<"\nBeta map\n";
    os<<me.a_;
    return os;
  }

private:
  std::map<T,Beta_Distribution> a_;
};

template <typename T>
T sample_rev_map(const std::map<double,T>& reverse_prior,std::mt19937_64& mt)
{
  std::uniform_real_distribution<> u;
  double r= u(mt);
  auto it=reverse_prior.lower_bound(r);
  return it->second;

}

template <class T>
std::pair<std::map<T,double>,double>
normalize_map(const std::map<T,double>& unnormalized_map)
{
  std::map<T,double> out(unnormalized_map);
  double Evidence=0;
  for (auto& e:out)
    {
      Evidence+=e.second;
    }
  for (auto& e:out) e.second/=Evidence;

  return {out,Evidence};
}
template <class T>
std::pair<std::map<T,double>,double>
logLik_to_p(const std::map<T,double>& logLikelihoods)
{
  std::map<T,double> out(logLikelihoods);
  double Evidence=0;
  double maxlog=out.begin()->second;
  for (auto& e:out)
    {
      if (e.second>maxlog) maxlog=e.second;
    }
  for (auto& e:out)
    {
      e.second=std::exp(e.second-maxlog);
      Evidence+=e.second;
    }

  for (auto& e:out) e.second/=Evidence;

  return {out,Evidence};
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

template <class T>
std::map<double,T>
cumulative_reverse_logmap(const std::map<T,double>& normalized_logmap)
{
  std::map<double,T> out;
  double sump=0;
  for (auto& e:normalized_logmap)
    {
      sump+=std::exp(e.second);
      out[sump]=e.first;
    }
  return out;
}


template<typename T>
class Probability_map
{
public:
  T operator()(std::mt19937_64& mt)const
  {
    return sample_rev_map(rev_,mt);
  }

  const std::map<T,double>& p() const
  {
    return p_;
  }

  Probability_map(const std::map<T,double>& myNormalized_map, double nsamples)
    :
      p_{myNormalized_map},rev_{cumulative_reverse_map(p_)}, nsamples_(nsamples)
  {

  }

  template<template<typename...>class V>
  Probability_map(const V<T>& x):p_(Uniform(x)),rev_(cumulative_reverse_map(p_)), nsamples_(0)
  {}
  Probability_map()=default;

  template<template<typename...>class V>
  static
  std::map<T,double> Uniform(const V<T>& x)
  {
    std::map<T,double> out;
    std::size_t n=x.size();
    double p=1.0/n;
    for (std::size_t i=0; i<n;++i )
      out[x[i]]+=p;
    return out;
  }


  void reduce(double nmax)
  {
    double f=nsamples_/nmax;
    if (f<1.0)
      {
        auto o=p_;
        for (auto& e:o) e.second=std::pow(e.second,f);
        *this=normalize(o,nmax).first;
      }

  }
  double nsamples()const {return nsamples_;}

  static std::pair<Probability_map,double> normalize(const std::map<T,double>& myposterior , double nsamples)
  {
    auto out= normalize_map(myposterior);
    return {Probability_map(out.first,nsamples),out.second};
  }

  friend
  std::istream& operator>>(std::istream& is, Probability_map& me)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>me.p_;
    me.rev_=cumulative_reverse_map(me.p_);
    std::getline(is,line);
    std::getline(is,line);
    is>>me.nsamples_;
    return is;
  }

  friend
  std::ostream& operator<<(std::ostream& os,const Probability_map& me)
  {
    os<<"\nProbability map\n";
    os<<me.p_;
    os<<"\nNumSamples\n";
    os<<me.nsamples();
    return os;
  }


private:

  std::map<T,double> p_;
  std::map<double,T> rev_;
  double nsamples_;
};

template<typename T>
class logLikelihood_map
{
public:
  T operator()(std::mt19937_64& mt)
  {
    return sample_rev_map(rev_,mt);
  }

  const std::map<T,double>& logLik() const
  {
    return logLik_;
  }

  std::map<T,double>const & p()const
  {
    return p_;
  }

  logLikelihood_map(const std::map<T,double>& mylogLikelihood_Map, double nsamples)
    :
      logLik_{mylogLikelihood_Map}
  {
    auto p=logLik_to_p(mylogLikelihood_Map);
    p_=std::move(p.first);
    Evidence_=p.second;
    rev_=cumulative_reverse_map(p_);
    nsamples_=nsamples;
  }

  logLikelihood_map()=default;
  void reduce(double nmax)
  {
    double f=nsamples_/nmax;
    if (f<1.0)
      {
        auto o=logLik_;
        for (auto& e:o) e.second*=f;
        *this=logLikelihood_map(o,nmax);
      }
  }


  double nsamples()const {return nsamples_;}
  friend
  std::istream& operator>>(std::istream& is, logLikelihood_map& me)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>me.logLik_;
    me.p_=logLik_to_p(me.logLik_).first;
    me.rev_=cumulative_reverse_map(me.p_);
    std::getline(is,line);
    std::getline(is,line);
    is>>me.nsamples_;
    return is;
  }

  friend
  std::ostream& operator<<(std::ostream& os,const logLikelihood_map& me)
  {
    os<<"\nlogLikelihood map\n";
    os<<me.logLik_;
    os<<"\nNumSamples\n";
    os<<me.nsamples();
    return os;
  }

private:
  std::map<T,double> logLik_;
  std::map<T,double>p_;
  std::map<double,T> rev_;
  double Evidence_;
  double nsamples_;

};




template<class Likelihood, class Data, typename T>
std::pair<Probability_map<T>,double>
Bayes_rule(const Likelihood& lik, const Data& data, const Probability_map<T>& prior)
{
  auto p=prior.p();
  for (auto& e:p)
    {
      double l=lik(e.first,data);
      e.second*=l;
    }
  double nsamples=prior.nsamples()+1;
  return Probability_map<T>::normalize(p,nsamples);
}

template<class logLikelihood, class Data, typename T>
logLikelihood_map<T>
logBayes_rule(const logLikelihood& loglik, const Data& data, const logLikelihood_map<T>& prior)
{
  auto logP=prior.logLik();
  for (auto& e:logP)
    {
      double logL=loglik(e.first,data);
      e.second+=logL;
    }
  double nsamples=prior.nsamples()+1;
  return logLikelihood_map<T>(logP,nsamples);
}


template< typename T, class F, class Likelihood,class P_map>
double Expectance(const F& f, const Likelihood& lik,const P_map& pm, const T& landa  )
{
  auto p=pm;
  double sum=0;
  for (auto& e:p)
    sum+=e.second*lik(e.first,landa)*f(landa);
  return sum;
}






#endif // DISTRIBUTIONS_H


