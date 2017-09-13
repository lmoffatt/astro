#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "Matrix.h"

#include<random>

#ifndef PI____
#define PI____
const double PI  =3.141592653589793238463;
#endif



class MultivariateGaussian
{
public:



  MultivariateGaussian(const M_Matrix<double>& mean,
                       const M_Matrix<double>& cov);


  MultivariateGaussian(const M_Matrix<double>& mean,
                       const M_Matrix<double>& cov,
                       const M_Matrix<double>& covInv
                       );

  MultivariateGaussian();

  double logP(const M_Matrix<double>& x)const ;
  double P(const M_Matrix<double>& x)const ;

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

  double chi2(const M_Matrix<double>& x)const;

  double operator()(const M_Matrix<double>& x)const
  {
    return logP(x);
  }
  M_Matrix<double> sample(std::mt19937_64& mt)const ;
  M_Matrix<double> operator()(std::mt19937_64& mt)const
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


  const M_Matrix<double>& Mean()const;
  M_Matrix<double> Cov()const;
  const M_Matrix<double>& CovInv()const;


  const M_Matrix<double>& Chol()const{return cho_cov_;}

  double logDetCov()const {return logDetCov_;}

  virtual M_Matrix<double> logPGradient(const M_Matrix<double>& x)const
  {
    return M_Matrix<double>(1u,x.size(),(x-Mean())*Cov());
  }
  virtual M_Matrix<double> logPHessian(const M_Matrix<double>& )const
  {
    return Cov();
  }


  MultivariateGaussian(const MultivariateGaussian& other)=default;
  MultivariateGaussian& operator=(const MultivariateGaussian& other)=default;


  ~MultivariateGaussian();

  std::size_t size() const;
private:
  M_Matrix<double> mean_;
  M_Matrix<double> cov_;
  M_Matrix<double> covinv_;
  M_Matrix<double> cho_cov_;
  double logDetCov_;

  // Distribution interface
public:
  virtual std::__cxx11::string myClass() const
  {
    return "MultivariateGaussian";
  }
};



inline
std::ostream& operator<<(std::ostream& os,const MultivariateGaussian& x)
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
            e.second.Parameters().first+=(1.-f)*0.5;
            e.second.Parameters().second*=f;
            e.second.Parameters().second+=(1.-f)*0.5;

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



#endif // DISTRIBUTIONS_H


