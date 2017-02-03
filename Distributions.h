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

  double operator()(const M_Matrix<double>& x)const
  {
    return logP(x);
  }
  M_Matrix<double> sample(std::mt19937_64& mt)const ;
  M_Matrix<double> operator()(std::mt19937_64& mt)const
  {
    return sample(mt);
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
  os<<"\nCovInv\n"<<x.CovInv();
  os<<"\nChol\n"<<x.Chol();
  os<<"\nlogDetCov\n"<<x.logDetCov();
  return os;
}


#endif // DISTRIBUTIONS_H


