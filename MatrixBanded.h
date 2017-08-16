#ifndef MATRIXBANDED_H
#define MATRIXBANDED_H

#include<vector>
#include<cassert>
#include "Matrix.h"
extern "C" void dgbsv_
( int* 	N,
  int*  	KL,
  int*  	KU,
  int*  	NRHS,
  double * AB,
  int*  	LDAB,
  int * IPIV,
  double *  B,
  int*  	LDB,
  int*  	INFO
  );


class MatrixBanded
{
public:
  MatrixBanded(std::size_t n, std::size_t KL, std::size_t KU)
    : n_{n},KL_(KL),KU_{KU},m_(2*KL+KU+1),d_(n*(2*KL+KU+1))
  {}
  double& operator()(std::size_t i, std::size_t j)
  {
    assert(std::max(int(j)-int(KU_),0)<=int(i));
    assert(int(i)<=std::min(int(n_)-1,int(j+KL_)));
     return d_[KL_+KU_+i-j+j*m_];
  }

  const double& operator()(std::size_t i, std::size_t j) const
  {
    assert(std::max(int(j)-int(KU_),0)<=int(i));
    assert(int(i)<=std::min(int(n_)-1,int(j+KL_)));
    return d_[KL_+KU_+i-j+j*m_];
  }

  std::vector<double> solve(const std::vector<double>& A)const
  {
    std::vector<double> AB(d_);
    std::vector<double> B(A);
    std::vector<int> IPIV(n_);
    int N=n_;
    int KL=KL_;
    int KU=KU_;
    int NRHS=1;
    int LDAB=m_;
    int LDB=n_;
    int INFO=0;
    dgbsv_(&N,&KL,&KU,&NRHS,&AB[0],&LDAB,&IPIV[0],&B[0],&LDB,&INFO);

    //auto C=(*this)*B;
    //auto D=C-A;
    return B;
  }

  MatrixBanded& operator*=(double dt)
  {
    for (std::size_t i=0; i<d_.size(); ++i)
      d_[i]*=dt;
    return *this;
  }
  std::size_t n()const {return n_;}
  std::size_t num_up()const {return KU_;}
  std::size_t num_low()const {return KL_;}

  std::vector<double> operator*(const std::vector<double> x)const
  {
    assert(x.size()==n_);
    std::vector<double> out(n_,0);
    for (std::size_t i=0; i<n_; ++i)
      {
        for (std::size_t j=std::max(n_+i-KL_,n_)-n_; j<std::min(i+KU_,n_); ++j)
          out[i]+=(*this)(i,j)*x[j];
      }
    return out;
  }


private:
  std::size_t  n_;std::size_t KL_;std::size_t KU_;
  std::size_t m_;
  std::vector<double> d_;
};



#endif // MATRIXBANDED_H
