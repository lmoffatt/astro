#ifndef MATRIX_H
#define MATRIX_H
/*!
 * @file Matrix.h


 */


#include <limits> // for std::numeric_limits
#include <cstddef> //  for std::size_t
#include <ostream>
#include <istream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cassert>
#include <random>

#include <iterator>     // std::iterator, std::input_iterator_tag

#include <iostream>

#include "mySerializer.h"
#include "myOrderOperators.h"

inline double logit(double x){return std::log(x/(1.0-x));}



inline std::pair<double,double> logit(double x,double sd){
  return {std::log(x/(1.0-x)),sd/(x*(1.0-x))};}


inline double logistic(double x){return 1.0/(1.0+std::exp(-x));}


template<typename T1, typename T2>
std::pair<T1,T2>& operator+=(std::pair<T1,T2>& x, const std::pair<T1,T2>& other)
{
  x.first+=other.first;
  x.second+=other.second;
  return x;
}

inline double average(double x, double y){return 0.5*(x+y);}

inline double sqr(double x){return x*x;}










namespace
{
  extern "C" void dgetrf_(int *M,
                          int* N,
                          double *A,
                          int* LDA,
                          int* IPIV,
                          int * INFO );

  extern "C" void dgetri_(int* n,
                          double *B,
                          int* dla,
                          int* ipiv,
                          double* work1,
                          int* lwork,
                          int* info);

  extern "C" void dgemm_(char * 	TRANSA,
                         char * 	TRANSB,
                         int * 	M,
                         int * 	N,
                         int * 	K,
                         double * ALPHA,
                         double * A,
                         int * 	LDA,
                         double * B,
                         int * 	LDB,
                         double * BETA,
                         double * C,
                         int * 	LDC
                         );


}



template <typename T>
std::vector<T> operator -(const std::vector<T>& x,const std::vector<T>& y)
{
  std::vector<T> out(x.size());
  for (std::size_t i=0; i<x.size(); ++i)
    out[i]=x[i]-y[i];
  return out;
}

template <typename T>
std::vector<T> operator +(const std::vector<T>& x,const std::vector<T>& y)
{
  std::vector<T> out(x.size());
  for (std::size_t i=0; i<x.size(); ++i)
    out[i]=x[i]+y[i];
  return out;
}


template <typename T>
std::vector<T> elemMult(const std::vector<T>& x,const std::vector<T>& y)
{
  std::vector<T> out(x.size());
  for (std::size_t i=0; i<x.size(); ++i)
    out[i]=x[i]*y[i];
  return out;
}

template <typename T>
std::vector<T> elemDiv(const std::vector<T>& x,const std::vector<T>& y)
{
  std::vector<T> out(x.size());
  for (std::size_t i=0; i<x.size(); ++i)
    out[i]=x[i]/y[i];
  return out;
}



template<typename T>
class M_Matrix
{
public:

  class MyConstIterator;
  class MyIterator : public std::iterator<std::input_iterator_tag, T>
  {
    std::size_t i_;
    std::size_t j_;
    M_Matrix<T>& m;
  public:
    friend class MyConstIterator;
    std::size_t iRow()const{return i_;}
    std::size_t jCol()const{return j_;}

    MyIterator(M_Matrix<T>& x) :m(x),i_(0),j_(0) {}
    MyIterator(M_Matrix<T>& x, std::size_t i, std::size_t j) :m(x),i_(i),j_(j) {}

    MyIterator(MyIterator& mit) : m(mit.m),i_(mit.i_),j_(mit.j_) {}
    MyIterator& operator++()
    {
      ++j_;
      if (j_>=ncols(m))
        {
          j_=0;
          ++i_;
        }
      return *this;}
    MyIterator operator++(int)
    {MyIterator tmp(*this); operator++(); return tmp;}
    bool operator==(const MyIterator& rhs)
    {
      if (i_!=rhs.i_)
        return false;
      else if (j_!=rhs.j_)
        return false;
      else return true;
    }
    bool operator!=(const MyIterator& rhs) {return ! (*this==rhs);}
    T& operator*() {return m(i_,j_);}
    const T& operator*() const {return m(i_,j_);}
  };


  class MyConstIterator : public std::iterator<std::input_iterator_tag, T>
  {
    std::size_t i_;
    std::size_t j_;
    const M_Matrix<T>& m;
  public:
    std::size_t iRow()const{return i_;}
    std::size_t jCol()const{return j_;}

    MyConstIterator(const M_Matrix<T>& x) :m(x),i_(0),j_(0) {}
    MyConstIterator(const M_Matrix<T>& x, std::size_t i, std::size_t j) :m(x),i_(i),j_(j) {}

    MyConstIterator(const MyConstIterator& mit) : m(mit.m),i_(mit.i_),j_(mit.j_) {}
    MyConstIterator(const MyIterator& mit) : m(mit.m),i_(mit.i_),j_(mit.j_) {}
    MyConstIterator& operator++()
    {
      ++j_;
      if (j_>=m.ncols())
        {
          j_=0;
          ++i_;
        }
      return *this;}
    MyConstIterator operator++(int)
    {MyConstIterator tmp(*this); operator++(); return tmp;}
    bool operator==(const MyConstIterator& rhs)
    {
      if (i_!=rhs.i_)
        return false;
      else if (j_!=rhs.j_)
        return false;
      else return true;
    }
    bool operator!=(const MyConstIterator& rhs) {return ! (*this==rhs);}
    const T& operator*() const {return m(i_,j_);}
  };



  typedef  MyIterator iterator;
  typedef  MyConstIterator const_iterator;


  iterator begin()
  {
    MyIterator out(*this);
    return out;
  }

  iterator end()
  {
    MyIterator out(*this,0,nrows(*this));
    return out;
  }

  const_iterator begin()const
  {
    MyConstIterator out(*this);
    return out;
  }

  const_iterator end() const
  {
    MyConstIterator out(*this,nrows(),0);
    return out;
  }


  M_Matrix():
    _nrows(std::size_t(0)),
    _ncols(std::size_t(0)),
    _ncells(std::size_t(0)),
    _data(0) {}  //default constructor


  M_Matrix (const M_Matrix<T> & sample)
    : _nrows(sample._nrows),
      _ncols(sample._ncols),
      _ncells(sample._ncells),
      //   _data(_ncells>0?new T[x.size()]:0)
      _data(sample._data)
  {
    /*  for (std::size_t i = 0; i < size(sample); ++i)
      (*this)[i] = sample[i];
 */
  }

  M_Matrix (std::size_t nrows,std::size_t ncols)
    : _nrows(nrows),
      _ncols(ncols),
      _ncells(nrows*ncols),
      //   _data(_ncells>0?new T[x.size()]:0)
      _data(nrows*ncols)
  {
    /*  for (std::size_t i = 0; i < size(sample); ++i)
      (*this)[i] = sample[i];
 */
  }


  template<typename S>
  M_Matrix (const M_Matrix<S> & sample)
    : _nrows(sample.nrows()),
      _ncols(sample.ncols()),
      _ncells(sample.size()),
      _data(sample.size())

  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      (*this)[i] = sample[i];

  }
  M_Matrix (const std::vector<std::vector<T>> & sample)
    : _nrows(sample.size()),
      _ncols(sample[0].size()),
      _ncells(sample.size()*sample[0].size())
    ,_data(sample.size()*sample[0].size())
  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      for (std::size_t j=0; j<sample[0].size(); ++j)
        (*this)(i,j) = sample[i][j];

  }


  template<typename S>
  M_Matrix (const std::vector<std::vector<S>> & sample)
    : _nrows(sample.size()),
      _ncols(sample[0].size()),
      _ncells(sample.size()*sample[0].size())
    ,_data(sample.size()*sample[0].size())
  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      (*this)[i] = sample[i];

  }



  M_Matrix<T>& operator=(const M_Matrix<T>& x)=default;




  ~M_Matrix()
  {
    //  if (_data>0)
    //	delete [] _data;
  }



  M_Matrix(std::size_t nrows_,std::size_t ncols_, std::vector<T> data):
    _nrows(nrows_),
    _ncols(ncols_),
    _ncells(nrows_*ncols_),
    //   _data(new T[_ncells])
    _data(_ncells)
  {

    for (size_t i=0; i<_ncells; ++i)
      _data[i]=data[i];
  }

  M_Matrix(std::size_t nrows_,std::size_t ncols_, const M_Matrix<T>& data):
    _nrows(nrows_),
    _ncols(ncols_),
    _ncells(nrows_*ncols_),
    //   _data(new T[_ncells])
    _data(_ncells)
  {

    for (size_t i=0; i<_ncells; ++i)
      _data[i]=data[i];
  }

  template <typename S>  M_Matrix(std::size_t nrows_,std::size_t ncols_, const M_Matrix<S>& data):
    _nrows(nrows_),
    _ncols(ncols_),
    _ncells(nrows_*ncols_),
    //   _data(new T[_ncells])
    _data(_ncells)
  {

    for (size_t i=0; i<_ncells; ++i)
      _data[i]=data[i];
  }


  M_Matrix(std::size_t nrows_,std::size_t ncols_,T data):
    _nrows(nrows_),
    _ncols(ncols_),
    _ncells(nrows_*ncols_),
    //   _data(new T[_ncells])
    _data(_ncells,data)
  {

  }


  template<class F>
  M_Matrix<T>
  apply(const F& f)
  {
    M_Matrix<T> out(nrows(),ncols());
    for (std::size_t i=0; i<size(); ++i)
      out[i]=f((*this)[i]);
    return out;
  }



  M_Matrix<T>& operator=(T X)
  {
    for (std::size_t i=0; i<_ncells; ++i) _data[i]=X;
    return *this;
  }



  size_t size()const
  {
    return this->_ncells;
  }

  size_t nrows()const
  {
    return _nrows;
  }

  size_t ncols()const
  {
    return _ncols;
  }

  T& operator[](std::size_t n)
  {
    return _data[n];
  }
  bool empty()const
  {
    return _data.empty();
  }

  const T& operator[](std::size_t n) const
  {
    return _data[n];
  }


  T&  operator() (std::size_t i,std::size_t j)
  {
    return (*this)[i*_ncols+j];

  }

  const  T&  operator() (std::size_t i,std::size_t j) const
  {
    return (*this)[i*_ncols+j];

  }


  M_Matrix<T>&  operator() (std::size_t i,
                            const std::string /*dummy*/,
                            const M_Matrix<T>& newValues);

  M_Matrix<T>&  operator() (const std::string /*dummy*/,
                            std::size_t jColumn,
                            const M_Matrix<T>& newValues);

  M_Matrix<T>  operator() (std::size_t iRow,
                           std::string /*dummy*/) const;


  M_Matrix<T>  operator() (std::string /*dummy*/,
                           std::size_t jColumn) const;

  void clear()
  {
    _nrows=0;
    _ncols=0;
    _ncells=0;
    _data.clear();
  }


  std::vector<T> toVector()const
  {
    return _data;
  }

  M_Matrix<T> toVector_of_Rows()const
  {
    return M_Matrix<T>(size(),1,_data);
  }
  M_Matrix<T> toVector_of_Cols()const
  {
    return M_Matrix<T>(1,size(),_data);
  }


  std::vector<std::vector<T>> toMatrix()const
  {
    std::vector<std::vector<T>> out(nrows(),std::vector<T>(ncols()));
    for (std::size_t i=0;i<nrows();++i)
      for (std::size_t j=0;j<ncols();++j)
        out[i][j]=(*this)(i,j);
    return out;
  }


private:
  std::size_t          _nrows;    /**< number of rows */
  std::size_t          _ncols;    /**< number of columns */
  std::size_t          _ncells;   /**< _nows*_ncells  */

  /** internal data
          @remarks can be a pointer or a vector
          @remarks now is a vector for debugging purposes
          */
  //	T                      *_data; /**< pointer to the data  */
  std::vector<T>        _data;

};



/**
   Returns a custom sized Matrix filled with ones
  @post (ones(n,m))(i,j)==T(1)
   */
template<typename T>
M_Matrix<T>  ones(size_t nrows_, size_t ncols_)
{
  M_Matrix<T> A(nrows_,ncols_);
  for (size_t i=0; i<A.size(); ++i)
    A[i]=1;
  return A;
}

/**
   Matrix filled with ones with the shape of the provided Matrix
   @post nrows(ones(x))==x.nrows()   ncols(ones(x))=x.ncols()
   */
template<typename T>
M_Matrix<T>  ones(const M_Matrix<T>& x)
{
  M_Matrix<T> A(x.nrows(),x.ncols());
  for (size_t i=0; i<A.size(); ++i)
    A[i]=1;

  return A;
}
/**
    Custom sized Matrix filled with zeros
  */
template<typename T>
M_Matrix<T>  zeros(std::size_t nrows_, std::size_t ncols_)
{
  M_Matrix<T> A(nrows_,ncols_);
  for (std::size_t i=0; i<A.size(); ++i)
    A[i]=0;
  return A;
}

/**
  Matrix of zeros with the shape of the provided Matrix
   @post nrows(ones(x))==x.nrows()   ncols(ones(x))=x.ncols()

  */
template<typename T>
M_Matrix<T>  zeros(const M_Matrix<T>& x)
{
  M_Matrix<T> A(x.nrows(),x.ncols());
  for (std::size_t i=0; i<A.size(); ++i)
    A[i]=0;

  return A;
}

/**
   Identity Matrix of the specified size
  */
template<typename T>
M_Matrix<T> eye(std::size_t n)
{
  M_Matrix<T> A=zeros<T>(n,n);
  for (size_t i=0; i<n; ++i)
    A(i,i)=T(1);
  return A;
}




template<typename T>
M_Matrix<T>  Rand(const M_Matrix<T>& x)
{
  std::normal_distribution<> normal;
  std::random_device rd;
  std::mt19937_64 sto(rd());
  auto out=zeros<T>(x);
  for (std::size_t i=0; i<out.size(); ++i)
    out[i]=normal(sto);
  return out;
}

template<typename T>
M_Matrix<T>  Rand(const M_Matrix<T>& x, std::mt19937_64& sto)
{
  std::normal_distribution<> normal;
  auto out=zeros<T>(x);
  for (std::size_t i=0; i<out.size(); ++i)
    out[i]=normal(sto);
  return out;
}


/**
  Blas

SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*/


/**
     Transpose the first and multiply by the second
     @post transpMult(x,y)==Transpose(x)*y
     @remarks It is faster, since we save copying matrices
    */
template<typename T>
M_Matrix<T> TranspMult(const M_Matrix<T>& x,const M_Matrix<T>& y)
{
  // First it has to find out if the last dimension of x matches the first of y


  // now we build the M_Matrix result
  M_Matrix<T> z=zeros<T> (x.ncols(),y.ncols());


  /***  as fortran uses the reverse order for matrices and we want to
              avoid a copying operation, we calculate
                  Transpose(Z)=Transpose(y)*Transpose(x)

                  Transpose(matrix)=just plain matrix in C++ format


              */
  char  	TRANSA='N';
  char  	TRANSB='T';
  int  	M=y.ncols();
  int  	N=x.ncols();
  int  	K=x.nrows();
  double  ALPHA=1.0;
  double*  A=const_cast<double*> (&y[0]);
  int  	LDA=M;
  double*  B=const_cast<double*> (&x[0]);
  int  	LDB=N;
  double BETA=0.0;
  double * C=&z[0];
  int  	LDC=M;



  try{
    dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
  }
  catch (...)
  {
    std::cerr<<" dgemm_error";
  }
  return z;
}

/**
     Multiply by the Transpose of the second matrix
     @post MultTransp(x,y)==x*Transpose(y)
     @remarks It is faster, since we save copying matrices
    */
template<typename T>
M_Matrix<T> multTransp(const M_Matrix<T>& x,const M_Matrix<T>& y)
{
  // First it has to find out if the last dimension of x matches the first of y
  //ASSERT_NE(x.size(),0);
  //ASSERT_EQ(x.ncols(),ncols(y));
  // now we build the M_Matrix result
  M_Matrix<T> z=zeros<T> (x.nrows(),y.nrows());
  char  	TRANSA='T';
  char  	TRANSB='N';
  int  	M=y.nrows();
  int  	N=x.nrows();
  int  	K=x.ncols();
  double  ALPHA=1.0;
  double*  A=const_cast<double*> (&y[0]);
  int  	LDA=K;
  double*  B=const_cast<double*> (&x[0]);
  int  	LDB=K;
  double BETA=0.0;
  double * C=&z[0];
  int  	LDC=M;

  try
  {
    dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
  }
  catch(...)
  {
    std::cerr<<" dgemm_ error";
  }

  return z;
}



/**
    Matrix multiplication.
    @pre \p x.ncols()==y.nrows()
    @returns z=x*y
    @post z.nrows()=rows(x), z.ncols()=ncols(y)
    @post z(i,j)= sum on k of x(i,k)*y(k,j)
    @post assert(x.ncols()==y.nrows())
    */
template<typename T>
M_Matrix<T> operator*(const M_Matrix<T>& x,const M_Matrix<T>& y)
{
  // First it has to find out if the last dimension of x matches the
  //first of y
  // now we build the M_Matrix result
  M_Matrix<T> z=zeros<T> (x.nrows(),y.ncols());

  /***  as fortran uses the reverse order for matrices and we want to
          avoid a copying operation, we calculate
              Transpose(Z)=Transpose(y)*Transpose(x)

              Transpose(matrix)=just plain matrix in C++ format


          */
  char  	TRANSA='N';
  char  	TRANSB='N';
  int  	M=y.ncols();
  int  	N=x.nrows();
  int  	K=x.ncols();
  double  ALPHA=1.0;
  double*  A=const_cast<double*> (&y[0]);
  int  	LDA=M;
  double*  B=const_cast<double*> (&x[0]);
  int  	LDB=K;
  double BETA=0.0;
  double * C=&z[0];
  int  	LDC=M;


  try
  {
    dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
  }
  catch (...)
  {
    assert(false);
  }
  return z;
}



inline
M_Matrix<double> operator*(const M_Matrix<double>& x,
                           const M_Matrix<std::size_t>& y)
{
  // First it has to find out if the last dimension of x matches the
  //first of y
  // now we build the M_Matrix result
  M_Matrix<double> z(x.nrows(),y.ncols(),0.0);
  // we build the dimensions std::vector of the result
  for (size_t i=0; i<z.nrows(); i++)
    for (size_t j=0; j<z.ncols(); j++)
      for (size_t k=0; k<x.ncols(); k++)
        z(i,j)+=x(i,k)*y(k,j);
  return z;
}
template<typename T>
M_Matrix<T> operator-(const M_Matrix<T>& x)
{
  M_Matrix<T> out(x.nrows(),x.ncols());
  for (std::size_t i=0; i<out.size(); ++i)
    out[i]=-x[i];
  return out;
}

template<typename T>
M_Matrix<T> operator-(M_Matrix<T>&& x)
{
  for (std::size_t i=0; i<x.size(); ++i)
    x[i]=-x[i];
  return x;
}





/** @name Aritmetic Assigment Operations between a Matrix and a scalar
    (single element)
      */
//@{
/**
     Scalar Adition assignment.
     @returns a reference to itself
     @post all the values of the matrix are summed up by the value x
     */

template<typename T>
M_Matrix<T>& operator+=(M_Matrix<T>& itself, T x)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]+=x;
  return itself;
}


/**
     Scalar Subtraction assignment.
     @returns a reference to itself
     @post all the values of the matrix are sustracted by the value x
     */
template<typename T>
M_Matrix<T>& operator-=(M_Matrix<T>& itself, T x)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]-=x;
  return itself;
}





/**
     Scalar Multiplication assignment.
     @returns a reference to itself
     @post all the values of the matrix are multiplied by the value x
     */
template<typename T>
M_Matrix<T>& operator*=(M_Matrix<T>& itself, T x)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]*=x;
  return itself;
}


/**
     Scalar Division assignment.
     @returns a reference to itself
     @post all the values of the matrix are divided by the value x
     */
template<typename T>
M_Matrix<T>& operator/=(M_Matrix<T>& itself, T x)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]/=x;
  return itself;
}

//@}
/** @name Aritmetic Assigment Operations between two Matrices
      */
//@{






/*!
     Matrix Addition assignment.
     @returns a reference to itself
     @pre  same number of rows and columns
     @post all the values of the matrix are summed up by the corresponing
           values of the other matrix\n\n
           assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))

     */
template<typename T>
M_Matrix<T>& operator+=(M_Matrix<T>& itself,
                        const M_Matrix<T>& other)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]+=other[i];
  return itself;
}

/*!
      Matrix Subtraction assignment.
      @returns a reference to itself
      @pre  same number of rows and columns
      @post all the values of the matrix are sustracted by the corresponing
            values of the other matrix\n\n
            assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))
      */
template<typename T>
M_Matrix<T>& operator-=(M_Matrix<T>& itself,
                        const M_Matrix<T>& other)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]-=other[i];
  return itself;
}

/*!
      Matrix Addition assignment with typecast.
      @warning will not compile unless typecast T(S) is defined
      @returns a reference to itself
      @pre  same number of rows and columns
      @post all the values of the matrix are summed up by the corresponing
            values of the other matrix\n\n
            assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))

      */
template<typename T,typename S>
M_Matrix<T>& operator+=(M_Matrix<T>& itself,
                        const M_Matrix<S>& other)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]+=T(other[i]);
  return itself;
}

/*!
       Matrix Subtraction assignment with typecast.
      @warning will not compile unless typecast T(S) is defined.
       @returns a reference to itself
       @pre  same number of rows and columns
       @post all the values of the matrix are sustracted by the corresponing
             values of the other matrix\n\n
             assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))
       */
template<typename T,typename S>
M_Matrix<T>& operator-=(M_Matrix<T>& itself,
                        const M_Matrix<S>& other)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]-=T(other[i]);
  return itself;
}




//@}









/** @name Aritmetic operation between a Matrix and a scalar (single element)
      */
//@{

/**
     Scalar Addition.
     @returns a copy of the matrix with its values summed by x
     */
template<typename T>
M_Matrix<T> operator+(const M_Matrix<T>& x,T t)
{    // we build the M_Matrix result
  M_Matrix<T> z(x);
  z+=t;
  return z;
}

/**
     Scalar Addition reverse order.
     */
template<typename T>
M_Matrix<T> operator+(T t,const M_Matrix<T>& x)
{
  return x+t;
}

/**
     Scalar Subtraction.
     @returns a copy of the matrix with its values substracted by x
     */
template<typename T>
M_Matrix<T> operator-(const M_Matrix<T>& x,T t)
{    // we build the M_Matrix result
  M_Matrix<T> z(x);
  z-=t;
  return z;
};

/**
     Scalar Subtraction reverse order.
     */
template<typename T>
M_Matrix<T> operator-(T t,const M_Matrix<T>& x)
{
  return x-t;
}

/**
     Scalar Multiplication.
     @returns a copy of the matrix with its values multiplied by the value x
     */
template<typename T>
M_Matrix<T> operator*(const M_Matrix<T>& x,T t)
{    // we build the M_Matrix result
  M_Matrix<T> z(x);
  z*=t;
  return z;
}

/**
     Scalar Multiplication reverse order.
     */
template<typename T>
M_Matrix<T> operator*(T t,const M_Matrix<T>& x)
{
  return x*t;
}


/**
     Scalar Division.
     @returns a copy of the matrix with its values divided by x
     @returns a matrix of real numbers
 */
template<typename T>
M_Matrix<double> operator/(const M_Matrix<T>& x,T t)
{    // we build the M_Matrix result
  M_Matrix<double> z(x);
  z/=double(t);
  return z;
};

/**
     Division by inhomogeneus types

     */

template<typename T,typename S>
M_Matrix<double> operator/(const M_Matrix<T>& x,S t)
{    // we build the M_Matrix result
  M_Matrix<double> z(x);
  z/=double(t);
  return z;
};









/**
     Scalar Division reverse order.
     */
template<typename T>
M_Matrix<double> operator/(T t,const M_Matrix<T>& x)
{
  M_Matrix<double> out(x.nrows(),x.ncols());
  for (std::size_t i=0;i<x.size();i++)
    out[i]=double(t)/double(x[i]);

  return out;
}


//@}
/**
 @name Aritmetic operations applied between two Matrices
  */
//@{

/**
 Matrix sum, element wise.
 @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
 @return z, where z.nrows()=rows(x), z.ncols()=ncols(y) and z(i,j)= sum
 on k of x(i,k)*y(k,j)
 @warning it \c assert the preconditions
 */
template<typename T>
M_Matrix<T> operator+(const M_Matrix<T>& x,const M_Matrix<T>& y)
{
  M_Matrix<T> z(x);
  for (size_t i=0; i<z.nrows(); i++)
    for (size_t j=0; j<z.ncols(); j++)
      z(i,j)+=y(i,j);
  return z;
}


/**
 Matrix sustraction, element wise.
 @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
 @return z, where z.nrows()=rows(x), z.ncols()=ncols(y) and
 z(i,j)= sum on k of x(i,k)*y(k,j)
 @warning it \c assert the preconditions
 */

template<typename T>
M_Matrix<T> operator-(const M_Matrix<T>& x,const M_Matrix<T>& y)
{
  if(x.size()!=y.size())
    assert(false);
  if (x.nrows()!=y.nrows())
    assert(false);
  M_Matrix<T> z(x.nrows(),x.ncols());
  for (size_t i=0; i<z.size(); i++)
    // for (size_t j=0; j<z.ncols(); j++)
    z[i]=x[i]-y[i];
  return z;
}




/**
 Multiplication of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
template<typename T>
M_Matrix<T> elemMult(const M_Matrix<T>& x,const M_Matrix<T>& y)
{
  assert(x.size()==y.size());
  assert(x.nrows()==y.nrows());
  M_Matrix<T> z(x);
  for (size_t i=0; i<z.nrows(); i++)
    for (size_t j=0; j<z.ncols(); j++)
      z(i,j)*=y(i,j);
  return z;
}

/**
 Division of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
template<typename T,typename S>
M_Matrix<double> elemDiv(const M_Matrix<T>& x,const M_Matrix<S>& y)
{
  M_Matrix<double> z(x);
  for (size_t i=0; i<z.nrows(); i++)
    for (size_t j=0; j<z.ncols(); j++)
      z(i,j)/=double(y(i,j));
  return z;
}

/**
 Safe Division of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
template<typename T,typename S>
M_Matrix<double> elemDivSafe(const M_Matrix<T>& x,const M_Matrix<S>& y)
{
  M_Matrix<double> z(x);
  for (size_t i=0; i<z.nrows(); i++)
    for (size_t j=0; j<z.ncols(); j++)
      if (y(i,j)!=0)
        z(i,j)/=double(y(i,j));
      else
        z(i,j)=0;
  return z;
}
// @}




//@}


template<typename T>
std::ostream& operator<<(std::ostream& os,const M_Matrix<T>& x)
{
  os<<"[";
  for (std::size_t i=0; i<x.nrows(); ++i)
    {
      for (std::size_t j=0; j<x.ncols(); ++j)
        os<<x(i,j)<<" ";
      os<<";";
    }
  os<<"]";
  return os;
}


template<typename T>
std::istream& operator>>(std::istream& is,M_Matrix<T>& x)
{
  std::vector<T> o;
  std::size_t nrows=0;
  char ch;
  while ((is>>ch)&&(ch!='[')){}
  while (ch!=']')
    {
      std::string s;
      while ((is.get(ch))&&((ch!=']')&&ch!=';'))
        {
          s.push_back(ch);
        }
      std::stringstream ss(s);
      T e;
      while (ss>>e) o.push_back(e);
      if (!o.empty()) ++nrows;
    }
  std::size_t ncols=o.size()/nrows;
  x=M_Matrix<T>(nrows,ncols,o);
  return is;

}

template<typename T>
T maxAbs(const M_Matrix<T>& x)
{
  T m=std::abs(x[0]);
  for (std::size_t i; i<x.size(); ++i)
    if (std::abs(x[i])>m) m=std::abs(x[i]);
  return m;
}



template<typename T>
T norm_inf(const M_Matrix<T>& x)
{
  T n(0);
  for (size_t i=0; i<x.nrows(); ++i)
    {
      T sum(0);
      for (size_t j=0; j<x.ncols(); ++j)
        if (x(i,j)>0)
          sum+=x(i,j);
        else
          sum-=x(i,j);

      n=std::max(n,sum);
    }
  return n;
}




template<typename T>
M_Matrix<T> TranspMult(const M_Matrix<T>& x,const M_Matrix<T>& y);


template<typename T>
M_Matrix<T> multTransp(const M_Matrix<T>& x,const M_Matrix<T>& y);

template<typename T>
M_Matrix<T> operator*(const M_Matrix<T>& x,const M_Matrix<T>& y);

inline
M_Matrix<double> operator*(const M_Matrix<std::size_t>& x,
                           const M_Matrix<double>& y)
{
  // First it has to find out if the last dimension of x matches the
  //first of y
  //ASSERT_EQ(x.ncols(),y.nrows());
  // now we build the M_Matrix result
  M_Matrix<double> z(x.ncols(),y.ncols(),0.0);
  // we build the dimensions std::vector of the result
  for (size_t i=0; i<z.ncols(); i++)
    for (size_t j=0; j<z.ncols(); j++)
      for (size_t k=0; k<x.ncols(); k++)
        z(i,j)+=x(i,k)*y(k,j);
  return z;
}



M_Matrix<double> operator*(const M_Matrix<double>& x,
                           const M_Matrix<std::size_t>& y);


inline
std::vector<double> operator*(const std::vector<double>& x,
                              const M_Matrix<double>& y)
{
  std::vector<double> out(y.ncols(),0);
  for (std::size_t i=0; i<x.size(); ++i)
    for (std::size_t j=0; j<y.ncols(); ++j)
      out[j]+=x[i]*y(i,j);
  return out;
}



inline
double xTSigmaX(const M_Matrix<double> &vector, const M_Matrix<double> &matrix)
{
  double sum=0;
  for (std::size_t i=0; i<matrix.nrows(); ++i)
    {
      sum+=vector[i]*matrix(i,i)*vector[i];
      for (std::size_t j=i+1; j<matrix.ncols();++j)
        sum+=2*vector[i]*matrix(i,j)*vector[j];
    }
  return sum;
}
inline
double xTSigmaX(const std::vector<double> &v, const M_Matrix<double> &matrix)
{
  double sum=0;
  for (std::size_t i=0; i<matrix.nrows(); ++i)
    {
      sum+=v[i]*matrix(i,i)*v[i];
      for (std::size_t j=i+1; j<matrix.ncols();++j)
        sum+=2*v[i]*matrix(i,j)*v[j];
    }
  return sum;
}






/**
        Exception class for matrix singularity (i.e, that do not have Matrix
      Inverse)
      */
class SingularMatrix_error: public std::runtime_error
{
public:
  SingularMatrix_error(const char* msg):std::runtime_error(msg){}
};



template<typename T>
M_Matrix<T> inv(const M_Matrix<T>& a)

{
  if ((a.size()>0)&&(a.nrows()==a.ncols()))
    {
      double *A;
      int info=0;
      //  char msg[101];
      int *ipiv;
      int lwork;
      int n =a.ncols();
      int m=n;
      M_Matrix<T> B(a);
      int dla=n;
      //A=new double[n*n];
      A= new double[n*n]; //more efficient code
      for (size_t k = 0; k < size_t(n*n); k++)
        *(A+k) = a[k];

      ipiv = new int[n];

      dgetrf_(&n, &m, A, &dla,ipiv,&info);

      lwork= n*n;
      double *work = new double[n*n];

      dgetri_(&n,A,&dla,ipiv,work,&lwork,&info);

      for (size_t k = 0; k < size_t(n*n); k++)
        B[k] = *(A+k);
      delete [] A;
      delete [] ipiv;
      delete [] work;
      if (info!=0)
        {
          throw SingularMatrix_error("cannot invert a singular matrix");
        }
      return B;
    }
  else
    return a;
}



inline
M_Matrix<double> invSafe(const M_Matrix<double>& matrix)

{
  M_Matrix<double> inverse;
  try
  {
    inverse=inv(matrix);
  }
  catch (SingularMatrix_error)
  {
    inverse={};
  }
  return inverse;
}


inline
bool isnan(const M_Matrix<double>& x)
{
  for (std::size_t i=0; i<x.size(); ++i)
    if (std::isnan(x[i]))
      return true;
  return false;
}





/**
      Transpose
      @post (Transpose(x))(i,j)==x(j,i)
      @returns the Transpose of the Matrix
      */

template<class T>
M_Matrix<T>  Transpose(const M_Matrix<T>& x)
{
  M_Matrix<T> tr(x.ncols(),x.nrows());
  for (size_t i=0; i<tr.nrows(); i++)
    for (size_t j=0; j<tr.ncols(); j++)
      tr(i,j)=x(j,i);
  return tr;
}



/**
       Diagonal of Matrix or Diagonal Matrix
       It has two behaviors:
       - If the input is a single column or a single row, it builds a diagonal
       Matrix with it
       - If the input is a Matrix, it returns the values of its diagonal

      */
template<typename T>
M_Matrix<T> diag(const M_Matrix<T>& x)
{
  size_t nr=x.nrows();
  size_t nc=x.ncols();
  if ((nr>1)&(nc>1))
    {
      std::size_t n=std::min(nr,nc);
      M_Matrix<T> diagM= zeros<T> ( 1,n);
      for (size_t i=0; i<n; ++i)
        diagM(0,i)=x(i,i);
      return diagM;
    }
  else
    {
      nr=std::max(nr,nc);
      M_Matrix<T> diagM=zeros<T>(nr,nr);
      for (size_t i=0; i<nr; ++i)
        diagM(i,i)=x[i];
      return diagM;
    }

}



template<typename T>
M_Matrix<T> diag_landa(const M_Matrix<T>& x,double landa)
{
  double landa1=landa+1;
  M_Matrix<T> diagM(x);
  for (size_t i=0; i<x.nrows(); ++i)
    diagM(i,i)*=landa1;
  return diagM;

}




/**
       Product of the Diagonal of a Matrix

      */
template<typename T>
T diagProduct(const M_Matrix<T>& x)
{
  size_t nr=x.nrows();
  size_t nc=x.ncols();
  double diagprod=1;
  std::size_t n=std::min(nr,nc);
  for (size_t i=0; i<n; ++i)
    diagprod*=x(i,i);
  return diagprod;


}




template<typename T>
M_Matrix<T> col_vector(const M_Matrix<T>& x)
{
  M_Matrix<T> colvec(x.size(),1);
  for (std::size_t i=0; i<x.size(); ++i)
    colvec[i]=x[i];
  return colvec;
}
template<typename T>
M_Matrix<T> row_vector(const M_Matrix<T>& x)
{
  M_Matrix<T> rowvec(1,x.size());
  for (std::size_t i=0; i<x.size(); ++i)
    rowvec[i]=x[i];
  return rowvec;
}








/**

   SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
  *
  *  -- LAPACK routine (version 3.3.1) --
  *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  *  -- April 2011                                                      --
  *
  *     .. Scalar Arguments ..
        CHARACTER          UPLO
        INTEGER            INFO, LDA, N
  *     ..
  *     .. Array Arguments ..
        DOUBLE PRECISION   A( LDA, * )
  *     ..
  *
  *  Purpose
  *  =======
  *
  *  DPOTRF computes the Cholesky factorization of a real symmetric
  *  positive definite matrix A.
  *
  *  The factorization has the form
  *     A = U**T * U,  if UPLO = 'U', or
  *     A = L  * L**T,  if UPLO = 'L',
  *  where U is an upper triangular matrix and L is lower triangular.
  *
  *  This is the block version of the algorithm, calling Level 3 BLAS.
  *
  *  Arguments
  *  =========
  *
  *  UPLO    (input) CHARACTER*1
  *          = 'U':  Upper triangle of A is stored;
  *          = 'L':  Lower triangle of A is stored.
  *
  *  N       (input) INTEGER
  *          The order of the matrix A.  N >= 0.
  *
  *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  *          N-by-N upper triangular part of A contains the upper
  *          triangular part of the matrix A, and the strictly lower
  *          triangular part of A is not referenced.  If UPLO = 'L', the
  *          leading N-by-N lower triangular part of A contains the lower
  *          triangular part of the matrix A, and the strictly upper
  *          triangular part of A is not referenced.
  *
  *          On exit, if INFO = 0, the factor U or L from the Cholesky
  *          factorization A = U**T*U or A = L*L**T.
  *
  *  LDA     (input) INTEGER
  *          The leading dimension of the array A.  LDA >= max(1,N).
  *
  *  INFO    (output) INTEGER
  *          = 0:  successful exit
  *          < 0:  if INFO = -i, the i-th argument had an illegal value
  *          > 0:  if INFO = i, the leading minor of order i is not
  *                positive definite, and the factorization could not be
  *                completed.
  *
  *  =====================================================================


    */
namespace{

  extern "C" void dpotrf_(char * 	UPLO,
                          int * N,
                          double * A,
                          int * LDA,
                          int * INFO);

  M_Matrix<double> UT(const M_Matrix<double>& x)
  {
    M_Matrix<double> y(x.nrows(),x.ncols());
    for (std::size_t i=0;i<x.nrows();i++)
      {
        for (std::size_t j=0;j<i; j++)
          y(j,i)=0;
        for (std::size_t j=i;j<x.ncols(); j++)
          y(j,i)=x(i,j);
      }
    return y;
  }

  M_Matrix<double> LT(const M_Matrix<double>& x)
  {
    M_Matrix<double> y(x.nrows(),x.ncols());
    for (std::size_t i=0;i<x.nrows();i++)
      {
        for (std::size_t j=0;j<i+1; j++)
          y(j,i)=x(i,j);
        for (std::size_t j=i+1;j<x.ncols(); j++)
          y(j,i)=0;
      }
    return y;
  }



}

inline
M_Matrix<double> chol(const M_Matrix<double>& x,const std::string& kind)
{

  if ((x.nrows()!=x.ncols())||x.size()==0)
    return M_Matrix<double>();
  char UPLO='L';
  M_Matrix<double> res;
  if (kind!="lower")
    {
      UPLO='U';
      res=UT(x);
    }
  else
    {
      res=LT(x);
    }
  int N=x.nrows();
  int LDA=N;
  int INFO;
  double* A=new double[x.size()];
  for (std::size_t i=0; i<x.size(); i++)
    A[i]=res[i];

  if (LDA==0)
    return M_Matrix<double>();
  try
  {
    dpotrf_(&UPLO,&N,A,&LDA,&INFO);
  }
  catch (...)
  {
    std::cerr<<" error";
  };
  for (std::size_t i=0; i<x.size(); i++)
    res[i]=A[i];
  res=Transpose(res);
  delete [] A;
  return res;

}




#endif // MATRIX_H
