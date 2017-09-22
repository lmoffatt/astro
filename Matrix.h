#ifndef MATRIX_H
#define MATRIX_H
/*!
 * @file Matrix.h


 */


#include <limits> // for std::numeric_limits
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
#include <algorithm>
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


template<typename T>
std::size_t i_max(const std::vector<T>& x)
{
  std::size_t j=0;
  for (std::size_t i=1; i<x.size(); ++i)
    if (x[i]>x[j]) j=i;
  return j;
}









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

  extern "C" void ssymm_ 	(char*  	SIDE,
                                 char*  	UPLO,
                                 int*  	M,
                                 int*  	N,
                                 double*  	ALPHA,
                                 double * /*, dimension(lda,*)*/  	A,
                                 int*  	LDA,
                                 double */*, dimension(ldb,*)*/  	B,
                                 int*  	LDB,
                                 double *  	BETA,
                                 double * /*, dimension(ldc,*) */ 	C,
                                 int*  	LDC
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
    const M_Matrix<T>& m;
    std::size_t i_;
    std::size_t j_;
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


  M_Matrix()=default;

  M_Matrix (const M_Matrix<T> & sample)=default;
  M_Matrix (M_Matrix<T> && sample)=default;

  M_Matrix& operator= (const M_Matrix<T> & sample)=default;
  M_Matrix& operator= (M_Matrix<T> && sample)=default;

  M_Matrix (std::size_t nrows,std::size_t ncols)
    : type_(FULL),
      _nrows(nrows),
      _ncols(ncols),
      _data(nrows*ncols)
  {
  }


  template<typename S>
  M_Matrix (const M_Matrix<S> & sample)
    : type_(TYPE(sample.type())),
      _nrows(sample.nrows()),
      _ncols(sample.ncols()),
      _data(sample.size())

  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      (*this)[i] = sample[i];

  }
  M_Matrix (const std::vector<std::vector<T>> & sample)
    : type_(FULL),
      _nrows(sample.size()),
      _ncols(sample[0].size()),
      _data(sample.size()*sample[0].size())
  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      for (std::size_t j=0; j<sample[0].size(); ++j)
        (*this)(i,j) = sample[i][j];

  }


  template<typename S>
  M_Matrix (const std::vector<std::vector<S>> & sample)
    :type_(FULL),
      _nrows(sample.size()),
      _ncols(sample[0].size()),
      _data(sample.size()*sample[0].size())
  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      (*this)[i] = sample[i];

  }


  ~M_Matrix()
  {
    //  if (_data>0)
    //	delete [] _data;
  }



  M_Matrix(std::size_t nrows_,std::size_t ncols_, std::vector<T> data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(data)
  {

  }

  M_Matrix(std::size_t nrows_,std::size_t ncols_, const M_Matrix<T>& data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(nrows_*ncols_)
  {
    for (size_t i=0; i<_data.size(); ++i)
      _data[i]=data[i];
  }

  template <typename S>  M_Matrix(std::size_t nrows_,std::size_t ncols_, const M_Matrix<S>& data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(nrows_*ncols_)
  //   _data(new T[_ncells])
  {

    for (size_t i=0; i<_data.size(); ++i)
      _data[i]=data[i];
  }


  M_Matrix(std::size_t nrows_,std::size_t ncols_,T data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    //   _data(new T[_ncells])
    _data(nrows_*ncols_,data)
  {

  }


  template<class F>
  M_Matrix<T>
  apply(const F& f)const
  {
    M_Matrix<T> out(nrows(),ncols());
    for (std::size_t i=0; i<size(); ++i)
      out[i]=f((*this)[i]);
    return out;
  }

  M_Matrix<T>& operator=(T X)
  {
    for (std::size_t i=0; i<size(); ++i) _data[i]=X;
    return *this;
  }



  size_t size()const
  {
    return _data.size();
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
    assert(i<nrows());
    assert(j<ncols());
    switch(type_){
      case FULL:
        return (*this)[i*_ncols+j];
      case SYMMETRIC:
        if(i>=j)
          return (*this)[i*(i+1)/2+j]; //stores at lower triangular portion
        else
          return (*this)[j*(j+1)/2+i];
      case DIAGONAL:
        assert(i==j);
        return (*this)[i];
      case SCALAR_DIAGONAL:
        assert(i==j);
        return (*this)[0];
      case SCALAR_FULL:
        return (*this)[0];
      case ZERO:
        assert(false);
        return (*this)[0];
      default:
        assert(false);
        return (*this)[0];
      }
  }

  T const&  operator() (std::size_t i,std::size_t j) const
  {
    assert(i<nrows());
    assert(j<ncols());
    switch(type_){
      case FULL:
        return (*this)[i*_ncols+j];
      case SYMMETRIC:
        if(i>=j)
          return (*this)[i*(i+1)/2+j]; //stores at lower triangular portion
        else
          return (*this)[j*(j+1)/2+i];
      case DIAGONAL:
        if(i==j)
          return (*this)[i];
        else
          return zero_;
      case SCALAR_DIAGONAL:
        if(i==j)
          return (*this)[0];
        else
          return zero_;
      case SCALAR_FULL:
        return (*this)[0];
      case ZERO:
        return zero_;
      default:
        return zero_;
      }
  }




  /** @name  Accesing all the values of a Row or Column at once
     */
  //@{

  /**
    Replacement of the ith Row.
    @param iRow the extracted row (0 is the first as always in C)
    @param newValues contains the numbers used to relace iRow
    @param dummy an ignored variable to differentiate from column
    extraction
    @pre iRow is smaller than nrows(*this)
    @pre size(newValues)==ncols(*this)
    @returns *this
    @post if the preconditions are met, the ith row is replaced by
    newValues
    @post newValues is treated as a vector, although is a Matrix. Its
    internal structure (i.e., ncols and nrows) is ignored.
    @post assert the precoditions
    */

  M_Matrix<T>&  operator() (std::size_t iRow,
                            std::string /*dummy*/,
                            const M_Matrix<T>& newValues)
  {
    assert(type_==FULL);
    assert(iRow<nrows());//number of rows
    assert(newValues.size()==ncols()); //number of columns
    for (std::size_t j=0; j<std::min(ncols(),size()); j++)
      this->operator()(iRow,j)=newValues[j];
    return *this;
  }





  /**
    Replacement of the jth Column.
    @param newValues contains the numbers used to relace jth Column
    @pre newValues is treated as a vector, although is a Matrix. Its
    internal structure (i.e., ncols and nrows) is ignored.
    @param jColumn the replaced column (0 is the first as always in C)
    @param dummy an ignored variable to differentiate from column
    extraction
    @pre jColumn is smaller than ncols(*this)
    @pre size(newValues)==nrows(*this)
    \returns *this
    @post if the preconditions are met, the jth Column is replaced by
    newValues
    @post assert the precoditions
    */

  template <typename V>
  M_Matrix<T>&  operator() (const std::string /*dummy*/,
                            std::size_t jColumn,
                            const V& newValues)
  {
    assert(type_==FULL);

    for (std::size_t i=0; i<std::min(nrows(),newValues.size()); i++)
      this->operator()(i,jColumn)=newValues[i];
    //  assert(ndim>1);
    //  assert(i<n[0]);//number of rows
    //  assert(j<n[1]); //number of columns
    return *this;
  }





  /**
    Copy of the ith Row
    @pre iRow is smaller than nrows(*this)
    @param iRow the extracted row (0 is the first as always in C)
    @param dummy an ignored variable to differentiate from column
    extraction
    \returns a 1-row ncols Matrix with the values of the ith row
    */

  M_Matrix<T>  operator() (std::size_t iRow,
                           const std::string /*dummy*/
                           ) const
  {
    assert(type_==FULL);
    M_Matrix<T> out(1,ncols());
    for (std::size_t j=0; j<ncols(); j++)
      out[j]=this->operator()(iRow,j);
    return out;
  }






  /**
    Copy of the jth Column
    @pre jColumn is smaller than ncols(*this)
    @param jColumn the extracted column (0 is the first as always in C)
    @param dummy is an ignored const string (like "") to differentiate
    from row extraction
    \returns a nrows 1-column Matrix with the values of the jth column
    */

  M_Matrix<T>  operator() (std::string /*dummy*/,
                           std::size_t jColumn
                           ) const
  {
    assert(type_==FULL);
    M_Matrix<T> out(nrows(),1);
    for (std::size_t i=0; i<nrows(); i++)
      out[i]=(*this)(i,jColumn);
    return out;
  }


  void clear()
  {
    type_=ZERO;
    _nrows=0;
    _ncols=0;
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

  enum TYPE{ZERO,FULL,SYMMETRIC,DIAGONAL,SCALAR_FULL,SCALAR_DIAGONAL};
  /**
     Returns a custom sized Matrix filled with ones
    @post (ones(n,m))(i,j)==T(1)
     */
  M_Matrix(std::size_t nrows_,std::size_t ncols_,TYPE t,T data):
    type_(t),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(getSize(t,nrows_,ncols_),data)
  {
  }
  M_Matrix(std::size_t nrows_,TYPE t,T data):
    type_(t),
    _nrows(nrows_),
    _ncols(nrows_),
    _data(getSize(t,nrows_),data)
  {
  }

  M_Matrix(std::size_t nrows_,std::size_t ncols_,TYPE t):
    type_(t),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(getSize(t,nrows_,ncols_))
  {
  }
  M_Matrix(std::size_t nrows_,TYPE t):
    type_(t),
    _nrows(nrows_),
    _ncols(nrows_),
    _data(getSize(t,nrows_))
  {
  }

  static
  M_Matrix
  unpackForLapack(const M_Matrix<double>& packedMatrix, char UPLO='L')
  {
    switch (packedMatrix.type())
      {
      case SYMMETRIC:
        {
          M_Matrix<double> out(packedMatrix.nrows(),packedMatrix.ncols());
          if (UPLO=='L')
            {
              for (std::size_t i=0; i<out.nrows(); ++i)
                for (std::size_t j=0; j<=i; ++j)
                  out(i,j)=packedMatrix(i,j);
              return out;
            }
          else
            {
              for (std::size_t i=0; i<out.nrows(); ++i)
                for (std::size_t j=i; j<=out.ncols(); ++j)
                  out(i,j)=packedMatrix(i,j);
              return out;
            }

        }
      default:
        return M_Matrix<double>(packedMatrix);
      }
  }


  M_Matrix
  full()const
  {
    M_Matrix out(nrows(),ncols());
    for (std::size_t i=0; i<nrows(); i++)
      for (std::size_t j=0; j<ncols(); ++j)
        out(i,j)=(*this)(i,j);
    return out;
  }

  static
  std::size_t getSize(TYPE t,std::size_t nrows,std::size_t ncols)
  {
    switch (t)
      {
      case ZERO: return 0;
      case FULL: return nrows*ncols;
      case SYMMETRIC:
        {
          assert(nrows==ncols);
          return (nrows*(ncols+1))/2;
        }
      case DIAGONAL:
        {
          return std::min(nrows,ncols);
          break;
        }
      case SCALAR_FULL:
      case SCALAR_DIAGONAL:
      default:
        return 1;
      }
  }
  std::size_t getSize(TYPE t,std::size_t nrows)
  {
    switch (t)
      {
      case ZERO: return 0;
      case FULL: return nrows*nrows;
      case SYMMETRIC:
        {
          return (nrows*(nrows+1))/2;
        }
      case DIAGONAL:
        {
          return nrows;
          break;
        }
      case SCALAR_FULL:
      case SCALAR_DIAGONAL:
      default:
        return 1;
      }
  }


  TYPE type()const {return type_;}


  bool isSymmetric()const {
    switch (type_)
      {
      case ZERO:
      case SYMMETRIC:
        return true;
      case SCALAR_DIAGONAL:
      case DIAGONAL:
      case SCALAR_FULL:
        return ncols()==nrows();
      case FULL:
        return false;
      }
  }

  /**
    Matrix of zeros with the shape of the provided Matrix
     @post nrows(ones(x))==x.nrows()   ncols(ones(x))=x.ncols()

    */

  template<typename S>
  friend
  auto
  TranspMult(const M_Matrix<T>& x,const M_Matrix<S>& y)->
  M_Matrix<decltype(operator*(std::declval<T>(),std::declval<S>()))>;
  M_Matrix(TYPE t, std::size_t nrows, std::size_t ncols, std::vector<T> data):
    type_(t),_nrows(nrows),_ncols(ncols),_data(data){}

  template<typename...Ts , template<typename...>class V>
  M_Matrix(TYPE t, std::size_t nrows, std::size_t ncols, V<Ts...> data):
    type_(t),_nrows(nrows),_ncols(ncols),_data(data.size()){
    for (std::size_t i=0; i<data.size(); ++i)
      _data[i]=data[i];
  }

  /*!
       Matrix Addition assignment.
       @returns a reference to itself
       @pre  same number of rows and columns
       @post all the values of the matrix are summed up by the corresponing
             values of the other matrix\n\n
             assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))

       */

  template<typename S>
  M_Matrix<T>& operator+=( const M_Matrix<S>& other)
  {
    return Aditive_Assigment_Operator([](T& it,const S& two){it+=two; return it;},
    *this,other);
  }
  M_Matrix& operator+=( const M_Matrix& other)
  {
    return Aditive_Assigment_Operator([](T& it,const T& two){it+=two; return it;},
    *this,other);
  }

  /*!
        Matrix Subtraction assignment.
        @returns a reference to itself
        @pre  same number of rows and columns
        @post all the values of the matrix are sustracted by the corresponing
              values of the other matrix\n\n
              assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))
        */
  template<typename S>
  M_Matrix<T>& operator-=( const M_Matrix<S>& other)
  {
    return Aditive_Assigment_Operator([](T& it,const S& two){it-=two; return it;},
    *this,other);
  }


private:
  TYPE type_;
  std::size_t          _nrows;    /**< number of rows */
  std::size_t          _ncols;    /**< number of columns */

  /** internal data
          @remarks can be a pointer or a vector
          @remarks now is a vector for debugging purposes
          */
  //	T                      *_data; /**< pointer to the data  */
  std::vector<T>        _data;
  T zero_=T();
};

template<typename T>

std::size_t nrows(const M_Matrix<T>& x){return x.nrows();}

template<typename T>
std::size_t ncols(const M_Matrix<T>& x){return x.ncols();}

template<typename T>
std::size_t size(const M_Matrix<T>& x){return x.size();}


template<typename T>
M_Matrix<T>  ones(size_t nrows, size_t ncols)
{
  return M_Matrix<T>(nrows,ncols,M_Matrix<T>::SCALAR_FULL,T(1));
}

template<typename T>
M_Matrix<T>  zeros(size_t nrows, size_t ncols)
{
  return M_Matrix<T>(nrows,ncols,M_Matrix<T>::ZERO);
}

/**
  Identity Matrix of the specified size
 */
template<typename T>
M_Matrix<T> eye(std::size_t n)
{
  return M_Matrix<T>(n,n,M_Matrix<T>::SCALAR_DIAGONAL,T(1));
}


template<typename T, class Predicate>
bool all(const M_Matrix<T>& x, const Predicate& p)
{
  for (std::size_t i=0; i< x.size(); ++i)
    if (!p(x[i]))
      return false;
  return true;
}

template<typename T, class Predicate>
bool any(const M_Matrix<T>& x, const Predicate& p)
{
  for (std::size_t i=0; i< x.size(); ++i)
    if (p(x[i]))
      return true;
  return false;
}

template<typename T>
M_Matrix<T>  Rand(const M_Matrix<T>& x)
{
  std::normal_distribution<> normal;
  std::random_device rd;
  std::mt19937_64 sto(rd());
  auto out=M_Matrix<T>(x);
  for (std::size_t i=0; i<out.size(); ++i)
    out[i]=normal(sto);
  return out;
}

template<typename T>
M_Matrix<T>  Rand(const M_Matrix<T>& x, std::mt19937_64& sto)
{
  std::normal_distribution<> normal;
  auto out=M_Matrix<T>(x);
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
      Transpose
      @post (Transpose(x))(i,j)==x(j,i)
      @returns the Transpose of the Matrix
      */

template<class T>
M_Matrix<T>  Transpose(const M_Matrix<T>& x)
{
  switch(x.type())
    {
    case M_Matrix<T>::FULL:
      {
        M_Matrix<T> out(x.ncols(),x.nrows());
        for (std::size_t i=0; i<x.nrows(); i++)
          for (std::size_t j=0; j<x.ncols(); ++j)
            out(j,i)=x(i,j);
        return out;
      }
    case M_Matrix<T>::ZERO:
    case M_Matrix<T>::SYMMETRIC:
    case M_Matrix<T>::DIAGONAL:
    case M_Matrix<T>::SCALAR_DIAGONAL:
    case M_Matrix<T>::SCALAR_FULL:
    default:
      {
        M_Matrix<T> out(x);
        return out;
      }
    }
}


/**
     Transpose the first and multiply by the second
     @post transpMult(x,y)==Transpose(x)*y
     @remarks It is faster, since we save copying matrices
    */
template<typename T, typename S>
auto
TranspMult(const M_Matrix<T>& x,const M_Matrix<S>& y)->
M_Matrix<decltype(operator*(std::declval<T>(),std::declval<S>()))>
{


  assert(nrows(x)==nrows(y));
  // First it has to find out if the last dimension of x matches the first of y
  // now we build the M_Matrix result
  typedef decltype(operator*(std::declval<T>(),std::declval<S>())) R;
  M_Matrix<R> z(x.ncols(),y.ncols(), R(0));
  for (std::size_t i=0; i<z.nrows(); ++i)
    for (std::size_t j=0; j<z.ncols(); ++j)
      for (std::size_t k=0; k<x.nrows(); ++k)
        z(i,j)+=x(k,i)*y(k,j);
  return z;
}

inline
M_Matrix<double> Lapack_Full_Product(const M_Matrix<double>& x,const M_Matrix<double>& y,
                                     bool transpose_x, bool transpose_y);

inline
M_Matrix<double> Lapack_Symmetric_Transpose_Product(const M_Matrix<double>& Sym,const M_Matrix<double>& Reg, bool SymRegT);



template<typename T, typename S>
auto
Matrix_Multiplication (const M_Matrix<T>& one, const M_Matrix<S>& other)
->M_Matrix<std::result_of<decltype(std::declval<T>(),std::declval<S>())>>
{
                                                                          assert(one.ncols()==other.nrows());
                                                                          typedef   std::result_of<decltype(std::declval<T>()*std::declval<S>())> R;
                                                                          if (one.type()==M_Matrix<T>::FULL)
{
                                                                          if (other.type()==M_Matrix<S>::FULL)
                                                                          return Full_Product(one,other);
                                                                          else if (other.type()==M_Matrix<S>::SYMMETRIC)
                                                                          return Full_Product(one,other);
                                                                          else if ((other.type()==M_Matrix<S>::DIAGONAL)
                                                                                   || (other.type()==M_Matrix<S>::SCALAR_DIAGONAL))
{
                                                                          M_Matrix<R> out(one.nrows(),other.ncols());
                                                                          for (std::size_t i=0; i<out.nrows(); ++i)
{
  for (std::size_t j=0; j<one.ncols(); ++j)
    out(i,j)=one(i,j)*other(j,j);
  for (std::size_t j=one.ncols(); j<out.ncols();++j)
    out(i,j)=0;
}
return out;
}
else if (other.type()==M_Matrix<S>::SCALAR_FULL)
{
  M_Matrix<R> out(one.nrows(),other.ncols());
  for (std::size_t i=0; i<out.nrows(); ++i)
    {
      R sum(0);
      for (std::size_t j=0; j<one.ncols(); ++j)
        sum+=one(i,j);
      sum*=other[0];
      for (std::size_t k=0; k<other.ncols(); ++k)
        out(i,k)=sum;
    }
  return out;
}
else // ZERO!
{
return M_Matrix<R>
(one.nrows(),other.ncols(),M_Matrix<R>::ZERO);

}

}
else   if (one.type()==M_Matrix<T>::SYMMETRIC)
{
  if (other.type()==M_Matrix<S>::FULL)
    {
      return Full_Product(one,other);
    }
  else if (other.type()==M_Matrix<S>::SYMMETRIC)
    {
      return Full_Product(one,other);
    }
  else if ((other.type()==M_Matrix<S>::DIAGONAL)
           || (other.type()==M_Matrix<S>::SCALAR_DIAGONAL))
    {
      M_Matrix<R> out
          (one.nrows(),other.ncols());
      for (std::size_t i=0; i<out.nrows(); ++i)
        {
          for (std::size_t j=0; j<=i; ++j)
            out(i,j)=one(i,j)*other(j,j);
        }
      return out;
    }
  else if (other.type()==M_Matrix<S>::SCALAR_FULL)
    {
      M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t i=0; i<out.nrows(); ++i)
        {
          R sum=one(i,i);
          for (std::size_t j=0; j<i; ++j)
            sum+=one(i,j);
          sum*=other[0];
          for (std::size_t k=0; k<other.ncols(); ++k)
            out(i,k)=sum;
        }
      return out;
    }
  else // ZERO!
    {
      return M_Matrix<R>
          (one.nrows(),other.ncols(),M_Matrix<R>::ZERO);

    }

}
else  if (one.type()==M_Matrix<T>::DIAGONAL)
{
  if (other.type()==M_Matrix<S>::FULL)
    {
      M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t i=0; i<other.nrows(); ++i)
        {
          for (std::size_t j=0; j<other.ncols(); ++j)
            out(i,j)=one(i,i)*other(i,j);
        }
      for (std::size_t i=one.ncols(); i<one.ncols();++i)
        for (std::size_t j=0; j<other.ncols(); ++j)
          out(i,j)=R(0);

      return out;
    }
  else if (other.type()==M_Matrix<S>::SYMMETRIC)
    {
      M_Matrix<R> out
          (one.nrows(),other.ncols());
      for (std::size_t i=0; i<other.nrows(); ++i)
        {
          for (std::size_t j=0; j<out.ncols(); ++j)
            out(i,j)=one(i,i)*other(i,j);
        }
      return out;
    }
  else if ((other.type()==M_Matrix<S>::DIAGONAL)
           || (other.type()==M_Matrix<S>::SCALAR_DIAGONAL))
    {
      M_Matrix<R> out
          (one.nrows(),other.ncols(),M_Matrix<R>::DIAGONAL);
      for (std::size_t i=0; i<out.nrows(); ++i)
        {
          out(i,i)=one(i,i)*other(i,i);
        }
      return out;
    }
  else if (other.type()==M_Matrix<S>::SCALAR_FULL)
    { M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t i=0; i<out.nrows(); ++i)
        {
          double d=one(i,i)*other(i,i);
          for (std::size_t j=0; j<other.ncols(); ++j)
            out(i,j)=d;
        }
      return out;
    }
  else // ZERO!
    {
      return M_Matrix<R>
          (one.nrows(),other.ncols(),M_Matrix<R>::ZERO);

    }

}

else  if (one.type()==M_Matrix<T>::SCALAR_DIAGONAL)
{
  if (other.type()==M_Matrix<S>::FULL)
    {
      M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t i=0; i<other.nrows(); ++i)
        {
          for (std::size_t j=0; j<other.ncols(); ++j)
            out(i,j)=one(i,i)*other(i,j);
        }
      for (std::size_t i=one.ncols(); i<one.ncols();++i)
        for (std::size_t j=0; j<other.ncols(); ++j)
          out(i,j)=0;

      return out;
    }
  else if (other.type()==M_Matrix<S>::SYMMETRIC)
    {
      M_Matrix<R> out
          (one.nrows(),other.ncols(),M_Matrix<R>::SYMMETRIC);
      for (std::size_t i=0; i<other.nrows(); ++i)
        {
          for (std::size_t j=0; j<=i; ++j)
            out(i,j)=one(i,i)*other(i,j);
        }
      return out;
    }
  else if (other.type()==M_Matrix<S>::DIAGONAL)

    {
      M_Matrix<R> out
          (one.nrows(),other.ncols(),M_Matrix<R>::DIAGONAL);
      for (std::size_t i=0; i<out.nrows(); ++i)
        {
          out(i,i)=one(i,i)*other(i,i);
        }
      return out;
    }
  else if ((other.type()==M_Matrix<S>::SCALAR_DIAGONAL))
    {
      return M_Matrix<R>
          (one.nrows(),other.ncols(),M_Matrix<R>::SCALAR_DIAGONAL,
           one[0]*other[0]);
    }
  else if (other.type()==M_Matrix<S>::SCALAR_FULL)
    {
      return  M_Matrix<R>
          (one.nrows(),other.ncols(),
           M_Matrix<R>::SCALAR_FULL,one[0]*other[0]);
    }
  else // ZERO!
    {
      return M_Matrix<R>
          (one.nrows(),other.ncols(),M_Matrix<R>::ZERO);

    }

}

else  if (one.type()==M_Matrix<T>::SCALAR_FULL)
{
  if (other.type()==M_Matrix<S>::FULL)
    {
      M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t k=0; k<other.ncols(); ++k)
        {
          R sum(0);
          for (std::size_t j=0; j<one.ncols(); ++j)
            sum+=other(j,k);
          sum*=one[0];
          for (std::size_t i=0; i<out.nrows(); ++i)
            out(i,k)=sum;
        }
      return out;
    }
  else if (other.type()==M_Matrix<S>::SYMMETRIC)
    {
      M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t i=0; i<other.nrows(); ++i)
        {
          R sum=other(i,i);
          for (std::size_t j=0; j<i; ++j)
            sum+=other(i,j);
          sum*=one[0];
          for (std::size_t k=0; k<one.nrows(); ++k)
            out(k,i)=sum;
        }
      return out;
    }
  else if (other.type()==M_Matrix<S>::DIAGONAL)
    {
      M_Matrix<R> out(one.nrows(),other.ncols());
      for (std::size_t i=0; i<out.nrows(); ++i)
        {
          R d=one(i,i)*other(i,i);
          for (std::size_t j=0; j<other.ncols(); ++j)
            out(i,j)=d;
        }
      return out;
    }
  else if ((other.type()==M_Matrix<S>::SCALAR_DIAGONAL))
    {
      return M_Matrix<R>
          (one.nrows(),other.ncols(),M_Matrix<R>::SCALAR_FULL,
           one[0]*other[0]);
    }
  else if (other.type()==M_Matrix<S>::SCALAR_FULL)
    {
      return  M_Matrix<R>
          (one.nrows(),other.ncols(),
           M_Matrix<R>::SCALAR_FULL,one[0]*other[0]*one.ncols());
    }
  else // ZERO!
    {
      return M_Matrix<R>
          (one.nrows(),other.ncols(),M_Matrix<R>::ZERO);

    }

}

else //ZERO
{
return M_Matrix<R>
(one.nrows(),other.ncols(),M_Matrix<R>::ZERO);
}
}



template<typename T, typename S>
auto
quadraticForm_Symmetric (const M_Matrix<T>& A, const M_Matrix<S>& B)
->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
{
  typedef   decltype(std::declval<T>()*std::declval<S>()) R;

  auto AB=A*B;
  M_Matrix<R> out(B.ncols(),B.ncols(),M_Matrix<R>::SYMMETRIC, R(0));
  for (std::size_t i=0; i<out.nrows(); ++i)
    for (std::size_t j=i; j<out.ncols(); ++j)
      for (std::size_t k=0; k<B.nrows(); ++k)
        out+=B(k,i)*AB(k,j);
  return out;

}


template<typename T, typename S>
auto
quadraticForm_BT_A_B (const M_Matrix<T>& A, const M_Matrix<S>& B)
->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
{
  assert(A.ncols()==B.nrows());
  assert(B.ncols()==A.nrows());
  typedef   decltype(std::declval<T>()*std::declval<S>()) R;
  if ((A.type()==M_Matrix<T>::ZERO)||
      (B.type()==M_Matrix<S>::ZERO))
    return M_Matrix<R>(B.ncols(),B.ncols(),M_Matrix<R>::ZERO);
  else if (A.isSymmetric())
    return quadraticForm_Symmetric(A,B);
  else
    return transpMult(B,A*B);
}






template<class E>
class FullPartition
{

  FullPartition(const M_Matrix<E>& x, std::size_t n):
    A_(n,n),
    B_(n,x.ncols()-n),
    C_(x.nrows()-n,n),
    D_(x.nrows()-n,x.ncols()-n)
  {
    for (std::size_t i=0; i<n; ++i)
      {
        for (std::size_t j=0; j<n; ++j)
          A_(i,j)=x(i,j);
        for (std::size_t j=n; j<x.ncols(); ++j)
          B_(i,j-n)=x(i,j);
      }
    for (std::size_t i=n; i<x.nrows(); ++i)
      {
        for (std::size_t j=0; j<n; ++j)
          C_(i-n,j)=x(i,j);
        for (std::size_t j=n; j<x.nrows(); ++j)
          D_(i-n,j-n)=x(i,j);
      }

  }

  FullPartition(const M_Matrix<E>& A,const M_Matrix<E>& B,const M_Matrix<E>& C,const M_Matrix<E>& D)
    :A_(A),B_(B),C_(C),D_(D){}
  FullPartition( M_Matrix<E>&& A, M_Matrix<E>&& B, M_Matrix<E>&& C, M_Matrix<E>&& D)
    :A_(A),B_(B),C_(C),D_(D){}

  M_Matrix<E> full()const
  {
    std::size_t N=A_.nrows()+D_.nrows();
    M_Matrix<E> out(N,N);
    std::size_t n=A_.nrows();
    for (std::size_t i=0; i<n; ++i)
      {
        for (std::size_t j=0; j<n; ++j)
          out(i,j)=A_(i,j);
        for (std::size_t j=n; j<out.ncols(); ++j)
          out(i,j)=B_(i,j-n);
      }
    for (std::size_t i=n; i<out.nrows(); ++i)
      {
        for (std::size_t j=0; j<n; ++j)
          out(i,j)=C_(i-n,j);
        for (std::size_t j=n; j<out.nrows(); ++j)
          out(i,j)=D_(i-n,j-n);
      }
    return out;

  }

  FullPartition inverse_by_A()const
  {
    auto invA=inv(A_);
    auto AinvB=invA*B_;
    auto Dinv=inv(D_-C_*AinvB);
    auto CAinv=C_*invA;
    auto Binv=-AinvB*Dinv;
    auto Cinv=-Dinv*CAinv;
    auto Ainv=invA+AinvB*Dinv*CAinv;
    return FullPartition(Ainv,Binv,Cinv,Dinv);
  }
  FullPartition inverse_by_D()const
  {
    auto invD=inv(D_);
    auto DinvC=invD*C_;
    auto Ainv=inv(A_-B_*DinvC);
    auto BDinv=B_*invD;
    auto Cinv=-DinvC*Ainv;
    auto Binv=-Ainv*BDinv;
    auto Dinv=invD+DinvC*Ainv*BDinv;
    return FullPartition(Ainv,Binv,Cinv,Dinv);
  }


private:
  M_Matrix<E>& A_;
  M_Matrix<E>& B_;
  M_Matrix<E>& C_;
  M_Matrix<E>& D_;
};

template<class E>
class SymmetricPartition
{

  SymmetricPartition(const M_Matrix<E>& x, std::size_t n):
    A_(n,n,M_Matrix<E>::SYMMETRIC),
    B_(n,x.ncols()-n),
    D_(x.nrows()-n,x.ncols()-n,M_Matrix<E>::SYMMETRIC)
  {
    for (std::size_t j=0; j<n; ++j)
      {
        for (std::size_t i=j; i<n; ++i)
          A_(i,j)=x(i,j);
      }
    for (std::size_t j=n; j<x.nrows(); ++j)
      {
        for (std::size_t i=0; i<n; ++i)
          B_(i,j-n)=x(i,j);
        for (std::size_t i=j; i<x.nrows(); ++i)
          D_(i-n,j-n)=x(i,j);
      }

  }

  SymmetricPartition(const M_Matrix<E>& A,const M_Matrix<E>& B,const M_Matrix<E>& D)
    :A_(A),B_(B),D_(D){}
  SymmetricPartition( M_Matrix<E>&& A, M_Matrix<E>&& B,  M_Matrix<E>&& D)
    :A_(A),B_(B),D_(D){}

  M_Matrix<E> full()const
  {
    std::size_t N=A_.nrows()+D_.nrows();
    M_Matrix<E> x(N,N, M_Matrix<E>::SYMMETRIC);
    std::size_t n=A_.nrows();
    for (std::size_t j=0; j<n; ++j)
      {
        for (std::size_t i=j; i<n; ++i)
          x(i,j)=A_(i,j);
      }
    for (std::size_t j=n; j<x.nrows(); ++j)
      {
        for (std::size_t i=0; i<n; ++i)
          x(i,j)=B_(i,j-n);
        for (std::size_t i=j; i<x.nrows(); ++i)
          x(i,j)=D_(i-n,j-n);
      }
    return x;

  }

  SymmetricPartition inverse_by_A()const
  {
    auto invA=inv(A_);
    auto BTAinvB=quadraticForm_BT_A_B(invA,B_);
    auto Dinv=inv(D_-BTAinvB);
    auto AinvB=invA*B_;
    auto Binv=-AinvB*Dinv;
    auto Ainv=invA
        +quadraticForm_BT_A_B(Dinv,Transpose(AinvB));
    return SymmetricPartition(Ainv,Binv,Dinv);
  }
  SymmetricPartition inverse_by_D()const
  {
    auto invD=inv(D_);
    auto BDinvBT=
        quadraticForm_BT_A_B(invD,Transpose(B_));
    auto Ainv=inv(A_-BDinvBT);
    auto BDinv=B_*invD;
    auto Binv=-Ainv*BDinv;
    auto Dinv=invD+
        +quadraticForm_BT_A_B(Ainv,BDinv);
    return SymmetricPartition(Ainv,Binv,Dinv);
  }


private:
  M_Matrix<E>& A_;
  M_Matrix<E>& B_;
  M_Matrix<E>& D_;
};



template<typename E>
M_Matrix<E>
Matrix_inverse(const M_Matrix<E> x)
{
  assert(x.nrows()==x.ncols());
  switch(x.type())
    {
    case M_Matrix<E>::FULL:
      {
        auto p=FullPartition<E>(x,x.nrows()/2);
        auto pinv=p.inverse_by_A();
        return pinv.full();
      }
    case M_Matrix<E>::SYMMETRIC:
      {
        auto p=SymmetricPartition<E>(x,x.nrows()/2);
        auto pinv=p.inverse_by_A();
        return pinv.full();
      }
    case M_Matrix<E>::DIAGONAL:
      {
        M_Matrix<E> out(x.nrows(),x.nrows(),M_Matrix<E>::DIAGONAL);
        for (std::size_t i=0; i<x.nrows(); ++i)
          out(i,i)=inv(x(i,i));
        return out;
      }
    case M_Matrix<E>::SCALAR_DIAGONAL:
      {
        M_Matrix<E> out(x.nrows(),x.nrows(),M_Matrix<E>::SCALAR_DIAGONAL);
        out(0,0)=inv(x(0,0));
        return out;
      }
    case M_Matrix<E>::SCALAR_FULL:
    case M_Matrix<E>::ZERO:
      {
        assert(false);
        return {};
      }


    }
}




inline
M_Matrix<double>
Matrix_Mult_double (const M_Matrix<double>& one, const M_Matrix<double>& other)
{
  assert(one.ncols()==other.nrows());
  if (one.type()==M_Matrix<double>::FULL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        return Lapack_Full_Product(one,other,false,false);
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          auto tr=Transpose(one);
          return Lapack_Symmetric_Transpose_Product(other,tr,false);
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              for (std::size_t j=0; j<one.ncols(); ++j)
                out(i,j)=one(i,j)*other(j,j);
              for (std::size_t j=one.ncols(); j<out.ncols();++j)
                out(i,j)=0;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double sum=0;
              for (std::size_t j=0; j<one.ncols(); ++j)
                sum+=one(i,j);
              sum*=other[0];
              for (std::size_t k=0; k<other.ncols(); ++k)
                out(i,k)=sum;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }
  else   if (one.type()==M_Matrix<double>::SYMMETRIC)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          auto tr=Transpose(other);
          return Lapack_Symmetric_Transpose_Product(one,tr,true);
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          auto full=other.full();
          return Lapack_Symmetric_Transpose_Product(one,full,true);
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out
              (one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              for (std::size_t j=0; j<out.ncols(); ++j)
                out(i,j)=one(i,j)*other(j,j);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double sum=one(i,i);
              for (std::size_t j=0; j<one.ncols(); ++j)
                sum+=one(i,j);
              sum*=other[0];
              for (std::size_t k=0; k<other.ncols(); ++k)
                out(i,k)=sum;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }
  else  if (one.type()==M_Matrix<double>::DIAGONAL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          for (std::size_t i=one.ncols(); i<one.ncols();++i)
            for (std::size_t j=0; j<other.ncols(); ++j)
              out(i,j)=0;

          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out
              (one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out
              (M_Matrix<double>::DIAGONAL,one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              out(i,i)=one(i,i)*other(i,i);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        { M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double d=one(i,i)*other(i,i);
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=d;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }

  else  if (one.type()==M_Matrix<double>::SCALAR_DIAGONAL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          for (std::size_t i=one.ncols(); i<one.ncols();++i)
            for (std::size_t j=0; j<other.ncols(); ++j)
              out(i,j)=0;

          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out
              (one.nrows(),other.ncols(),M_Matrix<double>::SYMMETRIC);
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::DIAGONAL)

        {
          M_Matrix<double> out
              (one.nrows(),other.ncols(),M_Matrix<double>::DIAGONAL);
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              out(i,i)=one(i,i)*other(i,i);
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::SCALAR_DIAGONAL,
               one[0]*other[0]);
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          return  M_Matrix<double>
              (one.nrows(),other.ncols(),
               M_Matrix<double>::SCALAR_FULL,one[0]*other[0]);
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }

  else  if (one.type()==M_Matrix<double>::SCALAR_FULL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t k=0; k<other.ncols(); ++k)
            {
              double sum=0;
              for (std::size_t j=0; j<one.ncols(); ++j)
                sum+=other(j,k);
              sum*=one[0];
              for (std::size_t i=0; i<out.nrows(); ++i)
                out(i,k)=sum;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              double sum=other(i,i);
              for (std::size_t j=0; j<i; ++j)
                sum+=other(i,j);
              sum*=one[0];
              for (std::size_t k=0; k<one.nrows(); ++k)
                out(k,i)=sum;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::DIAGONAL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double d=one(i,i)*other(i,i);
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=d;
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::SCALAR_FULL,
               one[0]*other[0]);
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          return  M_Matrix<double>
              (one.nrows(),other.ncols(),
               M_Matrix<double>::SCALAR_FULL,one[0]*other[0]*one.ncols());
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }

  else //ZERO
    {
      return M_Matrix<double>
          (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);
    }
}


inline
M_Matrix<double>
operator * (const M_Matrix<double>& one, const M_Matrix<double>& other)
{
  return Matrix_Mult_double(one,other);
}

inline
M_Matrix<double>
multTransp (const M_Matrix<double>& one, const M_Matrix<double>& other)
{
  assert(one.ncols()==other.ncols());
  if (one.type()==M_Matrix<double>::FULL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        return Lapack_Full_Product(one,other,false,true);
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          auto tr=Transpose(one);
          return Lapack_Symmetric_Transpose_Product(other,tr,false);
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              for (std::size_t j=0; j<one.ncols(); ++j)
                out(i,j)=one(i,j)*other(j,j);
              for (std::size_t j=one.ncols(); j<out.nrows();++j)
                out(i,j)=0;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double sum=0;
              for (std::size_t j=0; j<one.ncols(); ++j)
                sum+=one(i,j);
              sum*=other[0];
              for (std::size_t k=0; k<other.nrows(); ++k)
                out(i,k)=sum;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.nrows(),M_Matrix<double>::ZERO);

        }

    }
  else   if (one.type()==M_Matrix<double>::SYMMETRIC)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          return Lapack_Symmetric_Transpose_Product(one,other,true);
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          auto full=other.full();
          return Lapack_Symmetric_Transpose_Product(one,full,true);
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out
              (M_Matrix<double>::SYMMETRIC,one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,j)*other(j,j);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double sum=one(i,i);
              for (std::size_t j=0; j<i; ++j)
                sum+=one(i,j);
              sum*=other[0];
              for (std::size_t k=0; k<other.ncols(); ++k)
                out(i,k)=sum;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (M_Matrix<double>::FULL,one.nrows(),other.ncols());

        }

    }
  else  if (one.type()==M_Matrix<double>::DIAGONAL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t i=0; i<other.ncols(); ++i)
            {
              for (std::size_t j=0; j<other.nrows(); ++j)
                out(i,j)=one(i,i)*other(j,i);
            }
          for (std::size_t i=one.ncols(); i<one.ncols();++i)
            for (std::size_t j=0; j<other.nrows(); ++j)
              out(i,j)=0;

          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out
              (M_Matrix<double>::SYMMETRIC,one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out
              (M_Matrix<double>::DIAGONAL,one.nrows(),other.nrows());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              out(i,i)=one(i,i)*other(i,i);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        { M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double d=one(i,i)*other(i,i);
              for (std::size_t j=0; j<other.nrows(); ++j)
                out(i,j)=d;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (M_Matrix<double>::ZERO,one.nrows(),other.nrows());

        }

    }

  else  if (one.type()==M_Matrix<double>::SCALAR_DIAGONAL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t i=0; i<other.ncols(); ++i)
            {
              for (std::size_t j=0; j<other.nrows(); ++j)
                out(i,j)=one(i,i)*other(j,i);
            }
          for (std::size_t i=one.ncols(); i<one.ncols();++i)
            for (std::size_t j=0; j<other.nrows(); ++j)
              out(i,j)=0;

          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out
              (M_Matrix<double>::SYMMETRIC,one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::DIAGONAL)

        {
          M_Matrix<double> out
              (M_Matrix<double>::DIAGONAL,one.nrows(),other.nrows());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              out(i,i)=one(i,i)*other(i,i);
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          return M_Matrix<double>
              (one.nrows(),other.nrows(),M_Matrix<double>::SCALAR_DIAGONAL,
               one[0]*other[0]);
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          return  M_Matrix<double>
              (one.nrows(),other.nrows(),
               M_Matrix<double>::SCALAR_FULL,one[0]*other[0]);
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (M_Matrix<double>::ZERO,one.nrows(),other.nrows());

        }

    }

  else  if (one.type()==M_Matrix<double>::SCALAR_FULL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t k=0; k<other.nrows(); ++k)
            {
              double sum=0;
              for (std::size_t j=0; j<one.ncols(); ++j)
                sum+=other(k,j);
              sum*=one[0];
              for (std::size_t i=0; i<out.ncols(); ++i)
                out(i,k)=sum;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out(one.nrows(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              double sum=other(i,i);
              for (std::size_t j=0; j<i; ++j)
                sum+=other(i,j);
              sum*=one[0];
              for (std::size_t k=0; k<one.nrows(); ++k)
                out(k,i)=sum;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::DIAGONAL)
        {
          M_Matrix<double> out(one.nrows(),other.nrows());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double d=one(i,i)*other(i,i);
              for (std::size_t j=0; j<other.nrows(); ++j)
                out(i,j)=d;
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          return M_Matrix<double>
              (one.nrows(),other.nrows(),M_Matrix<double>::SCALAR_FULL,
               one[0]*other[0]);
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          return  M_Matrix<double>
              (one.nrows(),other.nrows(),
               M_Matrix<double>::SCALAR_FULL,one[0]*other[0]*one.ncols());
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.nrows(),other.nrows(),M_Matrix<double>::ZERO);

        }

    }

  else //ZERO
    {
      return M_Matrix<double>
          (one.nrows(),other.nrows(),M_Matrix<double>::ZERO);
    }
}


inline
M_Matrix<double>
transpMult (const M_Matrix<double>& one, const M_Matrix<double>& other)
{
  assert(one.nrows()==other.nrows());
  if (one.type()==M_Matrix<double>::FULL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        return Lapack_Full_Product(one,other,true,false);
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          return Lapack_Symmetric_Transpose_Product(other,one,false);
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              for (std::size_t j=0; j<one.nrows(); ++j)
                out(i,j)=one(j,i)*other(j,j);
              for (std::size_t j=one.nrows(); j<out.ncols();++j)
                out(i,j)=0;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double sum=0;
              for (std::size_t j=0; j<one.nrows(); ++j)
                sum+=one(j,i);
              sum*=other[0];
              for (std::size_t k=0; k<other.ncols(); ++k)
                out(i,k)=sum;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.ncols(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }
  else   if (one.type()==M_Matrix<double>::SYMMETRIC)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          auto tr=Transpose(other);
          return Lapack_Symmetric_Transpose_Product(one,tr,true);
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          auto full=other.full();
          return Lapack_Symmetric_Transpose_Product(one,full,true);
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out
              (M_Matrix<double>::SYMMETRIC,one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(j,i)*other(j,j);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double sum=one(i,i);
              for (std::size_t j=0; j<i; ++j)
                sum+=one(j,i);
              sum*=other[0];
              for (std::size_t k=0; k<other.ncols(); ++k)
                out(i,k)=sum;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (M_Matrix<double>::FULL,one.ncols(),other.ncols());

        }

    }
  else  if (one.type()==M_Matrix<double>::DIAGONAL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          for (std::size_t i=one.nrows(); i<one.nrows();++i)
            for (std::size_t j=0; j<other.ncols(); ++j)
              out(i,j)=0;

          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out
              (M_Matrix<double>::SYMMETRIC,one.ncols(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::DIAGONAL)
               || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          M_Matrix<double> out
              (M_Matrix<double>::DIAGONAL,one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              out(i,i)=one(i,i)*other(i,i);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        { M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double d=one(i,i)*other(i,i);
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=d;
            }
          return out;
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (M_Matrix<double>::ZERO,one.ncols(),other.ncols());

        }

    }

  else  if (one.type()==M_Matrix<double>::SCALAR_DIAGONAL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          for (std::size_t i=one.nrows(); i<one.nrows();++i)
            for (std::size_t j=0; j<other.ncols(); ++j)
              out(i,j)=0;

          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out
              (M_Matrix<double>::SYMMETRIC,one.ncols(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              for (std::size_t j=0; j<=i; ++j)
                out(i,j)=one(i,i)*other(i,j);
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::DIAGONAL)

        {
          M_Matrix<double> out
              (M_Matrix<double>::DIAGONAL,one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              out(i,i)=one(i,i)*other(i,i);
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          return M_Matrix<double>
              (one.ncols(),other.ncols(),M_Matrix<double>::SCALAR_DIAGONAL,
               one[0]*other[0]);
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          return  M_Matrix<double>
              (one.ncols(),other.ncols(),
               M_Matrix<double>::SCALAR_FULL,one[0]*other[0]);
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (M_Matrix<double>::ZERO,one.ncols(),other.ncols());

        }

    }

  else  if (one.type()==M_Matrix<double>::SCALAR_FULL)
    {
      if (other.type()==M_Matrix<double>::FULL)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t k=0; k<other.ncols(); ++k)
            {
              double sum=0;
              for (std::size_t j=0; j<one.nrows(); ++j)
                sum+=other(j,k);
              sum*=one[0];
              for (std::size_t i=0; i<out.nrows(); ++i)
                out(i,k)=sum;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::SYMMETRIC)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<other.nrows(); ++i)
            {
              double sum=other(i,i);
              for (std::size_t j=0; j<i; ++j)
                sum+=other(i,j);
              sum*=one[0];
              for (std::size_t k=0; k<one.ncols(); ++k)
                out(k,i)=sum;
            }
          return out;
        }
      else if (other.type()==M_Matrix<double>::DIAGONAL)
        {
          M_Matrix<double> out(one.ncols(),other.ncols());
          for (std::size_t i=0; i<out.nrows(); ++i)
            {
              double d=one(i,i)*other(i,i);
              for (std::size_t j=0; j<other.ncols(); ++j)
                out(i,j)=d;
            }
          return out;
        }
      else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
        {
          return M_Matrix<double>
              (one.ncols(),other.ncols(),M_Matrix<double>::SCALAR_FULL,
               one[0]*other[0]);
        }
      else if (other.type()==M_Matrix<double>::SCALAR_FULL)
        {
          return  M_Matrix<double>
              (one.ncols(),other.ncols(),
               M_Matrix<double>::SCALAR_FULL,one[0]*other[0]*one.nrows());
        }
      else // ZERO!
        {
          return M_Matrix<double>
              (one.ncols(),other.ncols(),M_Matrix<double>::ZERO);

        }

    }

  else //ZERO
    {
      return M_Matrix<double>
          (one.ncols(),other.ncols(),M_Matrix<double>::ZERO);
    }
}


/**
     Transpose the first and multiply by the second
     @post transpMult(x,y)==Transpose(x)*y
     @remarks It is faster, since we save copying matrices
    */
inline
M_Matrix<double>
Lapack_TranspMult(const M_Matrix<double>& x,const M_Matrix<double>& y)
{
  assert(nrows(x)==nrows(y));
  // First it has to find out if the last dimension of x matches the first of y
  // now we build the M_Matrix result
  M_Matrix<double> z(x.ncols(),y.ncols(),0.0);

  /***  as fortran uses the reverse order for matrices and we want to
              avoid a copying operation, we calculate
                  Transpose(Z)=Transpose(y)*Transpose(x)

                  Transpose(matrix)=just plain matrix in C++ format


              */
  char  TRANSA='N';
  char 	TRANSB='T';
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
inline
M_Matrix<double> Lapack_multTransp(const M_Matrix<double>& x,const M_Matrix<double>& y)
{
  // First it has to find out if the last dimension of x matches the first of y
  //ASSERT_NE(x.size(),0);
  //ASSERT_EQ(x.ncols(),ncols(y));
  // now we build the M_Matrix result
  M_Matrix<double> z(x.nrows(),y.nrows(),0.0);
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
inline M_Matrix<double> Lapack_Full_Product(const M_Matrix<double>& x,const M_Matrix<double>& y, bool transpose_x, bool transpose_y)
{
  std::size_t cols_i, rows_i, rows_e,cols_e;
  if (transpose_x)
    {
      rows_e=x.ncols();
      cols_i=x.nrows();
    }
  else
    {
      rows_e=x.nrows();
      cols_i=x.ncols();
    }

  if (transpose_y)
    {
      rows_i=y.ncols();
      cols_e=y.nrows();
    }
  else
    {
      rows_i=y.nrows();
      cols_e=y.ncols();
    }


  assert(rows_i==cols_i);
  // First it has to find out if the last dimension of x matches the
  //first of y
  // now we build the M_Matrix result
  M_Matrix<double> z(rows_e,cols_e,0.0);

  /***  as fortran uses the reverse order for matrices and we want to
          avoid a copying operation, we calculate
              Transpose(Z)=Transpose(y)*Transpose(x)

              Transpose(matrix)=just plain matrix in C++ format


          */
  char TRANSA;
  char TRANSB;

  if (transpose_y)
    TRANSA='T';
  else
    TRANSA='N';

  if (transpose_x)
    TRANSB='T';
  else
    TRANSB='N';

  int  	M=cols_e;
  int  	N=rows_e;
  int  	K=cols_i;

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



inline M_Matrix<double> Lapack_Symmetric_Transpose_Product(const M_Matrix<double>& Sym,const M_Matrix<double>& Reg, bool SymRegT)
{
  /**
   *
   * SSYMM

Purpose:

     SSYMM  performs one of the matrix-matrix operations

        C := alpha*A*B + beta*C,

     or

        C := alpha*B*A + beta*C,

     where alpha and beta are scalars,  A is a symmetric matrix and  B and
     C are  m by n matrices.

Parameters
    [in]	SIDE

              SIDE is CHARACTER*1
               On entry,  SIDE  specifies whether  the  symmetric matrix  A
               appears on the  left or right  in the  operation as follows:

                  SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,

                  SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,

    [in]	UPLO

              UPLO is CHARACTER*1
               On  entry,   UPLO  specifies  whether  the  upper  or  lower
               triangular  part  of  the  symmetric  matrix   A  is  to  be
               referenced as follows:

                  UPLO = 'U' or 'u'   Only the upper triangular part of the
                                      symmetric matrix is to be referenced.

                  UPLO = 'L' or 'l'   Only the lower triangular part of the
                                      symmetric matrix is to be referenced.


    [in]	LDA

              LDA is INTEGER
               On entry, LDA specifies the first dimension of A as declared
               in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
               LDA must be at least  max( 1, m ), otherwise  LDA must be at
               least  max( 1, n ).

    [in]	B

              B is REAL array, dimension ( LDB, N )
               Before entry, the leading  m by n part of the array  B  must
               contain the matrix B.

    [in]	LDB

              LDB is INTEGER
               On entry, LDB specifies the first dimension of B as declared
               in  the  calling  (sub)  program.   LDB  must  be  at  least
               max( 1, m ).

    [in]	BETA

              BETA is REAL
               On entry,  BETA  specifies the scalar  beta.  When  BETA  is
               supplied as zero then C need not be set on input.

    [in,out]	C

              C is REAL array, dimension ( LDC, N )
               Before entry, the leading  m by n  part of the array  C must
               contain the matrix  C,  except when  beta  is zero, in which
               case C need not be set on entry.
               On exit, the array  C  is overwritten by the  m by n updated
               matrix.

    [in]	LDC

              LDC is INTEGER
               On entry, LDC specifies the first dimension of C as declared
               in  the  calling  (sub)  program.   LDC  must  be  at  least
               max( 1, m ).


   *
   * */
  /**
     * @brief UPLO
      [in]	UPLO

                UPLO is CHARACTER*1
                 On  entry,   UPLO  specifies  whether  the  upper  or  lower
                 triangular  part  of  the  symmetric  matrix   A  is  to  be
                 referenced as follows:

                    UPLO = 'U' or 'u'   Only the upper triangular part of the
                                        symmetric matrix is to be referenced.

                    UPLO = 'L' or 'l'   Only the lower triangular part of the
                                        symmetric matrix is to be referenced.
     */
  char  	UPLO='L';

  M_Matrix<double>  S=M_Matrix<double>::unpackForLapack(Sym,UPLO);
  std::size_t rows_e;
  std::size_t cols_i;
  std::size_t rows_i;
  std::size_t cols_e;

  if (SymRegT)  // calculo Sym * Reg^T
    {
      rows_e=Sym.nrows();
      cols_i=Sym.ncols();
      rows_i=Reg.ncols();
      cols_e=Reg.nrows();

    }
  else   // calculo  Reg^T * Sym
    {
      rows_e=Reg.ncols();
      cols_i=Reg.nrows();
      rows_i=Sym.ncols();
      cols_e=Sym.nrows();
    }
  assert(rows_i==cols_i);
  // First it has to find out if the last dimension of x matches the
  //first of y
  // now we build the M_Matrix result
  M_Matrix<double> z(rows_e,cols_e,0.0);

  /***  as fortran uses the reverse order for matrices and we want to
          avoid a copying operation, we calculate
              Transpose(Z)=Transpose(y)*Transpose(x)

              Transpose(matrix)=just plain matrix in C++ format



         */
  /**
   * @brief SIDE
   *   [in]	SIDE

              SIDE is CHARACTER*1
               On entry,  SIDE  specifies whether  the  symmetric matrix  A
               appears on the  left or right  in the  operation as follows:

                  SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,

                  SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,


   */
  char  	SIDE;
  if (SymRegT)
    SIDE='R';
  else
    SIDE='L';

  /**
   * @brief M
    [in]	M

              M is INTEGER
               On entry,  M  specifies the number of rows of the matrix  C.
               M  must be at least zero.
   */

  int  	M=cols_e;

  /**
   * @brief N
    [in]	N

              N is INTEGER
               On entry, N specifies the number of columns of the matrix C.
               N  must be at least zero.

               como uso la transpuesta de C, es rows_e

   */

  int  	N=rows_e;

  /**
   * @brief ALPHA
    [in]	ALPHA

              ALPHA is REAL
               On entry, ALPHA specifies the scalar alpha.
   */
  double  	ALPHA=1.0;

  /**
   * @brief A
    [in]	A

              A is REAL array, dimension ( LDA, ka ), where ka is
               m  when  SIDE = 'L' or 'l'  and is  n otherwise.
               Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
               the array  A  must contain the  symmetric matrix,  such that
               when  UPLO = 'U' or 'u', the leading m by m upper triangular
               part of the array  A  must contain the upper triangular part
               of the  symmetric matrix and the  strictly  lower triangular
               part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
               the leading  m by m  lower triangular part  of the  array  A
               must  contain  the  lower triangular part  of the  symmetric
               matrix and the  strictly upper triangular part of  A  is not
               referenced.
               Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
               the array  A  must contain the  symmetric matrix,  such that
               when  UPLO = 'U' or 'u', the leading n by n upper triangular
               part of the array  A  must contain the upper triangular part
               of the  symmetric matrix and the  strictly  lower triangular
               part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
               the leading  n by n  lower triangular part  of the  array  A
               must  contain  the  lower triangular part  of the  symmetric
               matrix and the  strictly upper triangular part of  A  is not
               referenced.

   */

  double * /*, dimension(lda,*)*/ A=const_cast<double*> (&S[0]);
  int  	LDA=S.nrows();
  double*  B=const_cast<double*> (&Reg[0]);

  int  	LDB=Reg.ncols();
  double   	BETA=0;
  double * /*, dimension(ldc,*) */ 	C;
  C=&z[0];
  int  	LDC=M;


  try
  {
    ssymm_(&SIDE,&UPLO,&M,&N,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
  }
  catch (...)
  {
    assert(false);
  }
  return z;
}



template<typename T, typename S>
auto Full_Product(const M_Matrix<T>& x, const M_Matrix<S>& y)
->
M_Matrix<std::result_of<decltype(std::declval<T>()*std::declval<S>())>>
{
                                                                        assert(x.ncols()==y.nrows());
                                                                        typedef   std::result_of<decltype(std::declval<T>()*std::declval<S>())> R;
                                                                        M_Matrix<R> out(x.nrows(),y.ncols(),R(0));
                                                                        for (std::size_t i=0; i<x.nrows(); ++i)
for (std::size_t j=0; j<x.ncols(); ++j)
for (std::size_t k=0; k<y.ncols(); ++k)
out(i,k)+=x(i,j)*y(j,k);
return out;


}




template<typename T>
M_Matrix<T> operator-(const M_Matrix<T>& x)
{
  M_Matrix<T> out(x);
  for (std::size_t i=0; i<out.size(); ++i)
    out[i]=-out[i];
  return out;
}

template<typename T>
M_Matrix<T> operator-(M_Matrix<T>&& x)
{
  for (std::size_t i=0; i<x.size(); ++i)
    x[i]=-x[i];
  return x;
}


template<template<typename...>class M,typename... Ts>
M<Ts...>& operator+=(M<Ts...>& itself, const M<Ts...>&  x)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]+=x[i];
  return itself;
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

template<template<typename...>class M,typename T, typename... Ts>
M<Ts...>& operator+=(M<Ts...>& itself, T x)
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
/*
 * template<typename T>
M_Matrix<T>& operator-=(M_Matrix<T>& itself, T x)
{
  for (size_t i=0; i<itself.size(); i++)
    itself[i]-=x;
  return itself;
}

*/

/*
Matrix Equality Operator
*/
template<typename T>
bool operator==(const M_Matrix<T>& x,
                const M_Matrix<T>& y)
{
  if (x.type()!=y.type()) return false;
  if (x.size()!=y.size()) return false;
  else if (x.ncols()!=y.ncols()) return false;
  else for (std::size_t i=0; i< x.size(); ++i)
    if (x[i]!=y[i]) return false;
  return true;

}



/*
 Minor Operator based on a Lexicographic Comparison.
*/
template<typename T>
bool operator<(const M_Matrix<T>& x, const M_Matrix<T>& y)
{
  if (x.size()<y.size()) return true;
  else if (y.size()<x.size()) return false;
  else if (x.nrows()<y.nrows()) return true;
  else if (y.nrows()<x.nrows()) return false;
  else for (std::size_t i=0; i< x.size(); ++x)
    {
      if (x[i]<y[i]) return true;
      else if (y[i]<x[i]) return false;
    }
  return false;

}


/**
     Scalar Multiplication assignment.
     @returns a reference to itself
     @post all the values of the matrix are multiplied by the value x
     */
template<typename T>
M_Matrix<T>& operator*=(M_Matrix<T>& itself, T x)
{
  if (itself.type()==M_Matrix<T>::ZERO)
    return itself;
  else
    {
      for (size_t i=0; i<itself.size(); i++)
        itself[i]*=x;
      return itself;
    }
}


/**
     Scalar Division assignment.
     @returns a reference to itself
     @post all the values of the matrix are divided by the value x
     */
template<typename T>
M_Matrix<T>& operator/=(M_Matrix<T>& itself, T x)
{
  if (itself.type()!=M_Matrix<T>::ZERO)
    for (size_t i=0; i<itself.size(); i++)
      itself[i]/=x;
  return itself;
}

//@}
/** @name Aritmetic Assigment Operations between two Matrices
      */
//@{



template<typename T, typename S, class AsOp>
M_Matrix<T>& Aditive_Assigment_Operator(AsOp op,M_Matrix<T>& itself,
                                        const M_Matrix<S>& other)
{
  if (itself.type()==M_Matrix<T>::ZERO)
    {
      itself=other;
      return itself;
    }
  if (itself.type()==other.type())
    {
      assert(itself.size()==other.size());
      for (size_t i=0; i<itself.size(); i++)
        op(itself[i],other[i]);
      return itself;
    }
  else if(itself.type()==M_Matrix<T>::FULL)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::FULL:  // should not land here
        case M_Matrix<S>::SYMMETRIC:
        case M_Matrix<S>::SCALAR_FULL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              for (std::size_t j=0; j<itself.ncols(); ++j)
                op(itself(i,j),other(i,j));
            return itself;
          }
        case M_Matrix<S>::SCALAR_DIAGONAL:
        case M_Matrix<S>::DIAGONAL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              op(itself(i,i),other(i,i));
            return itself;
          }
        case M_Matrix<S>::ZERO:
          if (op(T(1),S(0))==T(1))
            return itself;
          else
            {
              itself=S(0);
              return itself;
            }

        }
    }
  else if(other.type()==M_Matrix<T>::FULL)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      M_Matrix<T> out(other);
      Aditive_Assigment_Operator(op,out,itself);
      itself=std::move(out);
      return itself;
    }
  else if(itself.type()==M_Matrix<T>::SYMMETRIC)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::FULL:  //should not land here
        case M_Matrix<S>::SYMMETRIC: // should not land here
        case M_Matrix<S>::SCALAR_FULL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              for (std::size_t j=0; j<=i; ++j)
                op(itself(i,j),other(i,j));
            return itself;
          }
        case M_Matrix<S>::SCALAR_DIAGONAL:
        case M_Matrix<S>::DIAGONAL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              op(itself(i,i),other(i,i));
            return itself;
          }
        case M_Matrix<S>::ZERO:
          return itself;

        }
    }
  else if(other.type()==M_Matrix<T>::SYMMETRIC)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      M_Matrix<T> out(other);
      Aditive_Assigment_Operator(op,out,itself);
      itself=std::move(out);
      return itself;
    }
  else if(itself.type()==M_Matrix<T>::DIAGONAL)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::FULL:  //should not land here
        case M_Matrix<S>::SYMMETRIC: // should not land here
        case M_Matrix<S>::SCALAR_FULL:
          {
            M_Matrix<T> out(M_Matrix<T>::SYMMETRIC,itself.nrows(),itself.ncols(),other[0]);
            for (std::size_t i=0; i<itself.nrows(); ++i)
              op(other(i,i),itself(i,i));
            itself=std::move(other);
            return itself;
          }
        case M_Matrix<S>::SCALAR_DIAGONAL:
        case M_Matrix<S>::DIAGONAL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              op(itself(i,i),other(i,i));
            return itself;
          }
        case M_Matrix<S>::ZERO:
          return itself;

        }
    }
  else if(itself.type()==M_Matrix<T>::SCALAR_FULL)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::FULL:  //should not land here
        case M_Matrix<S>::SYMMETRIC: // should not land here
        case M_Matrix<S>::SCALAR_FULL: // should not land here
        case M_Matrix<S>::SCALAR_DIAGONAL:
        case M_Matrix<S>::DIAGONAL:
          {
            M_Matrix<T> out(M_Matrix<T>::SYMMETRIC,itself.nrows(),itself.ncols(),itself[0]);
            for (std::size_t i=0; i<std::min(itself.nrows(), itself.ncols()); ++i)
              op(out(i,i),other(i,i));
            itself=std::move(out);
            return itself;
          }
        case M_Matrix<S>::ZERO:
          return itself;

        }
    }
  else //SCALAR DIAGONAL
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::FULL:  //should not land here
        case M_Matrix<S>::SYMMETRIC: // should not land here
        case M_Matrix<S>::SCALAR_FULL:  //could be
          {
            M_Matrix<T> out(M_Matrix<T>::SYMMETRIC,itself.nrows(),itself.ncols(),other[0]);
            for (std::size_t i=0; i<itself.nrows(); ++i)
              op(other(i,i),itself(i,i));
            itself=std::move(other);
            return itself;
          }
        case M_Matrix<S>::SCALAR_DIAGONAL:  // should not be
        case M_Matrix<S>::DIAGONAL:  //could
          {
            M_Matrix<T> out(other);
            for (std::size_t i=0; i<std::min(itself.nrows(), itself.ncols()); ++i)
              op(out(i,i),other(i,i));
            itself=std::move(out);
            return itself;
          }
        case M_Matrix<S>::ZERO:
          return itself;

        }
    }
}


template<typename T, typename S, class AsOp>
M_Matrix<T>& Element_Wise_Multiplicative_Assigment_Operator
(AsOp op,M_Matrix<T>& itself,
 const M_Matrix<S>& other)
{
  if (itself.type()==M_Matrix<T>::ZERO)
    {
      return itself;
    }
  else if (other.type()==M_Matrix<S>::ZERO)
    {
      itself=other;
      return itself;
    }
  else if (itself.type()==other.type())
    {
      assert(itself.size()==other.size());
      for (size_t i=0; i<itself.size(); i++)
        op(itself[i],other[i]);
      return itself;
    }
  else if ((itself.type()==M_Matrix<T>::DIAGONAL)
           ||(itself.type()==M_Matrix<T>::SCALAR_DIAGONAL))
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::ZERO: // not possible
        case M_Matrix<S>::FULL:
        case M_Matrix<S>::SYMMETRIC:
        case M_Matrix<S>::SCALAR_FULL:
        case M_Matrix<S>::SCALAR_DIAGONAL:
        case M_Matrix<S>::DIAGONAL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              op(itself(i,i),other(i,i));
            return itself;
          }

        }
    }

  else if ((other.type()==M_Matrix<T>::DIAGONAL)
           ||(other.type()==M_Matrix<T>::SCALAR_DIAGONAL))
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(itself.type())
        {
        M_Matrix<T> out(other);
        for (std::size_t i=0; i<std::min(out.nrows(), out.ncols()); ++i)
          op(out(i,i),itself(i,i));
        itself=std::move(other);
        return itself;

        }
    }
  else if(itself.type()==M_Matrix<T>::FULL)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::ZERO: // not
        case M_Matrix<S>::SCALAR_DIAGONAL:  //done
        case M_Matrix<S>::DIAGONAL: // done
        case M_Matrix<S>::FULL:
        case M_Matrix<S>::SYMMETRIC:
        case M_Matrix<S>::SCALAR_FULL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              for (std::size_t j=0; j<itself.ncols(); ++j)
                op(itself(i,j),other(i,j));
            return itself;
          }

        }
    }
  else if(other.type()==M_Matrix<T>::FULL)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      M_Matrix<T> out(other);
      for (std::size_t i=0; i<itself.nrows(); ++i)
        for (std::size_t j=0; j<itself.ncols(); ++j)
          op(out(i,j),itself(i,j));
      itself=std::move(out);
      return itself;
    }
  else if(itself.type()==M_Matrix<T>::SYMMETRIC)
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      switch(other.type())
        {
        case M_Matrix<S>::ZERO:  //not possible
        case M_Matrix<S>::SCALAR_DIAGONAL: //not possible
        case M_Matrix<S>::DIAGONAL: //not possible
        case M_Matrix<S>::FULL:  //should not land here
        case M_Matrix<S>::SYMMETRIC: // should not land here
        case M_Matrix<S>::SCALAR_FULL:
          {
            for (std::size_t i=0; i<itself.nrows(); ++i)
              for (std::size_t j=0; j<=i; ++j)
                op(itself(i,j),other(i,j));
            return itself;
          }
        }
    }
  else  // scalar_full and the other is symmetric
    {
      assert(itself.nrows()==other.nrows());
      assert(itself.ncols()==other.ncols());
      M_Matrix<T> out(other);
      for (std::size_t i=0; i<itself.nrows(); ++i)
        for (std::size_t j=0; j<=i; ++j)
          op(out(i,j),itself(i,j));
      itself=std::move(out);
      return itself;
    }
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
template<typename T, typename S>
M_Matrix<T>& elemMult(M_Matrix<T>& x,const M_Matrix<S>& y)
{
  return Element_Wise_Multiplicative_Assigment_Operator([](T& a, const S& b){a*=b; return a;},x,y);
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
M_Matrix<double> operator<<(const M_Matrix<double>& A, const M_Matrix<T>& B)
{
  M_Matrix<double> out
      (std::max(A.nrows(), B.nrows()), A.ncols()+B.ncols(), std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i=0; i<A.nrows();++i)
    {
      for (std::size_t j=0; j<A.ncols(); ++j)
        out(i,j)=A(i,j);
    }
  for (std::size_t i=0; i<B.nrows();++i)
    {
      for (std::size_t j=0; j<B.ncols(); ++j)
        out(i,j)=B(i,A.ncols()+j);
    }

  return out;

}




template<typename T>
std::istream& operator>>(std::istream& is,M_Matrix<T>& x)
{
  std::vector<T> o;
  std::size_t nrows=0;
  char ch;
  while ((is>>ch)&&(ch!='[')){}
  if(ch!='[')
    return is;
  else
    while (ch!=']')
      {
        std::string s;
        while ((is.get(ch))&&((ch!=']')&&ch!=';'))
          {
            s.push_back(ch);
          }
        std::stringstream ss(s);
        T e;
        std::size_t i=o.size();
        while (ss>>e) o.push_back(e);
        if (o.size()>i) ++nrows;
      }
  std::size_t ncols=o.size()/nrows;
  x=M_Matrix<T>(nrows,ncols,o);
  return is;

}

template<typename T>
T maxAbs(const M_Matrix<T>& x)
{
  T m=std::abs(x[0]);
  for (std::size_t i=0; i<x.size(); ++i)
    if (std::abs(x[i])>m) m=std::abs(x[i]);
  return m;
}

template<typename T, class Compare>
M_Matrix<T> sort(const M_Matrix<T>& x, Compare comp)
{
  std::vector<T> o=x.toVector();
  std::sort(o.begin(), o.end(), comp);
  return M_Matrix<T>(x.nrows(),x.ncols(),o);

}
template<typename T>
M_Matrix<T> sort(const M_Matrix<T>& x)
{
  std::vector<T> o=x.toVector();
  std::sort(o.begin(), o.end());
  return M_Matrix<T>(x.nrows(),x.ncols(),o);

}


template <class V>
double mean(const V& x)
{
  double sum=0;
  for (std::size_t i=0; i<x.size(); ++i)
    sum+=x[i];
  return sum/x.size();
}
template <class V>
double stddev(const V& x)
{
  double sum=0;
  double sum2=0;
  for (std::size_t i=0; i<x.size(); ++i)
    {
      sum+=x[i];sum2+=sqr(x[i]);
    }
  return std::sqrt(sum2/x.size()-sqr(sum/x.size()));
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





inline M_Matrix<double> xdiagXT(const M_Matrix<double>& x, const M_Matrix<double> Cdiag)
{
  M_Matrix<double> o(x.nrows(), x.nrows(),0.0);
  for (std::size_t i=0;  i<x.nrows(); ++i)
    for (std::size_t j=0; j<x.nrows(); ++j)
      for (std::size_t k=0; k<x.ncols(); ++k)
        o(i,j)+=Cdiag[k]*x(i,k)*x(j,k);
  return o;
}



inline M_Matrix<double> MultDiag(const M_Matrix<double> &x, const M_Matrix<double> d)
{
  M_Matrix<double> o(x.nrows(), x.ncols());
  for (std::size_t i=0;  i<x.nrows(); ++i)
    for (std::size_t j=0; j<x.ncols(); ++j)
      o(i,j)=x(i,j)*d[j];
  return o;
}


inline M_Matrix<double> DiagMult( const M_Matrix<double> d,const M_Matrix<double> &x)
{
  M_Matrix<double> o(x.nrows(), x.ncols());
  for (std::size_t i=0;  i<x.nrows(); ++i)
    for (std::size_t j=0; j<x.ncols(); ++j)
      o(i,j)=x(i,j)*d[i];
  return o;
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
      M_Matrix<T> diagM(n,n,M_Matrix<T>::DIAGONAL);
      for (size_t i=0; i<n; ++i)
        diagM(i,i)=x(i,i);
      return diagM;
    }
  else
    {
      nr=std::max(nr,nc);
      M_Matrix<T> diagM(nr,nr,M_Matrix<T>::DIAGONAL);
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
T logDiagProduct(const M_Matrix<T>& x)
{
  size_t nr=x.nrows();
  size_t nc=x.ncols();
  double diagprod=0;
  std::size_t n=std::min(nr,nc);
  for (size_t i=0; i<n; ++i)
    diagprod+=std::log(x(i,i));
  return diagprod;
}

template<typename T>
T Trace(const M_Matrix<T>& x)
{
  size_t nr=x.nrows();
  size_t nc=x.ncols();
  T out=0;
  std::size_t n=std::min(nr,nc);
  for (size_t i=0; i<n; ++i)
    out+=x(i,i);
  return out;


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



template<typename V>
M_Matrix<double> JTd2J(const M_Matrix<double>& J, const V& D2)
{
  assert(J.nrows()==D2.size());
  M_Matrix<double> out(J.ncols(),J.ncols());
  M_Matrix<double> Jc=Transpose(J);
  for (std::size_t i=0; i<Jc.nrows(); ++i)
    for (std::size_t j=0; j<Jc.ncols(); ++j)
      Jc(i,j)*=D2[j];
  out=Jc*J;
  return out;

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

  if (LDA==0)
    return M_Matrix<double>();
  try
  {
    dpotrf_(&UPLO,&N,&res[0],&LDA,&INFO);
  }
  catch (...)
  {
    std::cerr<<" error";
  };
  if (INFO!=0)
    {
      std::cerr<< "wrong cholesky!";
      return {};
    }
  else
    return Transpose(res);

}




#endif // MATRIX_H
