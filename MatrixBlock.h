#ifndef MATRIXBLOCK_H
#define MATRIXBLOCK_H
#include "Matrix.h"
#include "MatrixInverse.h"


template<typename T,class E>
class Matrix_Equal_Square_Blocks

{
public:
  Matrix_Equal_Square_Blocks(const M_Matrix<E>& v): e_(v){}
  size_t nrows()const
  {
    if (!empty())
    return e_.nrows()*e_[0].nrows();
    else return 0;
  }

  size_t ncols()const
  {
    if (!empty())
    return e_.ncols()*e_[0].nrows();
    else return 0;
  }

  std::size_t nblock_rows()const {return e_.nrows();}

  std::size_t nblock_cols()const {return e_.ncols();}


  std::size_t nPerBlock()const {if (!e_.empty()) return e_[0].size(); else return 0;}



  bool empty()const
  {
    return e_.empty();
  }



  T&  operator() (std::size_t i,std::size_t j)
  {
    assert(i<nrows());
    assert(j<ncols());
    std::size_t i_b=i/nPerBlock();
    std::size_t i_in_block=i-i_b*nPerBlock();

    std::size_t j_b=j/nPerBlock();
    std::size_t j_in_block=j-j_b*nPerBlock();

    return e_(i_b,j_b)(i_in_block,j_in_block);

  }

  const  T&  operator() (std::size_t i,std::size_t j) const
  {
    assert(i<nrows());
    assert(j<ncols());
    std::size_t i_b=i/nPerBlock();
    std::size_t i_in_block=i-i_b*nPerBlock();

    std::size_t j_b=j/nPerBlock();
    std::size_t j_in_block=j-j_b*nPerBlock();

    return e_(i_b,j_b)(i_in_block,j_in_block);

  }

  friend
  Matrix_Equal_Square_Blocks
  Transpose(const Matrix_Equal_Square_Blocks& me)
  {
    Matrix_Equal_Square_Blocks out(Transpose(me.e_));
  }

 private:
  M_Matrix<E> e_;
};


template<typename T>
class Diagonal_Matrix
{
public:
  Diagonal_Matrix(const std::vector<T>& d):diag_(d){}

  std::size_t nrows()const {return diag_.size();}

  std::size_t ncols()const {return nrows();}


  T& operator()(std::size_t i, std::size_t j)
  {
    assert(i==j);
    return diag_[i];
  }
  const T& operator()(std::size_t i, std::size_t j) const
  {
    if (i==j)
    return diag_[i];
    else return T(0);
  }
  template<class BiOp>
  Diagonal_Matrix& apply(BiOp biOp,const Diagonal_Matrix& other)
  {
    assert (nrows()==other.nrows());
    for (std::size_t i=0; i<nrows(); ++i)
      biOp((*this)(i,i),other(i,i));
    return *this;
  }

  Diagonal_Matrix& operator*=(const Diagonal_Matrix& other)
  {
    return apply([](T& me, const T& to){return me*=to;}, other);
  }

  Diagonal_Matrix& operator+=(const Diagonal_Matrix& other)
  {
    return apply([](T& me, const T& to){return me+=to;}, other);
  }


private:
  std::vector<T> diag_;
};


template<typename T>
M_Matrix<T>
operator*(const M_Matrix<T>& one, const Diagonal_Matrix<T>& two)
{
  assert(one.ncols()==two.nrows());
  M_Matrix<T> out(one);
  for (std::size_t i=0; i<one.nrows(); ++i)
    for (std::size_t j=0; j<one.ncols(); ++j)
      out(i,j)*=two(j,j);
  return out;
}

template<typename T>
M_Matrix<T>
operator*( const Diagonal_Matrix<T>& one,const M_Matrix<T>& two)
{
  assert(one.ncols()==two.nrows());
  M_Matrix<T> out(two);
  for (std::size_t i=0; i<two.nrows(); ++i)
    for (std::size_t j=0; j<two.ncols(); ++j)
      out(i,j)*=one(i,i);
  return out;
}
inline
double inv(double x){return 1.0/x;}

template<typename T>
Diagonal_Matrix<T>
inv(const Diagonal_Matrix<T>& me)
{
  Diagonal_Matrix<T> out(me);
  for (std::size_t i=0; i<me.nrows(); ++i)
    out(i,i)=inv(me(i,i));
}


template<typename T,class E>
class Diagonal_Equal_Blocks
{
public:
  Diagonal_Equal_Blocks(const std::vector<E>& v): e_(v){}
  size_t nrows()const
  {
    return e_.size()*e_[0].nrows();
  }

  size_t ncols()const
  {
    return nrows();
  }

  std::size_t nblocks()const {return e_.size();}

  std::size_t nPerBlock()const {if (!e_.empty()) return e_[0].size(); else return 0;}
  bool empty()const
  {
    return e_.empty();
  }


  E& diag(std::size_t i)
  {
    assert(i<nblocks());
    return e_[i];
  }

  E const & diag(std::size_t i)const
  {
    assert(i<nblocks());
    return e_[i];
  }


  T&  operator() (std::size_t i,std::size_t j)
  {
    std::size_t i_b=i/nPerBlock();
    std::size_t i_in_block=i-i_b*nPerBlock();

    std::size_t j_b=j/nPerBlock();
    std::size_t j_in_block=j-j_b*nPerBlock();

    assert (i_b==j_b);
    return e_[i_b](i_in_block,j_in_block);

  }

  const  T&  operator() (std::size_t i,std::size_t j) const
  {
    std::size_t i_b=i/nPerBlock();
    std::size_t i_in_block=i-i_b*nPerBlock();

    std::size_t j_b=j/nPerBlock();
    std::size_t j_in_block=j-j_b*nPerBlock();

    if  (i_b==j_b)
     return e_[i_b](i_in_block,j_in_block);
    else
      return T();

  }


  template<class BiOp>
  Diagonal_Equal_Blocks& apply(BiOp biOp,const Diagonal_Equal_Blocks& other)
  {
    assert (nblocks()==other.nblocks());
    assert (nPerBlock()==other.nPerBlock());
    for (std::size_t i=0; i<nblocks(); ++i)
      biOp(diag(i),other.diag(i));
    return *this;
  }
  Diagonal_Equal_Blocks& operator*=(const Diagonal_Equal_Blocks& other)
  {
    return apply([](E& me, const E& to){return me*=to;}, other);
  }

  Diagonal_Equal_Blocks& operator*=(const T& landa)
  {
    return apply([](E& me, const T& to){return me*=to;}, landa);
  }

  Diagonal_Equal_Blocks& operator+=(const Diagonal_Equal_Blocks& other)
  {
    return apply([](E& me, const E& to){return me+=to;}, other);
  }




 private:
  std::vector<E> e_;
};


template<typename T,class E>
Diagonal_Equal_Blocks<T,E>
operator+(const Diagonal_Equal_Blocks<T,E>& one,const Diagonal_Equal_Blocks<T,E>& two)
{
  Diagonal_Equal_Blocks<T,E> out(one);
  out+=two;
  return out;
}

template<typename T,class E>
Diagonal_Equal_Blocks<T,E>
operator*(const Diagonal_Equal_Blocks<T,E>& one,const Diagonal_Equal_Blocks<T,E>& two)
{
  Diagonal_Equal_Blocks<T,E> out(one);
  out*=two;
  return out;
}




template<typename T,class E>
Diagonal_Equal_Blocks<T,E>
inv(const Diagonal_Equal_Blocks<T,E>& one)
{
  std::vector<E> out;
  for (std::size_t i=0; i<one.nblocks(); ++i)
    out.push_back(inv(one.diag(i)));
  return Diagonal_Equal_Blocks<T,E>(out);
}



template<typename T, class E,
         template<typename, typename>class A,
         template<typename, typename>class B,
         template<typename, typename>class D
         >
class Symmetric_Square_Block
{
public:
  Symmetric_Square_Block(const A<T,E>& B11,const B<T,E>& B12,const D<T,E>& B22)
    :A_(B11),B_(B12),D_(B22){
    assert(A_.nrows()==A_.ncols());
    assert(D_.nrows()==D_.ncols());
    assert(B_.nrows()==A_.nrows());
    assert(B_.ncols()==D_.ncols());

  }

  std::size_t nrows()const{return A_.nrows()+D_.nrows();}
  std::size_t ncols()const {return nrows();}
  T& operator () (std::size_t i, std::size_t j)
  {
     if (i<A_.nrows())
       {
         if (j<A_.ncols())
           return A_(i,j);
         else
           return B_(i,j-A_.ncols());
       }
     else
       if (j<A_.ncols())
         return B_(j,i);
       else
         return D_(i-A_.nrows(),j-A_.ncols());
  }
  const T& operator () (std::size_t i, std::size_t j)const
  {
     if (i<A_.nrows())
       {
         if (j<A_.ncols())
           return A_(i,j);
         else
           return B_(i,j-A_.ncols());
       }
     else
       if (j<A_.ncols())
         return B_(j,i);
       else
         return D_(i-A_.nrows(),j-A_.ncols());
  }





private:
  A<T,E> A_;
  B<T,E> B_;
  D<T,E> D_;



};











#endif // MATRIXBLOCK_H
