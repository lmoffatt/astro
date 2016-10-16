#ifndef MATRIXINVERSE_H
#define MATRIXINVERSE_H
#include <vector>
#include <iostream>

std::vector< std::vector< double> >
inv(const std::vector< std::vector< double> >& matrix);

double det(const std::vector< std::vector< double> >& matrix);

std::vector< std::vector< double> >
chol(const std::vector< std::vector< double> >& matrix);

inline
std::vector<double>
mult(const std::vector<std::vector<double>>& A,
     const std::vector<double>& b)
{
  std::vector<double> o(A.size(),0);
  for (std::size_t i=0; i<A.size(); ++i)
    for (std::size_t j=0; j<b.size(); ++j)
      o[i]+=A[i][j]*b[j];
  return o;
}

inline
std::vector<double>
tr_mult(const std::vector<std::vector<double>>& A,
        const std::vector<double>& b)
{
  std::vector<double> o(A[0].size(),0);
  for (std::size_t i=0; i<A.size(); ++i)
    for (std::size_t j=0; j<b.size(); ++j)
      o[i]+=A[j][i]*b[j];
  return o;
}



inline
std::vector<std::vector<double>>
mult(const std::vector<std::vector<double>>& A,
     const std::vector<std::vector<double>>& B)
{
  std::vector<std::vector<double>> O(A.size(),std::vector<double>(B[0].size(),0));
  for (std::size_t i=0; i<A[0].size(); ++i)
    for (std::size_t j=0; j<B[0].size(); ++j)
      for (std::size_t k=0; k<B.size(); ++k)
        O[i][j]+=A[i][k]*B[k][j];
  return O;
}




template<typename T>
T min(const std::vector<T>& o)
{
  if (o.empty())
    return {};
  else
    {
      T s=o[0];
      for (std::size_t i=1; i<o.size(); ++i)
        if (o[i]<s) s=o[i];
      return s;
    }

}

template<typename T>
T max(const std::vector<T>& o)
{
  if (o.empty())
    return {};
  else
    {
      T s=o[0];
      for (std::size_t i=1; i<o.size(); ++i)
        if (o[i]>s) s=o[i];
      return s;
    }

}



std::ostream& operator<<(
    std::ostream& s,const std::vector< std::vector< double> >& matrix);


std::ostream& operator<<(
    std::ostream& s,const std::vector<  double> & aVector);


#endif // MATRIXINVERSE_H

