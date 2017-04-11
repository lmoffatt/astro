#include "MatrixInverse.h"
#include <cmath>

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

extern "C" void dpotrf_(char * 	UPLO,
			int * N,
			double * A,
			int * LDA,
			int * INFO);

std::vector< std::vector< double> > UT(const std::vector< std::vector< double> >& x)
{
    std::vector< std::vector< double> > y(x.size(),std::vector<double>(x[0].size(),0));
    for (std::size_t i=0;i<x.size();i++)
    {
        for (std::size_t j=0;j<i; j++)
            y[i][j]=0;
        for (std::size_t j=i;j<x.size(); j++)
            y[i][j]=x[i][j];
    }
    return y;
}

std::vector< std::vector< double> > LT(const std::vector< std::vector< double> >& x)
{
    std::vector< std::vector< double> > y(x.size(),std::vector<double>(x[0].size(),0));
    for (std::size_t i=0;i<x.size();i++)
    {
        for (std::size_t j=0;j<i+1; j++)
            y[i][j]=x[i][j];
        for (std::size_t j=i+1;j<x[0].size(); j++)
            y[i][j]=0;
    }
    return y;

   }

std::vector< std::vector< double> >
chol(const std::vector< std::vector< double> >& x)
{
std::string kind="lower";
   int n =x.size();
    char UPLO='L';
   std::vector< std::vector< double> > res;
    if (kind!="lower")
    {
        UPLO='U';
        res=UT(x);
    }
    else
    {
        res=LT(x);
    }
    int N=x.size();
    int LDA=N;
    int INFO;
    double* A=new double[x.size()*x[0].size()];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            *(A+i+n*j) = res[i][j];



    dpotrf_(&UPLO,&N,A,&LDA,&INFO);


    std::vector< std::vector<double> > result(n,std::vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i][j]=*(A+i+n*j);

    delete [] A;
    return result;

}

double det(const std::vector<std::vector<double> > &x)
{
    std::vector< std::vector< double> > ch=chol(x);
    double d=1;
    for (std::size_t i=0; i<x.size(); i++)
    {
        d*=(ch[i][i]);
    }
       return d*d;
}


