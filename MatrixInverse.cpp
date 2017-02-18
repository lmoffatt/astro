#include "MatrixInverse.h"

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


std::vector< std::vector< double> >
inv(const std::vector< std::vector< double> >& matrix)
{
        double *A, *work, *work1;
        int info=0;
        //  char msg[101];
        int *ipiv;
        int lwork;
        int n =matrix.size();
        int m=n;
        int dla=n;
        //A=new double[n*n];
        A= new double[n*n]; //more efficient code
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                *(A+i+n*j) = matrix[i][j];

        ipiv = new int[n];
        dgetrf_(&n, &m, A, &dla,ipiv,&info);
        lwork= -1;
        work1 = new double[2];
        dgetri_(&n,A,&dla,ipiv,work1,&lwork,&info);
        lwork = (int)(work1[0]);
        work = new double [2*lwork];
        dgetri_(&n,A,&dla,ipiv,work,&lwork,&info);

        std::vector< std::vector<double> > result(n,std::vector<double>(n));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i][j]=*(A+i+n*j);
        delete [] ipiv;
        delete [] work;
        delete [] work1;
        delete [] A;

        if (info>0)
            std::cerr<<"\n singular matrix\n";

        return result;
    }


