#include "Distributions.h"





double MultivariateGaussian::P(const M_Matrix<double>& x)const
{
    return exp(logP(x));
}


std::size_t MultivariateGaussian::size()const
{
    return mean_.size();
}

M_Matrix<double> MultivariateGaussian::sample(std::mt19937_64 &mt)const
{
    M_Matrix<double> r;
    std::normal_distribution<> normal;
    if (this->size()>0)
    {
        M_Matrix<double> z(mean_.nrows(),mean_.ncols());
        for (std::size_t i=0; i<size(); i++)
            z[i]=normal(mt);
        r=mean_+(z*cho_cov_);
    }
    return r;
}

MultivariateGaussian::MultivariateGaussian(const M_Matrix<double> &mean,
                                           const M_Matrix<double>& cov):
    mean_(mean),
    cov_(cov),
    covinv_(invSafe(cov)),
    cho_cov_(chol(cov,"upper")),
    logDetCov_(log(diagProduct(cho_cov_)))
{}

MultivariateGaussian::MultivariateGaussian(const M_Matrix<double>&mean
                                           , const M_Matrix<double> &cov
                                           , const M_Matrix<double> &covInv):
mean_(mean),
cov_(cov),
covinv_(covInv),
cho_cov_(chol(cov,"upper")),
logDetCov_(log(diagProduct(cho_cov_)))
{}

MultivariateGaussian::MultivariateGaussian():
    mean_(),
    covinv_(),
    cho_cov_(),
    logDetCov_()
{}

double MultivariateGaussian::logP(const M_Matrix<double> &x) const
{

     return -0.5*size()*log(PI)-logDetCov_-0.5*xTSigmaX(x-mean_,covinv_);
}


MultivariateGaussian::~MultivariateGaussian(){}


const M_Matrix<double>& MultivariateGaussian::Mean()const
{
    return mean_;
}
M_Matrix<double> MultivariateGaussian::Cov()const
{
    return inv(covinv_);
}

const M_Matrix<double>& MultivariateGaussian::CovInv()const
{
    return covinv_;
}



