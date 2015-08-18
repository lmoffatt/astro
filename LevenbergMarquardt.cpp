#include "LevenbergMarquardt.h"
#include "MatrixInverse.h"

#include <vector>
#include <limits>
#include <cmath>
#include <sstream>


LevenbergMarquardtParameters::LevenbergMarquardtParameters(
    ABC_function* f,
    const std::vector<double>& data,
    const Parameters& initialParam,
    const std::vector<double>& weigth, std::size_t numIterations):
    f_(f),
    data_(data),
    w_(weigth),
    initialParam_(initialParam),
    nPar_(initialParam_.size()),
    nData_(data.size()),
    dx_(1e-9),
    maxIter_(numIterations),
    maxFeval_(10000),
    minParamChange_(1e-9),
    minSSChange_(1e-9),
    minGradient_(0.0000001),
    maxLanda_(1e11),
    landa_(1000),
    v_(10),
    nIter_(0),
    nFeval_(0),
    currSS_(std::numeric_limits<double>::infinity()),
    newSSW0_(std::numeric_limits<double>::infinity()),
    currParam_(initialParam_),
    newParam_(currParam_),
    newParam0_(currParam_),
    currYfit_(std::vector<double>(nData_)),
    newYfit_(currYfit_),
    newYfit0_(currYfit_),
    J_(std::vector< std::vector< double> >(nData_,std::vector<double>(nPar_))),
    G_(currParam_.size()),
    JTWJ_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
    JTWJinv_(JTWJ_),
    d_(currParam_.size()),
    optimParam_(currParam_),
    surpassIter_(false),
    surpassFeval_(false),
    ParamChange_(std::numeric_limits<double>::infinity()),
    SSChange_(ParamChange_),
    NormGrad_(ParamChange_),
    smallParamChange_(false),
    smallSSChange_(false),
    smallGradient_(false){}



LevenbergMarquardtParameters::LevenbergMarquardtParameters(){}

LevenbergMarquardtParameters::LevenbergMarquardtParameters (const LevenbergMarquardtParameters& other):
    data_(other.data_),
    w_(other.w_),
    initialParam_(other.initialParam_),
    nPar_(other.nPar_),
    nData_(other.nData_),
    dx_(other.dx_),
    maxIter_(other.maxIter_),
    maxFeval_(other.maxFeval_),
    minParamChange_(other.minParamChange_),
    minSSChange_(other.minSSChange_),
    minGradient_(other.minGradient_),
    maxLanda_(other.maxLanda_),
    landa_(other.landa_),
    landa0_(other.landa0_),
    v_(other.v_),
    nIter_(other.nIter_),
    nFeval_(other.nFeval_),



    currSS_(other.currSS_),
    newSSW0_ (other.newSSW0_),
    currParam_ (other.currParam_),
    newParam_ (other.newParam_),
    newParam0_ (other.newParam0_),
            currYfit_(other.currYfit_),
            newYfit_(other.newYfit_),
    newYfit0_(other.newYfit0_),

            J_(other.J_),
            G_(other.G_),
            JTWJ_(other.JTWJ_),
            JTWJinv_(other.JTWJinv_),

            d_(other.d_),

            optimParam_(other.optimParam_),

            surpassIter_(other.surpassIter_),
            surpassFeval_(other.surpassFeval_),

            ParamChange_(other.ParamChange_),
            SSChange_(other.SSChange_),
            NormGrad_(other.NormGrad_),

            smallParamChange_(other.smallParamChange_),

            smallSSChange_(other.smallSSChange_),

            smallGradient_(other.smallGradient_)

    {}



LevenbergMarquardtParameters&
LevenbergMarquardtParameters::operator=(const LevenbergMarquardtParameters& other)
{
    if (this!=&other)
    {
        LevenbergMarquardtParameters tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(LevenbergMarquardtParameters& one, LevenbergMarquardtParameters& other)
{
    std::swap(one.data_,other.data_);
    std::swap(one.w_,other.w_);

    std::swap(one.initialParam_,other.initialParam_);
    std::swap(one.nPar_,other.nPar_);
    std::swap(one.nData_,other.nData_);
    std::swap(one.dx_,other.dx_);
    std::swap(one.maxIter_,other.maxIter_);
    std::swap(one.maxFeval_,other.maxFeval_);
    std::swap(one.minParamChange_,other.minParamChange_);
    std::swap(one.minSSChange_,other.minSSChange_);
    std::swap(one.minGradient_,other.minGradient_);
    std::swap(one.maxLanda_,other.maxLanda_);
    std::swap(one.landa_,other.landa_);
    std::swap(one.landa0_,other.landa0_);
    std::swap(one.v_,other.v_);
    std::swap(one.nIter_,other.nIter_);
    std::swap(one.nFeval_,other.nFeval_);



    std::swap(one.currSS_,other.currSS_);
    std::swap(one.newSSW0_,other.newSSW0_);
    std::swap(one.currParam_,other.currParam_);
    std::swap(one.newParam_,other.newParam_);
    std::swap(one.newParam0_,other.newParam0_);
    std::swap(one.currYfit_,other.currYfit_);
    std::swap(one.newYfit_,other.newYfit_);
    std::swap(one.newYfit_,other.newYfit_);
std::swap(one.J_,other.J_);
    std::swap(one.G_,other.G_);
    std::swap(one.JTWJ_,other.JTWJ_);
    std::swap(one.JTWJinv_,other.JTWJinv_);

    std::swap(one.d_,other.d_);

    std::swap(one.optimParam_,other.optimParam_);

    std::swap(one.surpassIter_,other.surpassIter_);
    std::swap(one.surpassFeval_,other.surpassFeval_);

    std::swap(one.ParamChange_,other.ParamChange_);
    std::swap(one.SSChange_,other.SSChange_);
    std::swap(one.NormGrad_,other.NormGrad_);

    std::swap(one.smallParamChange_,other.smallParamChange_);

    std::swap(one.smallSSChange_,other.smallSSChange_);

    std::swap(one.smallGradient_,other.smallGradient_);

    }

LevenbergMarquardtParameters& LevenbergMarquardtParameters::optimize()
{
    std::cout<<"nIter"<<"\t"<<"currSS"<<"\t\t"<<"landa"<<"\t";
    std::cout<<"ParamChange"<<"\t"<<"SSChange"<<"\t"<<"NormGrad"<<"\n";
    initialize();
    while (!meetConvergenceCriteria())
            iterate();

    std::cout<<report();

    optimParam_=currParam_;
    for (std::size_t i=0; i<nPar_; ++i)
        for (std::size_t j=0; j<nPar_; ++j)
        {
            JTWJ_[i][j]=0;
            for (std::size_t n=0; n<nData_; ++n)
            {
                JTWJ_[i][j]+=J_[n][i]*J_[n][j]*w_[n];
            }
        }

    JTWJinv_=inv(JTWJ_);
    evidence_=0;
    logDetPriorCov_=0;
    LogPostLik_=0;
    logDetPostStd_=0;
    for (std::size_t i=0; i<w_.size();i++)
    {
        LogPostLik_+=0.5*log(w_[i]);
    }
    for (std::size_t i=currParam_.size(); i<w_.size();i++)
    {
        logDetPriorCov_+=log(w_[i]);
    }
    for (std::size_t i=0;i<currParam_.size(); i++)
    {
        logDetPostStd_+=log(JTWJinv_[i][i]);
    }


    LogPostLik_+= -0.5*currSS_;
    logDetPostCov_=log(det(JTWJinv_));
    evidence_=LogPostLik_+0.5*log(det(JTWJinv_));





    optimParam_.setCovariance(JTWJinv_);
    return *this;
}

double LevenbergMarquardtParameters::getEvidence()const
{

return evidence_;
}
double LevenbergMarquardtParameters::getLogPostLik()const
{

return LogPostLik_;
}

double LevenbergMarquardtParameters::logDetPriorCov()const
{
    return logDetPriorCov_;
}
double LevenbergMarquardtParameters::logDetPostCov()const
{
    return logDetPostCov_;
}

double LevenbergMarquardtParameters::logDetPostStd()const
{
    return logDetPostStd_;
}

double LevenbergMarquardtParameters::SSdata()const
{
    return SSdata_;
}


void LevenbergMarquardtParameters::calculateCovariance()
{

}

void LevenbergMarquardtParameters::iterate()
{
    computeJacobian();
    computeSearchDirection();
    updateLanda();
    std::cout<<nIter_<<"\t"<<currSS_<<"\t\t"<<landa_<<"\t";
    std::cout<<ParamChange_<<"\t"<<SSChange_<<"\t"<<NormGrad_<<"\n";

    nIter_++;
}


void LevenbergMarquardtParameters::computeJacobian()
{

    for (std::size_t i=0; i<nPar_; i++)
    {
        Parameters x(currParam_);
        x[i]+=dx_;
        std::vector<double> yfit=f_->yfit(x);
        //std::cout<<nFeval_++<<" ";
        for (std::size_t n=0;n<nData_;++n)
        {
            double jj=(yfit[n]-currYfit_[n])/dx_;
            J_[n][i]=jj;
        }
    }
}


void LevenbergMarquardtParameters::computeSearchDirection()
{
    for (std::size_t i=0; i<nPar_; ++i)
        for (std::size_t j=0; j<nPar_; ++j)
        {
            JTWJ_[i][j]=0;
            for (std::size_t n=0; n<nData_; ++n)
            {
                JTWJ_[i][j]+=J_[n][i]*J_[n][j]*w_[n];
            }
        }



    for (std::size_t i=0; i<nPar_; ++i)
    {
        JTWJ_[i][i]*=1+landa_;
    }

    JTWJinv_=inv(JTWJ_);



    for (std::size_t i=0; i<nPar_; ++i)
    {
        G_[i]=0;
        for (std::size_t n=0; n<nData_;++n)
        {
            G_[i]+=(data_[n]-currYfit_[n])*w_[n]*J_[n][i];
        }

    }
    for (std::size_t i=0; i<nPar_; ++i)
    {
        d_[i]=0;
        for (std::size_t j=0; j<nPar_;++j)
        {
            d_[i]+=JTWJinv_[i][j]*G_[j];
        }

    }
    newParam_=currParam_;
    for (std::size_t i=0; i<nPar_; ++i)
        newParam_[i]+=d_[i];

    newYfit_=f_->yfit(newParam_);
    nFeval_++;

    newSSW_=0;
    for (std::size_t n=0; n<nData_; n++)
    {
        newSSW_+=(newYfit_[n]-data_[n])*(newYfit_[n]-data_[n])*w_[n];
    }
}


void LevenbergMarquardtParameters::updateLanda()
{
    std::size_t ifevalLoop=0;
    if ((newSSW_>=currSS_)||(newSSW_!=newSSW_))
    {
        while(((newSSW_>=currSS_)&&(nFeval_<maxFeval_))||(newSSW_!=newSSW_))
        {
            if (landa_*v_>=maxLanda_) break;
            landa0_=landa_;
            landa_=landa0_*v_;
            newSSW0_=newSSW_;
            newParam0_=newParam_;
            computeSearchDirection();
            ifevalLoop++;
         //   std::cerr<<landa_<<" ";
        }

    }
    else
    {
        landa0_=landa_;
        landa_=landa_/v_;
        newSSW0_=newSSW_;
        newParam0_=newParam_;
        newYfit0_=newYfit_;
        computeSearchDirection();
        ifevalLoop++;
        if ((newSSW_>=newSSW0_)||(newSSW_!=newSSW_))
        {
            landa_=landa0_;
            newParam_=newParam0_;
            newSSW_=newSSW0_;
            newYfit_=newYfit0_;
        }
    }
    if (newSSW_>=currSS_)
    {
        ParamChange_=0;
        SSChange_=0;

    }
    else
    {
        ParamChange_=0;
        for (std::size_t i=0; i<nPar_; ++i)
            ParamChange_+=(currParam_[i]-newParam_[i])*(currParam_[i]-newParam_[i]);
        ParamChange_=sqrt(ParamChange_);
        SSChange_=currSS_-newSSW_;
        currParam_=newParam_;
        currSS_=newSSW_;
        currYfit_=newYfit_;
    }
    NormGrad_=0;
    for (std::size_t i=0; i<nPar_; ++i)
        NormGrad_+=G_[i]*G_[i];
    NormGrad_=sqrt(NormGrad_);
 }


bool LevenbergMarquardtParameters::meetConvergenceCriteria()
{
    surpassIter_=(bool(nIter_>=maxIter_)&&(currSS_>100))||(nIter_>=maxIter_*5);
    surpassFeval_=nFeval_>=maxFeval_;
    surpassLanda_=landa_>=maxLanda_;
    smallParamChange_=ParamChange_<minParamChange_;
    smallSSChange_=SSChange_<minSSChange_;
    smallGradient_=NormGrad_<minGradient_;


    return surpassIter_||
            surpassFeval_||
            smallParamChange_||
            smallSSChange_||
            smallGradient_||
            surpassLanda_||
            currSS_!=currSS_;
}

void LevenbergMarquardtParameters::initialize()
{
    nIter_=0;
    nFeval_=0;
    ParamChange_=std::numeric_limits<double>::infinity();
    SSChange_=std::numeric_limits<double>::infinity();
    NormGrad_=std::numeric_limits<double>::infinity();
    currParam_=initialParam_;

    currYfit_=f_->yfit(currParam_);
    nFeval_++;
    currSS_=0;
    SSdata_=0;
    std::vector<double> ss(nData_);
    for (std::size_t n=0; n<nData_; ++n)
    {
       ss[n]=(currYfit_[n]-data_[n])*(currYfit_[n]-data_[n])*w_[n];
        currSS_+=ss[n];
        if (n<nData_-currParam_.size())
            SSdata_+=ss[n];
    }
    nFeval_++;
    nFeval_--;




}


Parameters LevenbergMarquardtParameters::OptimParameters()const
{
    return optimParam_;
}

std::size_t LevenbergMarquardtParameters::numEval()const
{
    return nFeval_;
}
std::size_t LevenbergMarquardtParameters::numIter()const
{
    return nIter_;
}
double LevenbergMarquardtParameters::SS()const
{
    return currSS_;
}
std::vector<double> LevenbergMarquardtParameters::Gradient()const
{
    return G_;
}
std::string LevenbergMarquardtParameters::report()const{
    std::stringstream output;

    output<<"Convergence critera: \t";
    if (surpassIter_)
        output<<nIter_<<" surpass number of iterations="<<maxIter_<<"\t";
    if (surpassFeval_)
        output<<nFeval_<<" surpass number of function evaluations="<<maxFeval_<<"\t";
    if (surpassLanda_)
        output<<landa_<<" surpass maximum landa value="<<maxLanda_<<"\t";
    if (smallParamChange_)
        output<<ParamChange_<<" surpass minium parameter change value="<<minParamChange_<<"\t";
    if (smallSSChange_)
        output<<SSChange_<<" surpass minium SS change value="<<minSSChange_<<"\t";
    if (smallGradient_)
        output<<NormGrad_<<" surpass minium gradient norm value="<<minGradient_<<"\t";
    if (currSS_!=currSS_)
        output<<currSS_<<" invalid value of the square sum"<<"\t";




    return output.str();



}


 std::ostream& operator<<(std::ostream& s, LevenbergMarquardtParameters& LM)
{
    s<<"LevenbergMarquardtParameters\n";
    s<<"Number of iterations \t"<<LM.numIter()<<"\t";
    s<<"Sum of squares"<<LM.currSS_<<"\n";
    s<<"Report \t"<<LM.report()<<"\n";
    s<<"Landa value \t"<<LM.landa_<<"\t";
    s<<"Parameter change \t"<<LM.minParamChange_<<"\t";
    s<<"SS change \t"<<LM.minSSChange_<<"\t";
    s<<"Gradient norm \t"<<LM.NormGrad_<<"\n";
    s<<"Evidence \t"<<LM.evidence_<<"\n";
    s<<"End\n";

}



