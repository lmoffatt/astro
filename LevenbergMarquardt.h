#ifndef LEVENBERGMARQUARDT
#define LEVENBERGMARQUARDT
#include "Parameters.h"


class ABC_function
{
public:
    virtual std::vector<double> yfit(const std::vector<double>& parameters)=0;
    virtual std::vector<double> yfit(const Parameters& parameters)const=0;
    virtual std::ostream& put(std::ostream& s,const Parameters& parameters)const=0;

};


class LevenbergMarquardtParameters
{
public:

    Parameters OptimParameters()const;

    std::size_t numEval()const;
    std::size_t numIter()const;
    double SS()const;
    std::vector<double> Gradient()const;


    LevenbergMarquardtParameters& optimize();


    LevenbergMarquardtParameters(ABC_function* f,
                       const std::vector<double>& data,
                       const Parameters& initialParam,
                       const std::vector<double>& weigth,
                                 std::size_t numIterations);




    double getEvidence()const;

    double getLogPostLik()const;

    double logDetPriorCov()const;
    double logDetPostCov()const;
    double logDetPostStd()const;
    double SSdata()const;


    LevenbergMarquardtParameters(const LevenbergMarquardtParameters& other);

    friend void swap(LevenbergMarquardtParameters& one, LevenbergMarquardtParameters& other);

    LevenbergMarquardtParameters& operator=(const LevenbergMarquardtParameters& other);

    LevenbergMarquardtParameters();

    ~LevenbergMarquardtParameters(){}
    std::string report()const;

   // void reset(const SimParameters& sp,const Treatment& tr);


    friend std::ostream& operator<<(std::ostream& s, LevenbergMarquardtParameters& LM);

private:
    void calculateCovariance();
    ABC_function* f_;
    std::vector<double> data_;

    std::vector<double> w_;

    Parameters initialParam_;

    std::size_t nPar_;
    std::size_t nData_;



    // parameters of the optimization
    /// delta x used for Jacobian approximation
    double dx_;
    std::size_t maxIter_;
    std::size_t maxFeval_;

    double minParamChange_;
    double minSSChange_;
    double minGradient_;

    double maxLanda_;

    // variables that change on each iteration

    double landa_;
    double landa0_;
    double v_;

    std::size_t nIter_;
    std::size_t nFeval_;



    double currSS_;
    double newSSW_;
    double newSSW0_;
    Parameters currParam_;
    Parameters newParam_;
    Parameters newParam0_;

    std::vector<double> currYfit_;
    std::vector<double> newYfit_;
    std::vector<double> newYfit0_;

    std::vector< std::vector< double> > J_;
    std::vector<double> G_;
    std::vector< std::vector<double> > JTWJ_;
    std::vector< std::vector<double> > JTWJinv_;

    double evidence_;
    double LogPostLik_;
    double logDetPriorCov_;
    double logDetPostCov_;
    double logDetPostStd_;
    double SSdata_;


    std::vector<double> d_;


    Parameters optimParam_;

    bool surpassIter_;
    bool surpassFeval_;
    bool surpassLanda_;

    double ParamChange_;
    double SSChange_;
    double NormGrad_;

    bool smallParamChange_;

    bool smallSSChange_;

    bool smallGradient_;


    bool meetConvergenceCriteria();

    void initialize();
    void iterate();

    void computeJacobian();
    void computeSearchDirection();
    void updateLanda();


    double SumSquare();




};


#endif // LEVENBERGMARQUARDT

