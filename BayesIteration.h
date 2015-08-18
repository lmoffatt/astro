#ifndef BAYESITERATION
#define BAYESITERATION


#include <vector>
#include <string>
#include <map>
#include "LevenbergMarquardt.h"


class ABC_data
{
public:
    virtual std::vector<double> getData()const=0;
    virtual std::vector<double> getDataStandardError()const=0;
    virtual std::vector<double> getDataWeigth();
    virtual ~ABC_data();
};


class ABC_model:public ABC_function
{
public:

    virtual ABC_model& setData(const ABC_data& experimentalData)=0;
    virtual ~ABC_model();

};




class BayesIteration:public ABC_model, public ABC_data
{
public:

    Parameters Posterior()const;

    Parameters Prior(std::size_t n=0)const;


    BayesIteration(const ABC_model* f,
                   Parameters prior,
                   const ABC_data* d,
                   const std::string &filename);


    std::map<double,Parameters> getRandomParameters(std::size_t num,double factor);

    std::map<double,Parameters> getRandomParameters(const Parameters& per,std::size_t num, double factor);


    double SumWeighedSquare(const Parameters& p);



    BayesIteration& addNewData(ABC_data* d);

    virtual ~BayesIteration();

    virtual void setFilename(const std::string filename);

    virtual std::vector<double> yfit (const std::vector<double>& param);
    virtual std::vector<double> yfit(const Parameters& parameters)const;

    virtual std::vector<double> getData()const;
    virtual std::vector<double> getDataStandardError()const;

    virtual ABC_model& setData(const ABC_data& experimentalData);

    BayesIteration(const BayesIteration& other);
    /*

    friend void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);
*/
    BayesIteration();

    BayesIteration& getPosterior();

    BayesIteration& getPosterior(const Parameters& startingPoint);

    BayesIteration& getPosterior(const Parameters& startingPoint,double factor, std::size_t numSeeds,double probParChange);

    Parameters getEvidence(const Parameters& maximumPostLik, std::size_t num);

    Parameters getHessian(const Parameters& MAP, double eps=1e-3);


    Parameters getHessianInterpol(const Parameters& MAP, double minP, double maxP);

    virtual std::ostream& put(std::ostream& s,const Parameters& parameters)const;



private:

    const ABC_model* m_;

    std::vector<const ABC_data*> data_;

    std::vector<Parameters> priors_;

    Parameters posterior_;
    std::size_t numSeeds_;

    std::vector<LevenbergMarquardtParameters> LM_;

    std::string filename_;


    };


#endif // BAYESITERATION

