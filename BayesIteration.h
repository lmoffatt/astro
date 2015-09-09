#ifndef BAYESITERATION
#define BAYESITERATION


#include <vector>
#include <string>
#include <map>
#include "LevenbergMarquardt.h"





class BayesIteration:public ABC_Multinomial_Model, public ABC_Freq_obs
{
public:

    Parameters Posterior()const;

    Parameters Prior(std::size_t n_obs=0)const;


    BayesIteration(const CortexLikelihood *f,
                   Parameters prior,
                   const ABC_Freq_obs* d,
                   const std::string &filename);


    std::map<double,Parameters> getRandomParameters(std::mt19937& mt,std::size_t num,double factor);

    std::map<double,Parameters> getRandomParameters(std::mt19937& mt,const Parameters& per,std::size_t num, double factor);


    double SumWeighedSquare(const Parameters& p);



    BayesIteration& addNewData(ABC_Freq_obs* d);

    virtual ~BayesIteration();

    virtual void setFilename(const std::string filename);

    virtual std::vector<std::vector<double>> f (const Parameters& parameters)const;

    virtual const ABC_Freq_obs& getData()const;

     BayesIteration(const BayesIteration& other);
    /*

    friend void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);
*/
    BayesIteration();

    BayesIteration& getPosterior();

    BayesIteration& getPosterior(const Parameters& startingPoint);

    BayesIteration& getPosterior(const Parameters& startingPoint, double factor, std::size_t numSeeds);

    Parameters getEvidence(const Parameters& maximumPostLik, std::size_t num);

    Parameters getHessian(const Parameters& MAP, double eps=1e-3);


    Parameters getHessianInterpol(const Parameters& MAP, double minP, double maxP);

    virtual std::ostream& put(std::ostream& s,const Parameters& parameters)const;



private:

    const CortexLikelihood* m_;

    std::vector<const ABC_Freq_obs*> data_;

    std::vector<Parameters> priors_;

    Parameters posterior_;
    std::size_t numSeeds_;

    std::vector<LevenbergMarquardtDistribution> LM_;

    std::string filename_;


    };


#endif // BAYESITERATION

