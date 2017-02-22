#ifndef CORTEXLIKELIHOOD
#define CORTEXLIKELIHOOD
#include "Models.h"
#include "Evidence.h"


class CortexLikelihood: public BaseAgent,virtual public ABC_Distribution_Model, public ABC_Freq_obs
{
  // BaseClass interface
public:
  static std::string ClassName();
  virtual std::string myClass() const override;

  // BaseObject interface
public:
  virtual CortexLikelihood* create() const override=0;
  virtual CortexSimulation simulate(const Parameters &parameters) const=0;
  virtual std::vector<std::vector<double> > g(const Parameters &parameters, const CortexSimulation &s) const=0;
  CortexLikelihood();
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override;
  virtual bool readBody(std::string &line, std::istream &s) override;

  std::vector<double> getNBins(const Experiment *e);
  // ABC_Multinomial_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override;

  virtual const ABC_Freq_obs &getData() const override;
  virtual const Parameters &getPrior() const override;

  CortexLikelihood(std::string id,
                   const Experiment* e,
                   const Parameters& prior
                   , double dx
                   , double dtmin,
                   std::size_t nPoints_per_decade,
                   double dtmax
                   , double tequilibrio);

  // ABC_Freq_obs interface
  virtual const std::vector<std::vector<double>> & n_obs()const override;
  virtual const std::vector<double> &n_obs(unsigned i) const override;
  virtual const std::vector<double>& bin_dens()const override;

  virtual const std::vector<double> &ntot_obs() const override;
  virtual std::ostream &put(std::ostream &s) const override;

  ~CortexLikelihood();

  std::ostream& write(std::ostream& s)const;


  CortexLikelihood(const CortexLikelihood& o)=delete;

  const Experiment* getExperiment()const;

  const BaseModel* getModel()const;

protected:
  void update();

  std::vector<std::vector<double>> getstate(const Experiment* e);

  std::vector<double> getNtot(const std::vector<std::vector<double>> nstate);

  Parameters prior_;
  const BaseModel* m_;
  const Experiment* e_;
  double dx_;
  double dtmin_;
  std::size_t nPoints_per_decade_;
  double dtmax_;
  double tequilibrio_;
  std::vector<std::vector<double> > nstate_;
  std::vector<double> ntot_;
  std::vector<double> bin_dens_;


};


class CortexMultinomialLikelihood:public ABC_Multinomial_Model, virtual public CortexLikelihood
{


  // BaseClass interface
public:

  CortexMultinomialLikelihood(std::string id,
                              const Experiment* e,
                              const Parameters& prior
                              ,double dx
                              ,double dtmin
                              ,std::size_t nPoints_per_decade
                              ,double dtmax
                              ,double tequilibrio):
    CortexLikelihood(id,e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio){}

  static std::string ClassName(){return "CortexMultinomialLikelihood";}
  virtual std::string myClass() const override{return ClassName();}

  // ABC_Distribution_Model interface
public:
  virtual std::vector<std::vector<double> > f(const Parameters &parameters) const override;

  // BaseObject interface
public:
  CortexMultinomialLikelihood(){}
  virtual CortexMultinomialLikelihood *create() const override
  {
    return new CortexMultinomialLikelihood;
  }

  // BaseObject interface

  // ABC_Distribution_Model interface
public:
  virtual std::vector<std::vector<double> > logLikCells(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::logLikCells(p);
  }
  virtual std::vector<double> logLikSamples(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::logLikSamples(p);
  }
  virtual double logLik(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::logLik(p);
  }
  virtual std::vector<double> epsilon(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::epsilon(p);
  }
  virtual const std::vector<double> weight(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::weight(p);
  }

  // ABC_Multinomial_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override
  {
    CortexLikelihood::setPrior(parameters);
  }
  virtual const ABC_Freq_obs &getData() const override
  {
    return CortexLikelihood::getData();
  }
  virtual const Parameters &getPrior() const override
  {
    return CortexLikelihood::getPrior();
  }

  // CortexLikelihood interface
public:
  virtual CortexSimulation simulate(const Parameters &parameters) const override
  {
  }
  virtual std::vector<std::vector<double> > g(const Parameters &parameters, const CortexSimulation &s) const override
  {
  }
};


class CortexPoisonLikelihood: public ABC_MultiPoison_Model, public CortexLikelihood
{
  // ABC_Multinomial_Model interface
public:


  // ABC_Multinomial_Model interface
  std::vector<std::vector<double> > f(const Parameters &parameters) const override;




  CortexPoisonLikelihood(std::string id,
                         const Experiment* e,
                         const Parameters& prior
                         ,double dx
                         ,double dtmin
                         ,std::size_t nPoints_per_decade
                         ,double dtmax
                         ,double tequilibrio):
    CortexLikelihood(id,e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio){}

  // ABC_Freq_obs interface

  ~CortexPoisonLikelihood(){}

public:
  static std::string ClassName(){return "CortexPoisonLikelihood";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual CortexPoisonLikelihood* create() const override
  {
    return new CortexPoisonLikelihood;
  }
  CortexPoisonLikelihood(){}


  // ABC_Distribution_Model interface
public:
  virtual std::vector<std::vector<double> > logLikCells(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::logLikCells(p);
  }
  virtual std::vector<double> logLikSamples(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::logLikSamples(p);
  }
  virtual double logLik(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::logLik(p);
  }
  virtual std::vector<double> epsilon(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::epsilon(p);
  }
  virtual const std::vector<double> weight(const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::weight(p);
  }


  // ABC_Distribution_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override
  {
    CortexLikelihood::setPrior(parameters);
  }
  virtual const ABC_Freq_obs &getData() const override
  {
    return CortexLikelihood::getData();
  }
  virtual const Parameters &getPrior() const override
  {
    return CortexLikelihood::getPrior();
  }
  CortexSimulation simulate(const Parameters &parameters) const;
  std::vector<std::vector<double> > g(const Parameters &parameters, const CortexSimulation &s) const;
};




class CortexMultinomialLikelihoodEvaluation: public BaseAgent
{
  // BaseClass interface
public:
  static std::string ClassName(){return "CortexLikelihoodEvaluation";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  CortexMultinomialLikelihoodEvaluation(){}
  ~CortexMultinomialLikelihoodEvaluation(){}
  virtual CortexMultinomialLikelihoodEvaluation *create() const override
  {
    return new CortexMultinomialLikelihoodEvaluation;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override{}
  virtual std::ostream &extract(std::ostream &s, const std::string & ="",
                                const std::string& ="") const override;



  virtual bool readBody(std::string &, std::istream &) override{return false;}

  CortexMultinomialLikelihoodEvaluation(const CortexLikelihood& CL,
                                        const Parameters& p)
    :
      CL_(&CL)
    ,p_(p)
  {}





protected:
  virtual void update() override{}
  const CortexLikelihood* CL_;
  Parameters p_;

  //--------Calculated variables-----------------//

  std::vector<std::vector<double> > p_exp_;


};





class CortexPoisonLikelihoodEvaluation: public BaseAgent
{
  // BaseClass interface
public:
  static std::string ClassName(){return "CortexPoisonLikelihoodEvaluation";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  CortexPoisonLikelihoodEvaluation(){}
  ~CortexPoisonLikelihoodEvaluation(){}
  virtual CortexPoisonLikelihoodEvaluation *create() const override
  {
    return new CortexPoisonLikelihoodEvaluation;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override{}
  virtual std::ostream &extract(std::ostream &s, const std::string & ="",
                                const std::string& ="") const override;



  virtual bool readBody(std::string &, std::istream &) override{return false;}

  CortexPoisonLikelihoodEvaluation(const CortexPoisonLikelihood& CL,
                                   const Parameters& p)
    :
      CL_(&CL)
    ,p_(p)
  {}





protected:
  virtual void update() override{}
  const CortexPoisonLikelihood* CL_;
  Parameters p_;

  //--------Calculated variables-----------------//

  std::vector<std::vector<double> > landa_;


};

template<typename T>
using myMatrix= std::vector<std::vector<T>> ;




class MyData
{
public:
  M_Matrix<std::size_t> operator()()const
  {

    std::size_t nrows= e_->n_obs().size();
    if (nrows>0)
      {
        std::size_t ncols=e_->n_obs()[0].size();
        M_Matrix<std::size_t> out(nrows,ncols);
        for (std::size_t i=0; i<nrows; ++i)
          for (std::size_t j=0; j<ncols; ++j)
            out(i,j)=e_->n_obs()[i][j];
        return out;
      }
    else return {};
  }

  MyData(ABC_Freq_obs *e ):e_(e){}

private:
  const ABC_Freq_obs *e_;
};

template<class Data=MyData>
class MyModel
{
public:
  M_Matrix<double> sample(std::mt19937_64& mt)const
  {
    Parameters sample=CL_->getPrior().randomSample(mt,1);
    return M_Matrix<double> (1,sample.size(),sample.trMeans());



  }
  mcmc_prior prior(const Data& data, M_Matrix<double> param)const
  {
    mcmc_prior out;
    out.param=param;

    auto& p=CL_->getPrior();
    std::size_t npar=param.size();

    Parameters Pa=p.toParameters(param.toVector());

    out.logPrior=p.logProb(Pa);
    out.D_prior.H=M_Matrix<double>(p.getInvCovariance());
    M_Matrix<double> d(param.nrows(),param.ncols(),p.trMeans());
    d=param - d;


    out.D_prior.G=d*out.D_prior.H;

    return out;

  }
  M_Matrix<double> f(const Data& data, M_Matrix<double> param)const{
    auto ff=CL_->f(param.toVector());
    std::size_t nrows= ff.size();
    if (nrows>0)
      {
        std::size_t ncols=ff[0].size();
        M_Matrix<double> out(nrows,ncols);
        for (std::size_t i=0; i<nrows; ++i)
          for (std::size_t j=0; j<ncols; ++j)
            out(i,j)=ff[i][j];
        return out;
      }
    else return {};

  }
  M_Matrix<double> logLanda(const Data& data, M_Matrix<double> param)const{
    M_Matrix<double> out=f(data,param);
    out=out.apply([](double x){return log10_guard(x);});
    return out;

  }

  MyModel(CortexLikelihood* CL):CL_(CL){}
 private:
   const CortexLikelihood* CL_;
};

inline
std::vector<std::tuple<double,std::size_t,std::size_t>>
getBeta(const M_Matrix<double>& b
        ,const M_Matrix<std::size_t>& samples
        , const M_Matrix<std::size_t>& skip)
{
  std::vector<std::tuple<double,std::size_t,std::size_t>>
 out(b.size());

  for (std::size_t i=0; i<b.size(); ++i)
     {
       std::get<0>(out[i])=b[i];
       if (samples.size()>1)
       std::get<1>(out[i])=samples[i];
       else
         std::get<1>(out[i])=samples[0];
      if (skip.size()>1)
       std::get<2>(out[i])=skip[i];
      else
        std::get<2>(out[i])=skip[0];
     }
  return out;

}





typedef
Thermodynamic_Integration_mcmc<
MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa,LevenbergMarquardt_step> TI;

typedef
Template_Tempering_mcmc<
MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa,LevenbergMarquardt_step> TT;





















#endif // CORTEXLIKELIHOOD

