#ifndef CORTEXLIKELIHOOD
#define CORTEXLIKELIHOOD
#include "Models.h"
//#include "Evidence.h"


class CortexLikelihood: public BaseAgent,virtual public ABC_Distribution_Model, public ABC_Freq_obs
{
  // BaseClass interface
public:
  static std::string ClassName();
  virtual std::string myClass() const override;

  // BaseObject interface
public:
  virtual CortexLikelihood* create() const override=0;
  virtual CortexSimulation simulate(const Experiment* e,const Parameters &parameters,const std::pair<std::vector<double>, std::vector<std::size_t>>& dts) const=0;
  virtual std::vector<std::vector<double> > g(const Experiment* e,const Parameters &parameters, const CortexSimulation &s) const=0;
  CortexLikelihood();
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override;
  virtual bool readBody(std::string &line, std::istream &s, std::ostream& logs) override;

  virtual std::ostream& writeYfitHeaderDataFrame(const Experiment*,std::ostream& os,const CortexSimulation &)const{return os;}


  std::vector<double> getNBins(const Experiment *e) const;
  // ABC_Multinomial_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override;

  virtual const ABC_Freq_obs &getData() const override;
  virtual const Parameters &getPrior() const override;

  CortexLikelihood(std::string id,
                   const Parameters& prior
                   , double dx
                   , double dtmin0,
                   double dtmin,
                   std::size_t nPoints_per_decade,
                   double dtmax
                   , double tequilibrio);

  CortexLikelihood(std::string id,
                   const Parameters& prior
                   , double dx
                   , double dtmin0,
                   double dtmin,
                   std::size_t nPoints_per_decade,
                   double dtmax
                   , double tequilibrio
                   , double maxlogError
                   , double dtinf);

  // ABC_Freq_obs interface
  virtual const std::vector<std::vector<double>>  n_obs(const Experiment* e)const override;
  virtual const std::vector<double> bin_dens(const Experiment* e)const override;

  virtual const std::vector<double> ntot_obs(const Experiment* e) const override;
  virtual std::ostream &put(std::ostream &s) const override;

  ~CortexLikelihood();

  std::ostream& write(std::ostream& s)const;


  CortexLikelihood(const CortexLikelihood& o)=delete;


  const BaseModel* getModel()const;




  CortexLikelihood(std::string id,
                   const Parameters &prior,
                   double dx,
                   double dtmin0,
                   double dtmin,
                   std::size_t nPoints_per_decade,
                   double dtmax,
                   double tequilibrio,
                   double maxlogError,
                   double f_maxlogError,
                   double dtinf,
                   std::size_t maxloop
                   , bool UseDerivative);

  CortexLikelihood(std::string id,
                   const Parameters &prior,
                   double dx,
                   double dtmin0,
                   double dtmin,
                   std::size_t nPoints_per_decade,
                   double dtmax,
                   double tequilibrio,
                   double t_maxlogError,
                   std::size_t maxloop);
protected:
  void update() override;

  std::vector<std::vector<double>> getstate(const Experiment* e) const;

  std::vector<double> getNtot(const std::vector<std::vector<double>> nstate) const;

  Parameters prior_;
  const BaseModel* m_;
  bool adapt_dt_;
  bool CrankNicholson_;
  double dx_;
  double dtmin0_;
  double dtmin_;
  std::size_t nPoints_per_decade_;
  double dtmax_;
  double tequilibrio_;
  double dtinf_;
  double maxlogError_;
  double f_maxlogErrorCN_;
  std::size_t maxloop_;
  bool UseDerivative_;


};


class CortexMultinomialLikelihood:public ABC_Multinomial_Model, virtual public CortexLikelihood
{


  // BaseClass interface
public:

  CortexMultinomialLikelihood(std::string id,
                              const Parameters& prior
                              ,double dx
                              ,double dtmin0
                              ,double dtmin
                              ,std::size_t nPoints_per_decade
                              ,double dtmax
                              ,double tequilibrio):
    CortexLikelihood(id,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio){}

  static std::string ClassName(){return "CortexMultinomialLikelihood";}
  virtual std::string myClass() const override{return ClassName();}

  // ABC_Distribution_Model interface
public:
  virtual std::vector<std::vector<double> > f(const Experiment *e, const Parameters &parameters, std::pair<std::vector<double>,std::vector<std::size_t>>& dts) const override;


  virtual std::vector<std::vector<double> > f
  (const Experiment* e,const Parameters &parameters,
   std::tuple<double,std::size_t,double>& dt_n_dtmax,
   std::size_t dt_max,
   std::pair<std::vector<double>, std::vector<std::size_t> > &) const override;

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
  virtual std::vector<std::vector<double> > logLikCells(const Experiment* e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::logLikCells(e,p);
  }
  virtual std::vector<double> logLikSamples(const Experiment* e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::logLikSamples(e,p);
  }
  virtual double logLik(const Experiment* e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::logLik(e,p);
  }
  virtual std::vector<double> epsilon(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::epsilon(e,p);
  }
  virtual const std::vector<double> weight(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_Multinomial_Model::weight(e,p);
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
  virtual CortexSimulation simulate(const Experiment* ,const Parameters &, const std::pair<std::vector<double>,std::vector<std::size_t>>&) const override
  {
    return {};
  }
  virtual std::vector<std::vector<double> > g(const Experiment*,const Parameters &, const CortexSimulation &) const override
  {
    return {};
  }
};


class CortexPoisonLikelihood: public ABC_MultiPoison_Model, public CortexLikelihood
{
  // ABC_Multinomial_Model interface
public:


  // ABC_Multinomial_Model interface

  virtual  std::vector<std::vector<double>> f(const Experiment* e, const Parameters& parameters, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const override;

  std::ostream& writeYfitHeaderDataFrame(const Experiment* e,std::ostream& os, const CortexSimulation &s)const override;



  CortexPoisonLikelihood(std::string id,
                         const Parameters& prior
                         ,double dx
                         ,double dtmin0
                         ,double dtmin
                         ,std::size_t nPoints_per_decade
                         ,double dtmax
                         ,double tequilibrio):
    CortexLikelihood(id,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio){}

  CortexPoisonLikelihood(std::string id,
                         const Parameters& prior
                         ,double dx
                         ,double dtmin0
                         ,double dtmin
                         ,std::size_t nPoints_per_decade
                         ,double dtmax
                         ,double tequilibrio
                         ,double maxlogError,
                         double dtinf):
    CortexLikelihood(id,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,maxlogError,dtinf){}



  CortexPoisonLikelihood(std::string id,
                        const Parameters &prior,
                         double dx,
                         double dtmin0,
                         double dtmin,
                         std::size_t nPoints_per_decade,
                         double dtmax,
                         double tequilibrio,
                         double maxlogError,
                         double f_maxlogError,
                         double dtinf,
                         std::size_t maxloop,
                         bool UseDerivative):
    CortexLikelihood(id,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,maxlogError,f_maxlogError,dtinf,maxloop,UseDerivative){}


  CortexPoisonLikelihood(std::string id,
                         const Parameters &prior,
                         double dx,
                         double dtmin0,
                         double dtmin,
                         std::size_t nPoints_per_decade,
                         double dtmax,
                         double tequilibrio,
                         double t_maxlogError,
                         std::size_t maxloop):
    CortexLikelihood(id,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,t_maxlogError,maxloop){}




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
  virtual std::vector<std::vector<double> > logLikCells(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::logLikCells(e,p);
  }
  virtual std::vector<double> logLikSamples(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::logLikSamples(e,p);
  }
  virtual double logLik(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::logLik(e,p);
  }
  virtual std::vector<double> epsilon(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::epsilon(e,p);
  }
  virtual const std::vector<double> weight(const Experiment *e,const std::vector<std::vector<double> > &p) const override
  {
    return ABC_MultiPoison_Model::weight(e,p);
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
  CortexSimulation simulate(const Experiment *e, const Parameters &parameters, const std::pair<std::vector<double>, std::vector<std::size_t> > &dts) const override;
  std::vector<std::vector<double> > g(const Experiment *e, const Parameters &parameters, const CortexSimulation &s) const override;

  std::vector<std::vector<double> > f
  (const Experiment* e,const Parameters &parameters,
   std::tuple<double, std::size_t, double> &dt_n_dtmax, std::size_t dts_max,
   std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const override;

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
  virtual std::ostream &extract(const Experiment *e, std::ostream &s,
                                const std::string& ="", const std::string& ="" ) const override;



  virtual bool readBody(std::string &, std::istream &, std::ostream& ) override{return false;}

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
  virtual std::ostream &extract(const Experiment* e,std::ostream &s, const std::string & ="",
                                const std::string& ="") const override;



  virtual bool readBody(std::string &, std::istream &, std::ostream& ) override{return false;}

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


class Experiment;

class MyData
{
public:
  M_Matrix<std::size_t> operator()()const
  {

    auto nobs=CL_->n_obs(e_);
    std::size_t nrows= nobs.size();
    if (nrows>0)
      {
        std::size_t ncols=nobs[0].size();
        M_Matrix<std::size_t> out(nrows,ncols);
        for (std::size_t i=0; i<nrows; ++i)
          for (std::size_t j=0; j<ncols; ++j)
            out(i,j)=nobs[i][j];
        return out;
      }
    else return {};
  }

  const Experiment* myExperiment()const {return e_;}
  MyData(const CortexLikelihood *CL, const Experiment * e):CL_(CL),e_(e){}

  MyData()=default;
  std::size_t size()const {return 1;}
private:
  const CortexLikelihood* CL_;
  const Experiment* e_;

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

























#endif // CORTEXLIKELIHOOD

