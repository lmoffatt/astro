#ifndef CORTEXLIKELIHOOD
#define CORTEXLIKELIHOOD
#include "Models.h"
#include "BayesIteration.h"


class CortexLikelihood: public BaseAgent,virtual public ABC_Distribution_Model, public ABC_Freq_obs
{
  // ABC_Multinomial_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override;

  virtual const ABC_Freq_obs &getData() const override;
  virtual const Parameters &getPrior() const override;



  CortexLikelihood(std::string id,
                   const Experiment* e,
                   const Parameters& prior
                   , double dx, double dt
                   , double tequilibrio);

  // ABC_Freq_obs interface
  virtual const std::vector<std::vector<double>> & n_obs()const override;
  virtual const std::vector<double> &n_obs(unsigned i) const override;
  virtual const std::vector<double>& bin_dens()const override;

  virtual const std::vector<double> &ntot_obs() const override;
  virtual std::ostream &put(std::ostream &s) const override;

  ~CortexLikelihood(){}

  std::ostream& write(std::ostream& s)const;


  const Experiment* getExperiment()const{return e_;}

  const BaseModel* getModel()const {return m_;}

protected:
  std::vector<std::vector<double>> getstate(const Experiment* e);

  std::vector<double> getNtot(const std::vector<std::vector<double>> nstate);

  Parameters prior_;
  BaseModel* m_;
  const Experiment* e_;
  double dx_;
  double dt_;
  double tequilibrio_;
  std::vector<std::vector<double> > nstate_;
  std::vector<double> ntot_;
  std::vector<double> bin_dens_;


  // BaseClass interface
public:
  static std::string ClassName(){return "CortexLikelihood";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual CortexLikelihood* create() const override=0;
  CortexLikelihood(){}
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"prior",prior_);
    writeField(s,"experimental_results",e_->id());
    writeField(s,"grid_legth",dx_);
    writeField(s,"sample_time",dt_);

    writeField(s,"tiempo_equilibrio",tequilibrio_);
    return s;
  }
  virtual void clear() override
  {
    prior_.clear();
    e_=nullptr;
    nstate_.clear();
    ntot_.clear();
  }
  virtual bool readBody(std::string &line, std::istream &s) override;

  std::vector<double> getNBins(const Experiment *e);
protected:
  void update(){
    m_=BaseModel::create(prior_);
    bin_dens_=getNBins(e_);
    nstate_=getstate(e_);
    ntot_=getNtot(nstate_);

  }
};


class CortexMultinomialLikelihood:public ABC_Multinomial_Model, virtual public CortexLikelihood
{


  // BaseClass interface
public:

  CortexMultinomialLikelihood(std::string id,
                              const Experiment* e,
                              const Parameters& prior
                              ,double dx
                              ,double dt
                              ,double tequilibrio):
    CortexLikelihood(id,e,prior,dx,dt,tequilibrio){}

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
                         ,double dt
                         ,double tequilibrio):
    CortexLikelihood(id,e,prior,dx,dt,tequilibrio){}

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






#endif // CORTEXLIKELIHOOD

