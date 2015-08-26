#ifndef CORTEXLIKELIHOOD
#define CORTEXLIKELIHOOD
#include "Models.h"
#include "BayesIteration.h"


class CortexModelLikelihood: public BaseAgent,public ABC_Multinomial_Model, public ABC_Freq_obs
{
  // ABC_Multinomial_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override;

  virtual const ABC_Freq_obs &getData() const override;
  virtual const Parameters &getPrior() const override;


  // ABC_Multinomial_Model interface
  virtual std::vector<std::vector<double> > p_exp(const Parameters &parameters) const override;



  CortexModelLikelihood(std::string id,
                        const Experiment* e,
                        const Parameters& prior
                        ,double dt
                        ,double tequilibrio);

  // ABC_Freq_obs interface
  virtual const std::vector<std::vector<double>> & n_obs()const override;
  virtual const std::vector<double> &n_obs(unsigned i) const override;
  virtual const std::vector<double> &ntot_obs() const override;
  virtual std::ostream &put(std::ostream &s) const override;

  ~CortexModelLikelihood(){}

  std::ostream& write(std::ostream& s)const;


  const Experiment* getExperiment()const{return e_;}


private:
  std::vector<std::vector<double>> getstate(const Experiment* e);

  std::vector<double> getNtot(const std::vector<std::vector<double>> nstate);

  Parameters prior_;
  BaseModel* m_;
  const Experiment* e_;
  double dt_;
  double tequilibrio_;
  std::vector<std::vector<double> > nstate_;
  std::vector<double> ntot_;


  // BaseClass interface
public:
  static std::string ClassName(){return "CortexLikelihood";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual CortexModelLikelihood* create() const override
  {
    return new CortexModelLikelihood;
  }
  CortexModelLikelihood(){}
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"prior",prior_);
    writeField(s,"experimental_results",e_->id());
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

protected:
  void update(){
    m_=BaseModel::create(prior_);
    nstate_=getstate(e_);
    ntot_=getNtot(nstate_);

  }
};





class CortexLikelihoodEvaluation: public BaseAgent
{
  // BaseClass interface
public:
  static std::string ClassName(){return "CortexLikelihoodEvaluation";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  CortexLikelihoodEvaluation(){}
  ~CortexLikelihoodEvaluation(){}
  virtual CortexLikelihoodEvaluation *create() const override
  {
    return new CortexLikelihoodEvaluation;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override{}
  virtual std::ostream &extract(std::ostream &s, const std::string & s1="" ,
                                const std::string& s2="") const override;



  virtual bool readBody(std::string &line, std::istream &s) override{}

  CortexLikelihoodEvaluation(const CortexModelLikelihood& CL,
                             const Parameters& p)
    :
      CL_(&CL)
    ,p_(p)
  {}





protected:
  virtual void update() override{}
  const CortexModelLikelihood* CL_;
  Parameters p_;

  //--------Calculated variables-----------------//

  std::vector<std::vector<double> > p_exp_;




};







#endif // CORTEXLIKELIHOOD

