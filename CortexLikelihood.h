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
  virtual CortexSimulation simulate(const Parameters &parameters,const std::pair<std::vector<double>, std::vector<std::size_t>>& dts) const=0;
  virtual std::vector<std::vector<double> > g(const Parameters &parameters, const CortexSimulation &s) const=0;
  CortexLikelihood();
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override;
  virtual bool readBody(std::string &line, std::istream &s, std::ostream& logs) override;

  virtual std::ostream& writeYfitHeaderDataFrame(std::ostream& os,const CortexSimulation &)const{return os;}


  std::vector<double> getNBins(const Experiment *e);
  // ABC_Multinomial_Model interface
public:
  virtual void setPrior(const Parameters &parameters) override;

  virtual const ABC_Freq_obs &getData() const override;
  virtual const Parameters &getPrior() const override;

  CortexLikelihood(std::string id,
                   Experiment* e,
                   const Parameters& prior
                   , double dx
                   , double dtmin0,
                   double dtmin,
                   std::size_t nPoints_per_decade,
                   double dtmax
                   , double tequilibrio);

  CortexLikelihood(std::string id,
                   Experiment *e,
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


  virtual void setSimulation();

  virtual void setMeasure();

  CortexLikelihood(std::string id,
                   Experiment *e,
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
                   Experiment *e,
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

  std::vector<std::vector<double>> getstate(const Experiment* e);

  std::vector<double> getNtot(const std::vector<std::vector<double>> nstate);

  Parameters prior_;
  const BaseModel* m_;
  Experiment* e_;
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
  std::vector<std::vector<double> > nstate_;
  std::vector<double> ntot_;
  std::vector<double> bin_dens_;


};


class CortexMultinomialLikelihood:public ABC_Multinomial_Model, virtual public CortexLikelihood
{


  // BaseClass interface
public:

  CortexMultinomialLikelihood(std::string id,
                              Experiment* e,
                              const Parameters& prior
                              ,double dx
                              ,double dtmin0
                              ,double dtmin
                              ,std::size_t nPoints_per_decade
                              ,double dtmax
                              ,double tequilibrio):
    CortexLikelihood(id,e,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio){}

  static std::string ClassName(){return "CortexMultinomialLikelihood";}
  virtual std::string myClass() const override{return ClassName();}

  // ABC_Distribution_Model interface
public:
  virtual std::vector<std::vector<double> > f(const Parameters &parameters, std::pair<std::vector<double>,std::vector<std::size_t>>& dts) const override;

  virtual std::vector<std::vector<double> > f
  (const Parameters &parameters,
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
  virtual CortexSimulation simulate(const Parameters &, const std::pair<std::vector<double>,std::vector<std::size_t>>&) const override
  {
    return {};
  }
  virtual std::vector<std::vector<double> > g(const Parameters &, const CortexSimulation &) const override
  {
    return {};
  }
};


class CortexPoisonLikelihood: public ABC_MultiPoison_Model, public CortexLikelihood
{
  // ABC_Multinomial_Model interface
public:


  // ABC_Multinomial_Model interface
  virtual  std::vector<std::vector<double>> f(const Parameters& parameters, std::pair<std::vector<double>,std::vector<std::size_t>>& dts) const override;


  std::ostream& writeYfitHeaderDataFrame(std::ostream& os, const CortexSimulation &s)const override;



  CortexPoisonLikelihood(std::string id,
                         Experiment* e,
                         const Parameters& prior
                         ,double dx
                         ,double dtmin0
                         ,double dtmin
                         ,std::size_t nPoints_per_decade
                         ,double dtmax
                         ,double tequilibrio):
    CortexLikelihood(id,e,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio){}

  CortexPoisonLikelihood(std::string id,
                         Experiment* e,
                         const Parameters& prior
                         ,double dx
                         ,double dtmin0
                         ,double dtmin
                         ,std::size_t nPoints_per_decade
                         ,double dtmax
                         ,double tequilibrio
                         ,double maxlogError,
                         double dtinf):
    CortexLikelihood(id,e,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,maxlogError,dtinf){}



  CortexPoisonLikelihood(std::string id,
                         Experiment *e,
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
    CortexLikelihood(id,e,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,maxlogError,f_maxlogError,dtinf,maxloop,UseDerivative){}


  CortexPoisonLikelihood(std::string id,
                         Experiment *e,
                         const Parameters &prior,
                         double dx,
                         double dtmin0,
                         double dtmin,
                         std::size_t nPoints_per_decade,
                         double dtmax,
                         double tequilibrio,
                         double t_maxlogError,
                         std::size_t maxloop):
    CortexLikelihood(id,e,prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,t_maxlogError,maxloop){}




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
  CortexSimulation simulate(const Parameters &parameters, const std::pair<std::vector<double>, std::vector<std::size_t> > &dts) const override;
  std::vector<std::vector<double> > g(const Parameters &parameters, const CortexSimulation &s) const override;

  std::vector<std::vector<double> > f
  (const Parameters &parameters,
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
  virtual std::ostream &extract(std::ostream &s, const std::string & ="",
                                const std::string& ="") const override;



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
  virtual std::ostream &extract(std::ostream &s, const std::string & ="",
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
  mcmc_prior<double> prior(const Data& , M_Matrix<double> param)const
  {
    mcmc_prior<double> out;
    out.param=param;

    auto& p=CL_->getPrior();
    //    std::size_t npar=param.size();

    Parameters Pa=p.toParameters(param.toVector());

    out.logPrior=p.logProb(Pa);
    out.D_prior.H=M_Matrix<double>(p.getInvCovariance());
    M_Matrix<double> d(param.nrows(),param.ncols(),p.trMeans());
    d=param - d;


    out.D_prior.G=d*out.D_prior.H;

    return out;

  }
  M_Matrix<double> f(const Data& , M_Matrix<double> param, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)
  const{
    auto ff=CL_->f(param.toVector(),dts);
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

  M_Matrix<double> f(const Data& , M_Matrix<double> param,
                     double& dtmin,std::size_t& nper10,double& dtmax,
                     std::size_t dts_max,
                     std::pair<std::vector<double>,std::vector<std::size_t>>& dts)
  const{
    std::tuple<double,std::size_t,double> dtmin_n_dtmax{dtmin,nper10,dtmax};
    auto ff=CL_->f(param.toVector(),dtmin_n_dtmax,dts_max,dts);
    dtmin=std::get<0>(dtmin_n_dtmax);
    nper10=std::get<1>(dtmin_n_dtmax);
    dtmax=std::get<2>(dtmin_n_dtmax);

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



  M_Matrix<double> logLanda(const Data& data, M_Matrix<double> param, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const{
    M_Matrix<double> out=f(data,param,dts);
    out=out.apply([](double x){return log10_guard(x);});
    return out;

  }

  std::string id()const {return CL_->getPrior().id();}


  Parameters getPrior()const
  {
    return CL_->getPrior();
  }


  Parameters getParameter(const M_Matrix<double>& par)const
  {
    return CL_->getPrior().toParameters(par.toVector());
  }

  CortexSimulation getSimulation(const Parameters& par,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const
  {
    return CL_->simulate(par,dts);
  }

  CortexSimulation getSimulation(const M_Matrix<double>& par,std::pair<std::vector<double>, std::vector<std::size_t>> dts)const
  {
    auto CLc= const_cast<CortexLikelihood*> (CL_);
    CLc->setSimulation();
    auto out= CLc->simulate(getParameter(par),dts);
    CLc->setMeasure();
    return out;
  }

  const CortexLikelihood& getLikelihood()const
  {
    return *CL_;
  }


  MyModel(CortexLikelihood* CL):CL_(CL){}
private:
  const CortexLikelihood* CL_;
};


template<class Data=std::vector<MyData>>
class MyRandomEffectsModel
{
public:
  typedef  M_Matrix<M_Matrix<double>>   E;

  M_Matrix<E> sample(std::mt19937_64& mt, std::size_t n)const
  {
    Parameters hyperParameter=CL_->getPrior().randomHiperSample(mt,1);

    M_Matrix<E> out(1,2);
    out[0]=hyperParameter.hyperParameters();

    E   samples(1,n);
    for (std::size_t i=0; i<n; ++i)
      {
        auto sample=hyperParameter.randomSample(mt,1.0);
        samples[i]=M_Matrix<double>(1,sample.size(),sample.trMeans());

      }
    out[1]=std::move(samples);
    return out;
  }
  mcmc_prior<E> prior(const Data& , M_Matrix<E> param)const
  {
    mcmc_prior<E> out;
    out.param=param;

    auto& p=CL_->getPrior();
    std::size_t npar=param.size();

    std::size_t nrep=param[1].size();
    Parameters Pa=p.toHyperParameters(param[0]);
    double logPrior=p.logHiperProb(Pa);
    out.D_prior.H=M_Matrix<E>(2,2,M_Matrix<E>::SYMMETRIC);
    out.D_prior.G=M_Matrix<E>(1,2,M_Matrix<E>::FULL);

    auto dm=-M_Matrix<double>(1,npar,p.trMeans());
    dm+=param[0][0];

    auto ds=-M_Matrix<double>(1,npar,p.mean_of_log_std());
    ds+=param[0][1];

    out.D_prior.H(0,0)=p.getHyperHessian();

    out.D_prior.G[0][0]=dm*out.D_prior.H(0,0)(0,0);
    out.D_prior.G[0][1]=ds*out.D_prior.H(0,0)(1,1);


    out.D_prior.H(1,1)=E(nrep,nrep,E::SCALAR_DIAGONAL);
    out.D_prior.G[1]=E(1,nrep);




    out.D_prior.H(1,1)[0]=Pa.getHessian();

    out.D_prior.H(0,1)=E(2,nrep,E::FULL);



    double logPriorEff=0;
    for (std::size_t i=0; i<nrep; ++i)
      {
        Parameters P_i=p.toParameters(param[1][i].toVector());
        logPriorEff+=Pa.logProb(P_i);
        M_Matrix<double> d=param[1][i]-param[0][0];
        out.D_prior.G[1][i]=d*out.D_prior.H(1,1)[0];
        out.D_prior.G[0][0]-=out.D_prior.G[1][i];
        out.D_prior.G[0][1]-=diag(d)*out.D_prior.G[1][i]+eye<double>(npar);

        out.D_prior.H(0,0)(0,0)+=out.D_prior.H(1,1)[0];

        out.D_prior.H(0,1)(0,i)=-out.D_prior.H(1,1)[0];
        out.D_prior.H(0,1)(1,i)=
            d*out.D_prior.H(0,1)(0,i)*2.0;


        out.D_prior.H(0,0)(0,1)-=out.D_prior.H(0,1)(1,i);
        out.D_prior.H(0,0)(1,1)-=diag(d)*out.D_prior.H(0,1)(1,i);


      }

    out.logPrior=logPrior+logPriorEff;

    return out;

  }
  M_Matrix<E> f(const Data& , M_Matrix<E> param, std::vector<std::pair<std::vector<double>,std::vector<std::size_t>>>& dts)
  const
  {
    std::size_t nrep=param[1].size();
    M_Matrix<E> out(1,1,E(nrep,1));
    for (std::size_t ir=0; ir<nrep; ++ir)
      {
        auto ff=CL_->f(param[1][ir].toVector(),dts[ir]);
        std::size_t nrows= ff.size();
        if (nrows>0)
          {
            std::size_t ncols=ff[0].size();
            out[0][ir]=M_Matrix<double>(nrows,ncols);
            for (std::size_t i=0; i<nrows; ++i)
              for (std::size_t j=0; j<ncols; ++j)
                out[0][ir](i,j)=ff[i][j];
          }
        else return {};
      }
    return out;
  }




  M_Matrix<double> logLanda(const Data& data, M_Matrix<double> param, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const{
    M_Matrix<double> out=f(data,param,dts);
    out=out.apply([](double x){return log10_guard(x);});
    return out;

  }

  std::string id()const {return CL_->getPrior().id();}


  Parameters getPrior()const
  {
    return CL_->getPrior();
  }


  Parameters getParameter(const M_Matrix<double>& par)const
  {
    return CL_->getPrior().toParameters(par.toVector());
  }

  CortexSimulation getSimulation(const Parameters& par,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const
  {
    return CL_->simulate(par,dts);
  }

  CortexSimulation getSimulation(const M_Matrix<double>& par,std::pair<std::vector<double>, std::vector<std::size_t>> dts)const
  {
    auto CLc= const_cast<CortexLikelihood*> (CL_);
    CLc->setSimulation();
    auto out= CLc->simulate(getParameter(par),dts);
    CLc->setMeasure();
    return out;
  }

  const CortexLikelihood& getLikelihood()const
  {
    return *CL_;
  }


  MyRandomEffectsModel(CortexLikelihood* CL):CL_(CL){}
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




template <class Ad>
using TI=
Thermodynamic_Integration_mcmc<Ad,
MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa,LevenbergMarquardt_step>;

template <class Ad>
using TT=
Template_Tempering_mcmc<Ad,
MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa,LevenbergMarquardt_step> ;





















#endif // CORTEXLIKELIHOOD

