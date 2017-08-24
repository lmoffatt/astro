#ifndef MODELS
#define MODELS

#include "CortexMeasure.h"
#include "CortexSimulation.h"
#include "myTuples.h"
//#include "BayesIteration.h"
#include "MatrixBanded.h"
#include <map>
#include <vector>
#include <numeric>
#include <iostream>
inline
double MultinomialLikelihood(const std::vector<double>& n, const std::vector<double>& p, double ntot=0)
{
  if (n.size()!=p.size())
    return std::numeric_limits<double>::quiet_NaN();
  if (ntot==0)
    for (double ni:n)
      ntot+=ni;
  double lik=0;
  for (unsigned i=0; i< n.size(); ++i)
    {
      double pLik=0;
      if (n[i]>0)
        pLik=n[i]*log(ntot*p[i]/n[i]);
      //  std::cout<<"i="<<i<<" obs="<<n[i]<<" exp="<<ntot*p[i]<<" pLik="<<pLik<<"\n";
      lik+=pLik;
    }
  return lik;
}





class SimplestModel
{
public:

  class Param
  {
  public:
    double Dpsi_;
    double Domega_;
    double epsilon_;

    double kon_psi_;
    double kcat_psi_;
    double kon_omega_;
    double kcat_omega_;

    std::vector<double> Keq_gmax_psi_;
    std::vector<double> Keq_gmax_omega_;



    std::vector<double> ksig_omega_;

    std::vector<double> ksig_max_omega_;
    std::vector<double> ksig_max_psi_;

    std::vector<double> g_left_;
    std::vector<double> g_rigth_;

    std::vector<double> g_max_omega_;
    std::vector<double> g_max_psi_;

    std::vector<double> a_;

    std::vector<double> a_omega_;
    std::vector<double> a_psi_;


    std::vector<double> N_;
    std::vector<double> M_;

    double dens_Astr_;
    double dens_Neur_;
    double dens_Microglia_;

    double inj_width_;
    double DAMP_psi_ratio_;
    double DAMP_omega_ratio_;

    double prot_concentration_;
    double DAMP_MW_;



  };

  SimplestModel();



  CortexState init(const Param &p, const CortexExperiment &s) const;

  CortexSimulation simulate(const Parameters& par,const Param &p, const CortexExperiment &sp, double dt)const;

  CortexSimulation simulate(Parameters par,  Param p, const Experiment &sp, double dx, double dtmin0,  double dtmin, std::size_t nPoints_per_decade, double dtmax, double tequilibrio)const;

  std::pair<CortexSimulation, std::vector<double> >
  simulate_Adapted(Parameters par,  Param p, const Experiment &sp, double dx, double dtmin0,  double dtmin, std::size_t nPoints_per_decade, double dtmax, double tequilibrio, double maxLogError, double dtinf)const;



  CortexState init(const SimplestModel::Param &p, const Experiment &s, double dx) const;

  std::pair<CortexState, double> nextEuler_Adapt(const SimplestModel::Param &p, const std::pair<CortexState, double> &c, double dt, double dtmin, double maxlogError, std::vector<double> &dts) const;

  void nextEuler(const SimplestModel::Param &p, CortexState &c, double dt) const;

  std::pair<Di<CortexState>, double> nextEuler_Adapt_i(const Di<SimplestModel::Param> &p, const std::pair<Di<CortexState>, double> &c, double dt, double dtmin, double maxlogError) const;


  Di<CortexSimulation> simulate_Adapted_i(const Di<Parameters> &par, const Di<Param> &p, const Experiment &sp, double dx, double desiredError, double dt0, std::size_t nPoints_per_decade, double dtmax,  double dtmin, double tequilibrio) const;
  CortexSimulation simulate(Parameters par, Param p, const Experiment &sp, double dx, const std::pair<std::vector<double>,std::vector<std::size_t>> &dts, double tequilibrio) const;

  void CrankNicholsonStep(const CortexState::Der &d, const CortexState::dDer &dd, const SimplestModel::Param &p, CortexState &c, double dt) const;

  std::pair<CortexState, std::pair<std::size_t,double> > CrankNicholsonStep(const CortexState::Der &d, const CortexState::dDer &dd, const SimplestModel::Param &p, const CortexState &c, double dt, double maxlogError, std::size_t maxloop) const;

  std::pair<CortexState, double> nextCrankNicholson_Adapt(const SimplestModel::Param &p, const std::pair<CortexState, double> &c, double dt,  double dtmin, double maxlogError, double maxlogErrorCN, std::pair<std::vector<double>, std::vector<std::size_t> > &dts, std::size_t maxloop, bool useDerivative) const;
  std::vector<double> dPsi_Bound_dPsi(const SimplestModel::Param &p, const CortexState &c) const;
  CortexState::dDer ddStep(const SimplestModel::Param &p, const CortexState &c) const;
  std::vector<std::vector<double> > dPsi_Bound_dRho(const SimplestModel::Param &p, const CortexState &c) const;
  std::vector<double> dOmega_Bound_dOmega(const SimplestModel::Param &p, const CortexState &c) const;
  std::vector<std::vector<double> > domega_Bound_dRho(const SimplestModel::Param &p, const CortexState &c) const;



  std::pair<CortexSimulation, std::pair<std::vector<double>,std::vector< std::size_t>>> simulate_CrankNicholson_Adapted(Parameters par, SimplestModel::Param p, const Experiment &sp, double dx, double dtmin0,  double dtmin, std::size_t nPoints_per_decade, double dtmax, double tequilibrio, double maxLogError, double maxlogErrorCN, double dtinfratio, std::size_t maxloop, bool UseDerivative) const;

  std::pair<CortexState, double> CrankNicholson_Adapt_2(const CortexState &one, const CortexState::Der &d, const CortexState::dDer &dd, const SimplestModel::Param &p, const std::pair<CortexState, double> &c, double dt, double dtmin, double maxlogError, double maxlogErrorCN, std::size_t maxloop, std::pair<std::vector<double>, std::vector<std::size_t> > &dts) const;

  std::pair<CortexState, std::pair<double, std::size_t> > nextCrankNicholson(const SimplestModel::Param &p, const CortexState &c, double dt, double maxlogErrorCN, std::size_t maxloop, bool UseDerivative) const;

  std::pair<CortexSimulation, std::pair<std::vector<double>, std::vector<std::size_t> > >  simulate_CrankNicholson(Parameters par, SimplestModel::Param p, const Experiment &sp, double dx, double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double tequilibrio, double maxlogErrorCN, std::size_t maxloop,bool UseDerivative) const;

  CortexSimulation simulate_CrankNicholson_dt
  (Parameters par, Param p, const Experiment &sp, double dx, const std::pair<std::vector<double>, std::vector<std::size_t> > &dts, double tequilibrio, bool UseDerivative) const;

  std::pair<CortexState, std::pair<std::size_t, double> > CrankNicholsonStep(const CortexState::Der &d, const SimplestModel::Param &p, const CortexState &c, double dt, double maxlogError, std::size_t maxloop) const;

  std::pair<CortexState, double> CrankNicholson_Adapt_2(const CortexState &one, const CortexState::Der &d, const SimplestModel::Param &p, const std::pair<CortexState, double> &c, double dt, double dtmin, double maxlogError, double maxlogErrorCN, std::size_t maxloop, std::pair<std::vector<double>, std::vector<std::size_t> > &dts) const;
  std::pair<CortexSimulation, std::pair<std::vector<double>, std::vector<std::size_t> > > simulate_Euler(Parameters par, Param p, const Experiment &sp, double dx, double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, std::size_t dts_max,double tequilibrio) const;
private:
  void addDamp(CortexState& c,const Param& p)const;

  std::vector<double> dPsi_dt(const Param &p, const CortexState& c) const;


  void Psi_Bound(const Param &p, CortexState &c) const;


  std::vector<double> dOmega_dt(const Param &p, const CortexState &c)const;

  void Omega_Bound(const Param &p,  CortexState& c) const;


  std::vector<std::vector<double> > dRho_dt(const Param &p, const CortexState &c, bool hasOmega)const;
  std::vector<CortexState> nextEuler_i(const std::vector<SimplestModel::Param> &p, const std::vector<CortexState> &c, double dt) const;

  std::pair<CortexState, double> nextEuler_Adapt_2(const CortexState &one, const CortexState::Der &d, const SimplestModel::Param &p, const std::pair<CortexState, double> &c, double dt, double dtmin, double maxlogError, std::vector<double> &dts) const;


  std::pair<Di<CortexState>, double> nextEuler_Adapt_2_i(CortexState one, const CortexState::Der &d, const Di<SimplestModel::Param> &p, const std::pair<Di<CortexState>, double> &c, double dt,  double dtmin, double maxlogError) const;

  void EulerStep(const CortexState::Der &d, const SimplestModel::Param &p, CortexState &c, double dt) const;

  CortexState::Der dStep(const SimplestModel::Param &p, const CortexState &c) const;


};





class BaseModel
{
public:
  virtual std::string id()const=0;

  virtual Parameters getParameters()const=0;

  virtual void loadParameters(const Parameters& p)=0;

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,  double tequilibrio) const=0;




  virtual const SimplestModel& myModel()const=0;
  virtual SimplestModel& myModel()=0;

  virtual SimplestModel::Param toModelParameters()const=0;

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const=0;


  std::pair<CortexSimulation, std::pair<std::vector<double>,std::vector< std::size_t>>> run_CN_adapt(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double maxlogErrorCN,double dtinf, std::size_t maxloop, bool Usederivative) const
  {
    return myModel().simulate_CrankNicholson_Adapted(getParameters(),toModelParameters(),e,dx,dtmin0,dtmin,
                                                     nPoints_per_decade,dtmax,teq,maxlogError,maxlogErrorCN,dtinf,maxloop,Usederivative);

  }

  virtual CortexSimulation run_CN_dt(const Experiment& e,double dx,const std::pair<std::vector<double>, std::vector<std::size_t>>& dts,double tequilibrio, bool UseDerivative) const
  {
    return myModel().simulate_CrankNicholson_dt(getParameters(),toModelParameters(),e,dx,dts,tequilibrio,UseDerivative);

  }
  std::pair<CortexSimulation, std::pair<std::vector<double>, std::vector<std::size_t> > > run_CN(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,  double tequilibrio,double maxlogError,std::size_t maxloop, bool UseDerivative) const
  {
    return myModel().simulate_CrankNicholson(getParameters(),toModelParameters(),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,maxlogError,maxloop,UseDerivative);
  }

  std::pair<CortexSimulation, std::pair<std::vector<double>, std::vector<std::size_t> > >
  run_Euler(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,  std::size_t dts_max,double tequilibrio) const
  {
    return myModel().simulate_Euler
        (getParameters(),toModelParameters(),e,
         dx,dtmin0,dtmin,nPoints_per_decade,dtmax,dts_max,tequilibrio);

  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,  double tequilibrio,double maxlogError,double dtinf) const=0;

  virtual BaseModel* clone()const =0;

  static BaseModel* create(const Parameters& p);

  virtual ~BaseModel(){}

  virtual std::vector<double> getObservedProbFromModel(const std::vector<double>& modelRho)const;

  virtual std::vector<double> getObservedNumberFromModel(const std::vector<double>& modelRho)const;

  virtual std::vector<double> getNumberAtInjuryFromModel(const std::vector<double>& modelRho,
                                                         double f)const;

  virtual std::vector<double> getObservedNumberFromData(const std::vector<double>& modelRho)const;


  virtual std::vector<std::string> getModelStateLabels()const;

  virtual std::vector<std::string> getObservedStateLabels()const;

  virtual std::vector<std::string> getApoptoticStatesAtInjuryLabels()const;


  //  double likelihood(const Experiment& m, const CortexSimulation& s)const
  //  {
  //    unsigned is=0;
  //    double sumLogLik=0;

  //    for (unsigned ie=0; ie<m.numMeasures(); ie++)
  //      {
  //        const CortexMeasure* cm=m.getMeasure(ie);
  //        double t=cm->dia()*24*60*60;
  //        while (s.t_[is]<t
  //               &&is<s.t_.size())
  //          ++is;
  //        std::vector<double> rho_sim;
  //        std::vector<double> rho_meas;
  //        for (unsigned ix=0; ix<cm->dx().size(); ++ix)
  //          {
  //            rho_sim=getObservedProbFromModel(s.rho_[is][ix]);
  //            auto ob=cm->meanAstro(ix);
  //            rho_meas=getObservedNumberFromData(ob);
  //            sumLogLik+=MultinomialLikelihood(rho_meas,rho_sim);
  //          }

  //      }
  //    return sumLogLik;

  //  }






  virtual std::size_t getNumberOfObservedStates() const;

  virtual std::size_t getNumberOfSimulatedStatesAtInjury() const;
private:
  static std::map<double,BaseModel*> models_;
  static std::map<double,BaseModel*> getModels();
};





class MicrogliaModel: public BaseModel
{
public:

  virtual ~MicrogliaModel(){}

  std::vector<std::string> getModelStateLabels() const override
  {
    return {"Neuron","M2","M1","type0","typeI","typeII","typeIII","typeIV","typeV"};
  }

  virtual std::vector<std::string> getApoptoticStatesAtInjuryLabels()const override
  {
    return {"Neuron","type0","typeI","typeII","typeIII","typeIV","typeV"};

  }


  virtual std::vector<double> getObservedProbFromModel(const std::vector<double>& modelRho)const override
  {
    std::vector<double>  v(5);
    double n=modelRho[4]+modelRho[5]+modelRho[6]+modelRho[7]+modelRho[8];
    v[0]=modelRho[4]/n;
    v[1]=modelRho[5]/n;
    v[2]=modelRho[6]/n;
    v[3]=modelRho[7]/n;
    v[4]=modelRho[8]/n;
    return v;
  }
  virtual std::vector<double> getObservedNumberFromModel(const std::vector<double>& modelRho)const override
  {
    std::vector<double>  v(5);
    v[0]=modelRho[4];
    v[1]=modelRho[5];
    v[2]=modelRho[6];
    v[3]=modelRho[7];
    v[4]=modelRho[8];
    return v;

  }

  std::vector<double> getNumberAtInjuryFromModel(const std::vector<double> &modelRho, double f) const override
  {
    std::vector<double>  v(7);
    v[0]=modelRho[0]*f;
    v[1]=modelRho[3]*f;
    v[2]=modelRho[4]*f;
    v[3]=modelRho[5]*f;
    v[4]=modelRho[6]*f;
    v[5]=modelRho[7]*f;
    v[6]=modelRho[8]*f;
    return v;
  }



  virtual std::vector<double> getObservedNumberFromData(const std::vector<double>& modelRho)const override
  {
    std::vector<double>  v(5);
    v[0]=modelRho[0];
    v[1]=modelRho[1];
    v[2]=modelRho[2];
    v[3]=modelRho[3];
    v[4]=modelRho[4];
    return v;
  }


  virtual std::size_t getNumberOfSimulatedStatesAtInjury() const override
  {
    return 7;
  }

};















class Model00:public BaseModel
{
  SimplestModel m;

  Model00* clone()const override { return new Model00;}



  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;


    double g_23_;
    double g_32_;

    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);

    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;
    s.g_left_[4]=p.g_32_;
    s.g_left_[5]=p.g_32_;
    s.g_left_[6]=p.g_32_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model00(){}
  ~Model00(){}
  virtual std::string id()  const override
  {
    return "Model 0.0";
  }
  static double number()
  {
    return 0;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );
    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters  (const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");


    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");


  }



  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }




  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }
  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run
  (const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade,
   double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model00(const Parameters& p)
  {
    loadParameters(p);
  }



};



class Model00m:public MicrogliaModel
{
  SimplestModel m;

  Model00m* clone()const  override { return new Model00m;}



  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_M1M2_;
    double g_M2M1_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;


    double g_23_;
    double g_32_;


    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;


    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(9,0);

    s.ksig_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_psi_=std::vector<double>(9,0);

    s.g_left_=std::vector<double> (9,0);

    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_32_;
    s.g_left_[8]=p.g_32_;




    s.g_rigth_=std::vector<double> (9,0);

    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_23_;
    s.g_rigth_[7]=p.g_23_;

    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_;

    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;
    s.g_max_psi_[6]=p.g_max_;
    s.g_max_psi_[7]=p.g_max_;


    s.a_=std::vector<double> (9,0);
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_2_*p.a_factor_;
    s.a_[7]=s.a_[6]*p.a_factor_;
    s.a_[8]=s.a_[7]*p.a_factor_;

    s.a_omega_=std::vector<double> (9,0);
    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_N_;
    s.N_[2]=p.N_N_;

    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_0_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_2_*1.5;
    s.N_[7]=p.N_2_*3;
    s.N_[8]=p.N_2_*6;



    s.M_=std::vector<double> (9,0);



    s.dens_Astr_=p.N_Astr_;

    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model00m(){}
  ~Model00m(){}
  virtual std::string id() const override
  {
    return "Model_10.0";
  }
  static double number()
  {
    return 10.0;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_M2M1",p_.g_M2M1_);

    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );
    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );

    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);






    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_M1M2_=p.get("g_M1M2");
    p_.g_M2M1_=p.get("g_M2M1");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");


    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");

    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }




  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }



  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model00m(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model011:public BaseModel
{
  SimplestModel m;

  Model011* clone()const override { return new Model011;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_34_;
    double g_45_;

    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model011(){}
  ~Model011(){}
  virtual std::string id() const override
  {
    return "Model 0.11";
  }
  static double number()
  {
    return 0.11;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );

    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );




    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }


  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model011(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};





class Model012:public BaseModel
{
  SimplestModel m;

  Model012* clone()const override { return new Model012;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model012(){}
  ~Model012(){}
  virtual std::string id() const  override
  {
    return "Model 0.12";
  }
  static double number()
  {
    return 0.12;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model012(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};



class Model012_22:public BaseModel
{
  SimplestModel m;

  Model012_22* clone()const override { return new Model012_22;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model012_22(){}
  ~Model012_22(){}
  virtual std::string id() const  override
  {
    return "Model 0.12022";
  }
  static double number()
  {
    return 0.12022;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );




    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");



    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model012_22(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};



class Model012_51:public BaseModel
{
  SimplestModel m;

  Model012_51* clone()const override { return new Model012_51;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;

    double Kg_2_;
    double Kg_3_;
    double Kg_4_;
    double Kg_5_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);
    s.Keq_gmax_psi_[2]=p.Kg_2_;
    s.Keq_gmax_psi_[3]=p.Kg_3_;
    s.Keq_gmax_psi_[4]=p.Kg_4_;
    s.Keq_gmax_psi_[5]=p.Kg_5_;


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model012_51(){}
  ~Model012_51(){}
  virtual std::string id() const  override
  {
    return "Model 0.12051";
  }
  static double number()
  {
    return 0.12051;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("Kg_2",p_.Kg_2_);
    out.push_back("Kg_3",p_.Kg_3_);
    out.push_back("Kg_4",p_.Kg_4_);
    out.push_back("Kg_5",p_.Kg_5_);

    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.Kg_2_=p.get("Kg_2");
    p_.Kg_3_=p.get("Kg_3");
    p_.Kg_4_=p.get("Kg_4");
    p_.Kg_5_=p.get("Kg_5");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }



  Model012_51(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model013:public BaseModel
{
  SimplestModel m;

  Model013* clone()const override { return new Model013;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;


    double g_23_;
    double g_34_;
    double g_45_;

    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model013(){}
  ~Model013(){}
  virtual std::string id() const  override
  {
    return "Model 0.13";
  }
  static double number()
  {
    return 0.13;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }



  Model013(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model013_23:public BaseModel
{
  SimplestModel m;

  Model013_23* clone()const override { return new Model013_23;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;


    double g_23_;
    double g_34_;
    double g_45_;

    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_Neuron_;

    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);


    s.a_psi_[0]=p.a_max_Neuron_;


    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model013_23(){}
  ~Model013_23(){}
  virtual std::string id() const  override
  {
    return "Model 0.13023";
  }
  static double number()
  {
    return 0.13023;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );


    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );




    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.a_max_Neuron_=p.get("a_max_Neuron");


    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model013_23(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};




class Model013_23_31:public BaseModel
{
  SimplestModel m;

  Model013_23_31* clone()const override { return new Model013_23_31;}

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;
    double g_23_;
    double g_32_;


    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_Neuron_;

    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;
    s.g_left_[4]=p.g_32_;
    s.g_left_[5]=p.g_32_;
    s.g_left_[6]=p.g_32_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);


    s.a_psi_[0]=p.a_max_Neuron_;


    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_1_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_3_;
    s.N_[5]=p.N_4_;
    s.N_[6]=p.N_5_;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model013_23_31(){}
  ~Model013_23_31(){}
  virtual std::string id() const  override
  {
    return "Model 0.13023031";
  }
  static double number()
  {
    return 0.13023031;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );


    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );




    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.a_max_Neuron_=p.get("a_max_Neuron");


    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }




  Model013_23_31(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model013_23_31m:public MicrogliaModel
{
  SimplestModel m;

  Model013_23_31m* clone()const override { return new Model013_23_31m;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;

    double g_M1M2_;
    double g_M2M1_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;
    double g_23_;
    double g_32_;

    double g_max_M2_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;


    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_Neuron_;

    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(9,0);


    s.ksig_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_psi_=std::vector<double>(9,0);

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;
    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_32_;
    s.g_left_[8]=p.g_32_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;


    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_23_;
    s.g_rigth_[7]=p.g_23_;

    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_M2_;
    s.g_max_psi_[4]=p.g_max_2_;
    s.g_max_psi_[5]=p.g_max_3_;
    s.g_max_psi_[6]=p.g_max_4_;
    s.g_max_psi_[7]=p.g_max_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_omega_=std::vector<double> (9,0);
    s.a_psi_=std::vector<double> (9,0);


    s.a_psi_[0]=p.a_max_Neuron_;


    s.a_psi_[4]=p.a_max_1_;
    s.a_psi_[5]=p.a_max_2_;
    s.a_psi_[6]=p.a_max_3_;
    s.a_psi_[7]=p.a_max_4_;
    s.a_psi_[8]=p.a_max_5_;




    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;

    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;



    s.M_=std::vector<double> (9,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model013_23_31m(){}
  ~Model013_23_31m(){}
  virtual std::string id() const  override
  {
    return "Model 0.13023031m";
  }
  static double number()
  {
    return 10.13023031;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );

    out.push_back("g_max_M2",p_.g_max_M2_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_M2",p_.N_M2_ );
    out.push_back("N_M1",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );


    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );




    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);






    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_max_M2_=p.get("g_max_M2");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.a_max_Neuron_=p.get("a_max_Neuron");


    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");


    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }



  Model013_23_31m(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model013_51:public BaseModel
{
  SimplestModel m;

  Model013_51* clone()const override { return new Model013_51;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;

    double Kg_2_;
    double Kg_3_;
    double Kg_4_;
    double Kg_5_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;


    double g_23_;
    double g_34_;
    double g_45_;

    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);
    s.Keq_gmax_psi_[2]=p.Kg_2_;
    s.Keq_gmax_psi_[3]=p.Kg_3_;
    s.Keq_gmax_psi_[4]=p.Kg_4_;
    s.Keq_gmax_psi_[5]=p.Kg_5_;


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model013_51(){}
  ~Model013_51(){}
  virtual std::string id() const  override
  {
    return "Model 0.13051";
  }
  static double number()
  {
    return 0.13051;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);

    out.push_back("Kg_2",p_.Kg_2_);
    out.push_back("Kg_3",p_.Kg_3_);
    out.push_back("Kg_4",p_.Kg_4_);
    out.push_back("Kg_5",p_.Kg_5_);

    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.Kg_2_=p.get("Kg_2");
    p_.Kg_3_=p.get("Kg_3");
    p_.Kg_4_=p.get("Kg_4");
    p_.Kg_5_=p.get("Kg_5");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }





  Model013_51(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model021:public BaseModel
{
  SimplestModel m;

  Model021* clone()const override { return new Model021;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;



    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;


    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;






    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model021(){}
  ~Model021(){}
  virtual std::string id() const  override
  {
    return "Model 0.21";
  }
  static double number()
  {
    return 0.21;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );


    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_=p.get("a_max");
    p_.a_max_Neuron_=p.get("a_max_Neuron");


    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }





  Model021(const Parameters& p)
  {
    loadParameters(p);
  }


};





class Model022:public BaseModel
{
  SimplestModel m;

  Model022* clone()const override { return new Model022;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;

    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);


    s.a_psi_[0]=p.a_max_Neuron_;


    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model022(){}
  ~Model022(){}
  virtual std::string id() const  override
  {
    return "Model 0.22";
  }
  static double number()
  {
    return 0.22;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );


    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");


    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }





  Model022(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model023:public BaseModel
{
  SimplestModel m;

  Model023* clone()const override { return new Model023;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_Neuron_;

    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);


    s.a_psi_[0]=p.a_max_Neuron_;


    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model023(){}
  ~Model023(){}
  virtual std::string id() const  override
  {
    return "Model 0.23";
  }
  static double number()
  {
    return 0.23;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );


    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.a_max_Neuron_=p.get("a_max_Neuron");


    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");

    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }





  Model023(const Parameters& p)
  {
    loadParameters(p);
  }


};





class Model031:public BaseModel
{
  SimplestModel m;

  Model031* clone()const override { return new Model031;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);
    s.Keq_gmax_omega_=std::vector<double>(7,0);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;




    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_1_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_3_;
    s.N_[5]=p.N_4_;
    s.N_[6]=p.N_5_;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model031(){}
  ~Model031(){}
  virtual std::string id() const  override
  {
    return "Model 0.31";
  }
  static double number()
  {
    return 0.31;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }





  Model031(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};


class Model051:public BaseModel
{
  SimplestModel m;

  Model051* clone()const override { return new Model051;}





  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;

    double Kg_2_;
    double Kg_3_;
    double Kg_4_;
    double Kg_5_;


    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;


    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double inj_width_3_;
    double inj_width_7_;

  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_);

    s.Keq_gmax_psi_[2]=p.Kg_2_;
    s.Keq_gmax_psi_[3]=p.Kg_3_;
    s.Keq_gmax_psi_[4]=p.Kg_4_;
    s.Keq_gmax_psi_[5]=p.Kg_5_;


    s.Keq_gmax_omega_=std::vector<double>(7,0);

    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);



    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    return s;
  }


  myParameters p_;


  // BaseModel interface
public:
  Model051(){}
  ~Model051(){}
  virtual std::string id() const  override
  {
    return "Model 0.51";
  }
  static double number()
  {
    return 0.51;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);

    out.push_back("Kg_2",p_.Kg_2_);
    out.push_back("Kg_3",p_.Kg_3_);
    out.push_back("Kg_4",p_.Kg_4_);
    out.push_back("Kg_5",p_.Kg_5_);


    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );
    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);






    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq_psi");
    p_.kcat_=p.get("kcat_psi");
    p_.Kg_2_=p.get("Kg_2");
    p_.Kg_3_=p.get("Kg_3");
    p_.Kg_4_=p.get("Kg_4");
    p_.Kg_5_=p.get("Kg_5");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");


    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");


  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }





  Model051(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};



class Model10:public BaseModel
{
  SimplestModel m;

  Model10* clone()const override { return new Model10;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;
    s.g_left_[4]=p.g_32_;
    s.g_left_[5]=p.g_32_;
    s.g_left_[6]=p.g_32_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model10(){}
  ~Model10(){}
  virtual std::string id() const  override
  {
    return "Model 1.0";
  }
  static double number()
  {
    return 1;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");


    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model10(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};



class Model100m:public MicrogliaModel
{
  SimplestModel m;

  Model100m* clone()const override { return new Model100m;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;
    double g_M1M2_;
    double g_M2M1_;

    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);
    s.ksig_omega_[1]=p.k_sig_;
    s.ksig_omega_[5]=p.k_sig_;
    s.ksig_omega_[6]=p.k_sig_*1.5;
    s.ksig_omega_[7]=p.k_sig_*3;
    s.ksig_omega_[8]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_psi_=std::vector<double>(9,0);

    s.ksig_max_psi_[1]=p.k_sig_max_;
    s.ksig_max_psi_[5]=p.k_sig_max_;
    s.ksig_max_psi_[6]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[7]=p.k_sig_max_*3;
    s.ksig_max_psi_[8]=p.k_sig_max_*6;

    s.ksig_max_omega_[1]=p.k_sig_max_;
    s.ksig_max_omega_[5]=p.k_sig_max_;
    s.ksig_max_omega_[6]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[7]=p.k_sig_max_*3;
    s.ksig_max_omega_[8]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_32_;
    s.g_left_[8]=p.g_32_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_23_;
    s.g_rigth_[7]=p.g_23_;

    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_;

    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;
    s.g_max_psi_[6]=p.g_max_;
    s.g_max_psi_[7]=p.g_max_;

    s.g_max_omega_[1]=p.g_max_;

    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;
    s.g_max_omega_[6]=p.g_max_;
    s.g_max_omega_[7]=p.g_max_;


    s.a_=std::vector<double> (9,0);
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_2_*p.a_factor_;
    s.a_[7]=s.a_[4]*p.a_factor_;
    s.a_[8]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (9,0);
    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[5]=p.a_max_;
    s.a_psi_[6]=p.a_max_*p.a_factor_;
    s.a_psi_[7]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[8]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[5]=p.a_max_;
    s.a_omega_[6]=p.a_max_*p.a_factor_;
    s.a_omega_[7]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[8]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_N_;
    s.N_[2]=p.N_N_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_0_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_2_*1.5;
    s.N_[7]=p.N_2_*3;
    s.N_[8]=p.N_2_*6;



    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_N_;
    s.M_[2]=p.N_N_;
    s.M_[3]=p.N_0_;
    s.M_[4]=p.N_0_;
    s.M_[5]=p.N_2_;
    s.M_[6]=p.N_2_*1.5;
    s.M_[7]=p.N_2_*3;
    s.M_[8]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;



    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model100m(){}
  ~Model100m(){}
  virtual std::string id() const  override
  {
    return "Model 11.0";
  }
  static double number()
  {
    return 11;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);

    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");
    p_.g_max_=p.get("g_max");

    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model100m(const Parameters& p)
  {
    loadParameters(p);
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


};



class Model111:public BaseModel
{
  SimplestModel m;

  Model111* clone()const override { return new Model111;}

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_12_;
    double g_23_;
    double g_34_;
    double g_45_;

    double g_10_;
    double g_21_;

    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model111(){}
  ~Model111(){}
  virtual std::string id() const  override
  {
    return "Model 1.11";
  }
  static double number()
  {
    return 1.11;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model111(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model112:public BaseModel
{
  SimplestModel m;

  Model112* clone()const override { return new Model112;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model112(){}
  ~Model112(){}
  virtual std::string id() const  override
  {
    return "Model 1.12";
  }
  static double number()
  {
    return 1.12;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model112(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model112_22:public BaseModel
{
  SimplestModel m;

  Model112_22* clone()const override { return new Model112_22;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[2]=p.a_max_1_;
    s.a_omega_[3]=p.a_max_2_;
    s.a_omega_[4]=p.a_max_3_;
    s.a_omega_[5]=p.a_max_4_;
    s.a_omega_[6]=p.a_max_5_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model112_22(){}
  ~Model112_22(){}
  virtual std::string id() const  override
  {
    return "Model 1.12022";
  }
  static double number()
  {
    return 1.12022;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model112_22(const Parameters& p)
  {
    loadParameters(p);
  }


};





class Model112_22_31:public BaseModel
{
  SimplestModel m;

  Model112_22_31* clone()const override { return new Model112_22_31;}


  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[2]=p.a_max_1_;
    s.a_omega_[3]=p.a_max_2_;
    s.a_omega_[4]=p.a_max_3_;
    s.a_omega_[5]=p.a_max_4_;
    s.a_omega_[6]=p.a_max_5_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model112_22_31(){}
  ~Model112_22_31(){}
  virtual std::string id() const  override
  {
    return "Model 1.12022031";
  }
  static double number()
  {
    return 1.12022031;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model112_22_31(const Parameters& p)
  {
    loadParameters(p);
  }


};

class Model112_51:public BaseModel
{
  SimplestModel m;

  Model112_51* clone()const override { return new Model112_51;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double Kg_2_;
    double Kg_3_;
    double Kg_4_;
    double Kg_5_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_psi_[2]=p.Kg_2_;
    s.Keq_gmax_psi_[3]=p.Kg_3_;
    s.Keq_gmax_psi_[4]=p.Kg_4_;
    s.Keq_gmax_psi_[5]=p.Kg_5_;

    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);
    s.Keq_gmax_omega_[2]=p.Kg_2_;
    s.Keq_gmax_omega_[3]=p.Kg_3_;
    s.Keq_gmax_omega_[4]=p.Kg_4_;
    s.Keq_gmax_omega_[5]=p.Kg_5_;



    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model112_51(){}
  ~Model112_51(){}
  virtual std::string id() const  override
  {
    return "Model 1.12051";
  }
  static double number()
  {
    return 1.12051;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("Kg_2",p_.Kg_2_);
    out.push_back("Kg_3",p_.Kg_3_);
    out.push_back("Kg_4",p_.Kg_4_);
    out.push_back("Kg_5",p_.Kg_5_);

    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.Kg_2_=p.get("Kg_2");
    p_.Kg_3_=p.get("Kg_3");
    p_.Kg_4_=p.get("Kg_4");
    p_.Kg_5_=p.get("Kg_5");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model112_51(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model112_52:public BaseModel
{
  SimplestModel m;

  Model112_52* clone()const override { return new Model112_52;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double Kg_psi_2_;
    double Kg_psi_3_;
    double Kg_psi_4_;
    double Kg_psi_5_;

    double Kg_omega_2_;
    double Kg_omega_3_;
    double Kg_omega_4_;
    double Kg_omega_5_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_psi_[2]=p.Kg_psi_2_;
    s.Keq_gmax_psi_[3]=p.Kg_psi_3_;
    s.Keq_gmax_psi_[4]=p.Kg_psi_4_;
    s.Keq_gmax_psi_[5]=p.Kg_psi_5_;

    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);
    s.Keq_gmax_omega_[2]=p.Kg_omega_2_;
    s.Keq_gmax_omega_[3]=p.Kg_omega_3_;
    s.Keq_gmax_omega_[4]=p.Kg_omega_4_;
    s.Keq_gmax_omega_[5]=p.Kg_omega_5_;

    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model112_52(){}
  ~Model112_52(){}
  virtual std::string id() const  override
  {
    return "Model 1.12052";
  }
  static double number()
  {
    return 1.12052;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("Kg_psi_2",p_.Kg_psi_2_);
    out.push_back("Kg_psi_3",p_.Kg_psi_3_);
    out.push_back("Kg_psi_4",p_.Kg_psi_4_);
    out.push_back("Kg_psi_5",p_.Kg_psi_5_);

    out.push_back("Kg_omega_2",p_.Kg_omega_2_);
    out.push_back("Kg_omega_3",p_.Kg_omega_3_);
    out.push_back("Kg_omega_4",p_.Kg_omega_4_);
    out.push_back("Kg_omega_5",p_.Kg_omega_5_);


    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.Kg_psi_2_=p.get("Kg_psi_2");
    p_.Kg_psi_3_=p.get("Kg_psi_3");
    p_.Kg_psi_4_=p.get("Kg_psi_4");
    p_.Kg_psi_5_=p.get("Kg_psi_5");

    p_.Kg_omega_2_=p.get("Kg_omega_2");
    p_.Kg_omega_3_=p.get("Kg_omega_3");
    p_.Kg_omega_4_=p.get("Kg_omega_4");
    p_.Kg_omega_5_=p.get("Kg_omega_5");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model112_52(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model113:public BaseModel
{
  SimplestModel m;

  Model113* clone()const override { return new Model113;}

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_12_;
    double g_23_;
    double g_34_;
    double g_45_;

    double g_10_;
    double g_21_;

    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model113(){}
  ~Model113(){}
  virtual std::string id() const  override
  {
    return "Model 1.13";
  }
  static double number()
  {
    return 1.13;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model113(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model113_42:public BaseModel
{
  SimplestModel m;

  Model113_42* clone()const override { return new Model113_42;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_12_;
    double g_23_;
    double g_34_;
    double g_45_;

    double g_10_;
    double g_21_;

    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_1;
    double k_sig_max_2;
    double k_sig_max_3;
    double k_sig_max_4;
    double k_sig_max_5;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_[2]=p.k_sig_max_1;
    s.ksig_max_omega_[3]=p.k_sig_max_2;
    s.ksig_max_omega_[4]=p.k_sig_max_3;
    s.ksig_max_omega_[5]=p.k_sig_max_4;
    s.ksig_max_omega_[6]=p.k_sig_max_5;

    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[2]=p.k_sig_max_1;
    s.ksig_max_psi_[3]=p.k_sig_max_2;
    s.ksig_max_psi_[4]=p.k_sig_max_3;
    s.ksig_max_psi_[5]=p.k_sig_max_4;
    s.ksig_max_psi_[6]=p.k_sig_max_5;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model113_42(){}
  ~Model113_42(){}
  virtual std::string id() const  override
  {
    return "Model 1.13042";
  }
  static double number()
  {
    return 1.13042;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max_1",p_.k_sig_max_1);
    out.push_back("k_sig_max_2",p_.k_sig_max_2);
    out.push_back("k_sig_max_3",p_.k_sig_max_3);
    out.push_back("k_sig_max_4",p_.k_sig_max_4);
    out.push_back("k_sig_max_5",p_.k_sig_max_5);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_1=p.get("k_sig_max_1");
    p_.k_sig_max_2=p.get("k_sig_max_2");
    p_.k_sig_max_3=p.get("k_sig_max_3");
    p_.k_sig_max_4=p.get("k_sig_max_4");
    p_.k_sig_max_5=p.get("k_sig_max_5");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model113_42(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model114:public BaseModel
{
  SimplestModel m;

  Model114* clone()const override { return new Model114;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model114(){}
  ~Model114(){}
  virtual std::string id() const  override
  {
    return "Model 1.14";
  }
  static double number()
  {
    return 1.14;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model114(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model114_24_44:public BaseModel
{
  SimplestModel m;

  Model114_24_44* clone()const override { return new Model114_24_44;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_psi_Neuron_;
    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_1;
    double k_sig_max_psi_2;
    double k_sig_max_psi_3;
    double k_sig_max_psi_4;
    double k_sig_max_psi_5;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_[2]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[3]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[2]=p.k_sig_max_psi_1;
    s.ksig_max_psi_[3]=p.k_sig_max_psi_2;
    s.ksig_max_psi_[4]=p.k_sig_max_psi_3;
    s.ksig_max_psi_[5]=p.k_sig_max_psi_4;
    s.ksig_max_psi_[6]=p.k_sig_max_psi_5;



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[2]=p.a_max_psi_1_;
    s.a_psi_[3]=p.a_max_psi_2_;
    s.a_psi_[4]=p.a_max_psi_3_;
    s.a_psi_[5]=p.a_max_psi_4_;
    s.a_psi_[6]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (7,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[2]=p.a_max_omega_1_;
    s.a_omega_[3]=p.a_max_omega_2_;
    s.a_omega_[4]=p.a_max_omega_3_;
    s.a_omega_[5]=p.a_max_omega_4_;
    s.a_omega_[6]=p.a_max_omega_5_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model114_24_44(){}
  ~Model114_24_44(){}
  virtual std::string id() const  override
  {
    return "Model 1.14024044";
  }
  static double number()
  {
    return 1.14024044;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_1",p_.k_sig_max_psi_1);
    out.push_back("k_sig_max_psi_2",p_.k_sig_max_psi_2);
    out.push_back("k_sig_max_psi_3",p_.k_sig_max_psi_3);
    out.push_back("k_sig_max_psi_4",p_.k_sig_max_psi_4);
    out.push_back("k_sig_max_psi_5",p_.k_sig_max_psi_5);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_1=p.get("k_sig_max_psi_1");
    p_.k_sig_max_psi_2=p.get("k_sig_max_psi_2");
    p_.k_sig_max_psi_3=p.get("k_sig_max_psi_3");
    p_.k_sig_max_psi_4=p.get("k_sig_max_psi_4");
    p_.k_sig_max_psi_5=p.get("k_sig_max_psi_5");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model114_24_44(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model114_24_32_44:public BaseModel
{
  SimplestModel m;

  Model114_24_32_44* clone()const override { return new Model114_24_32_44;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_psi_Neuron_;
    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_1;
    double k_sig_max_psi_2;
    double k_sig_max_psi_3;
    double k_sig_max_psi_4;
    double k_sig_max_psi_5;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_[2]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[3]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[2]=p.k_sig_max_psi_1;
    s.ksig_max_psi_[3]=p.k_sig_max_psi_2;
    s.ksig_max_psi_[4]=p.k_sig_max_psi_3;
    s.ksig_max_psi_[5]=p.k_sig_max_psi_4;
    s.ksig_max_psi_[6]=p.k_sig_max_psi_5;



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;
    s.g_left_[4]=p.g_32_;
    s.g_left_[5]=p.g_32_;
    s.g_left_[6]=p.g_32_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[2]=p.a_max_psi_1_;
    s.a_psi_[3]=p.a_max_psi_2_;
    s.a_psi_[4]=p.a_max_psi_3_;
    s.a_psi_[5]=p.a_max_psi_4_;
    s.a_psi_[6]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (7,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[2]=p.a_max_omega_1_;
    s.a_omega_[3]=p.a_max_omega_2_;
    s.a_omega_[4]=p.a_max_omega_3_;
    s.a_omega_[5]=p.a_max_omega_4_;
    s.a_omega_[6]=p.a_max_omega_5_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_1_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_3_;
    s.N_[5]=p.N_4_;
    s.N_[6]=p.N_5_;


    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_0_;
    s.M_[2]=p.M_1_;
    s.M_[3]=p.M_2_;
    s.M_[4]=p.M_3_;
    s.M_[5]=p.M_4_;
    s.M_[6]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model114_24_32_44(){}
  ~Model114_24_32_44(){}
  virtual std::string id() const  override
  {
    return "Model 1.14024032044";
  }
  static double number()
  {
    return 1.14024032044;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("M_0",p_.M_0_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_1",p_.k_sig_max_psi_1);
    out.push_back("k_sig_max_psi_2",p_.k_sig_max_psi_2);
    out.push_back("k_sig_max_psi_3",p_.k_sig_max_psi_3);
    out.push_back("k_sig_max_psi_4",p_.k_sig_max_psi_4);
    out.push_back("k_sig_max_psi_5",p_.k_sig_max_psi_5);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_1=p.get("k_sig_max_psi_1");
    p_.k_sig_max_psi_2=p.get("k_sig_max_psi_2");
    p_.k_sig_max_psi_3=p.get("k_sig_max_psi_3");
    p_.k_sig_max_psi_4=p.get("k_sig_max_psi_4");
    p_.k_sig_max_psi_5=p.get("k_sig_max_psi_5");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model114_24_32_44(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model114_24_32_44m:public MicrogliaModel
{
  SimplestModel m;

  Model114_24_32_44m* clone()const override { return new Model114_24_32_44m;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;


    double g_max_psi_M2_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;

    double g_max_omega_M2_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;
    double a_2_;
    double a_factor_;
    double a_max_psi_Neuron_;
    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_omega_M2;
    double k_sig_max_omega_M1;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;
    double k_sig_max_psi_1;
    double k_sig_max_psi_2;
    double k_sig_max_psi_3;
    double k_sig_max_psi_4;
    double k_sig_max_psi_5;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);
    s.ksig_omega_[1]=p.k_sig_;
    s.ksig_omega_[2]=p.k_sig_;
    s.ksig_omega_[5]=p.k_sig_;
    s.ksig_omega_[6]=p.k_sig_*1.5;
    s.ksig_omega_[7]=p.k_sig_*3;
    s.ksig_omega_[8]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_[1]=p.k_sig_max_omega_M2;
    s.ksig_max_omega_[2]=p.k_sig_max_omega_M1;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[7]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[8]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;
    s.ksig_max_psi_[4]=p.k_sig_max_psi_1;
    s.ksig_max_psi_[5]=p.k_sig_max_psi_2;
    s.ksig_max_psi_[6]=p.k_sig_max_psi_3;
    s.ksig_max_psi_[7]=p.k_sig_max_psi_4;
    s.ksig_max_psi_[8]=p.k_sig_max_psi_5;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_32_;
    s.g_left_[8]=p.g_32_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_23_;
    s.g_rigth_[7]=p.g_23_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M2_;

    s.g_max_psi_[4]=p.g_max_psi_2_;
    s.g_max_psi_[5]=p.g_max_psi_3_;
    s.g_max_psi_[6]=p.g_max_psi_4_;
    s.g_max_psi_[7]=p.g_max_psi_5_;

    s.g_max_omega_[1]=p.g_max_omega_M2_;

    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_2_*p.a_factor_;
    s.a_[7]=s.a_[4]*p.a_factor_;
    s.a_[8]=s.a_[5]*p.a_factor_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[4]=p.a_max_psi_1_;
    s.a_psi_[5]=p.a_max_psi_2_;
    s.a_psi_[6]=p.a_max_psi_3_;
    s.a_psi_[7]=p.a_max_psi_4_;
    s.a_psi_[8]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model114_24_32_44m(){}
  ~Model114_24_32_44m(){}
  virtual std::string id() const  override
  {
    return "Model 11.14024032044m";
  }
  static double number()
  {
    return 11.14024032044;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );

    out.push_back("g_max_psi_M2",p_.g_max_psi_M2_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );

    out.push_back("g_max_omega_M2",p_.g_max_omega_M2_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_0",p_.N_M2_ );
    out.push_back("N_0",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_0",p_.M_M2_ );
    out.push_back("M_0",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max_omega_M2",p_.k_sig_max_omega_M2);
    out.push_back("k_sig_max_omega_M1",p_.k_sig_max_omega_M1);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);
    out.push_back("k_sig_max_psi_1",p_.k_sig_max_psi_1);
    out.push_back("k_sig_max_psi_2",p_.k_sig_max_psi_2);
    out.push_back("k_sig_max_psi_3",p_.k_sig_max_psi_3);
    out.push_back("k_sig_max_psi_4",p_.k_sig_max_psi_4);
    out.push_back("k_sig_max_psi_5",p_.k_sig_max_psi_5);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }
  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");

    p_.g_max_psi_M2_=p.get("g_max_psi_M2");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_M2_=p.get("g_max_omega_M2");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_omega_M2=p.get("k_sig_max_omega_M2");
    p_.k_sig_max_omega_M1=p.get("k_sig_max_omega_M1");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");
    p_.k_sig_max_psi_1=p.get("k_sig_max_psi_1");
    p_.k_sig_max_psi_2=p.get("k_sig_max_psi_2");
    p_.k_sig_max_psi_3=p.get("k_sig_max_psi_3");
    p_.k_sig_max_psi_4=p.get("k_sig_max_psi_4");
    p_.k_sig_max_psi_5=p.get("k_sig_max_psi_5");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model114_24_32_44m(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model_01:public MicrogliaModel
{
  SimplestModel m;

  Model_01* clone()const override { return new Model_01;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_M2;
    double k_sig_M1;



    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);

    s.ksig_omega_[1]=p.k_sig_M2;
    s.ksig_omega_[2]=p.k_sig_M1;



    s.ksig_max_omega_=std::vector<double>(9,0);

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;



    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_01(){}
  ~Model_01(){}
  virtual std::string id() const  override
  {
    return "Model_01";
  }
  static double number()
  {
    return 100;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );

    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_M2",p_.k_sig_M2);
    out.push_back("k_sig_M1",p_.k_sig_M1);


    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }
  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_M2=p.get("k_sig_M2");
    p_.k_sig_M1=p.get("k_sig_M1");



    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model_01(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model_11:public MicrogliaModel
{
  SimplestModel m;

  Model_11* clone()const override { return new Model_11;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;
    double g_max_psi_1_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_psi_0_;

    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_M2;
    double k_sig_M1;



    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);

    s.ksig_omega_[1]=p.k_sig_M2;
    s.ksig_omega_[2]=p.k_sig_M1;



    s.ksig_max_omega_=std::vector<double>(9,0);

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;

    s.g_max_psi_[1]=p.g_max_psi_1_;
    s.g_max_psi_[4]=p.g_max_psi_2_;
    s.g_max_psi_[5]=p.g_max_psi_3_;
    s.g_max_psi_[6]=p.g_max_psi_4_;
    s.g_max_psi_[7]=p.g_max_psi_5_;


    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[4]=p.a_max_psi_1_;
    s.a_psi_[5]=p.a_max_psi_2_;
    s.a_psi_[6]=p.a_max_psi_3_;
    s.a_psi_[7]=p.a_max_psi_4_;
    s.a_psi_[8]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_11(){}
  ~Model_11(){}
  virtual std::string id() const  override
  {
    return "Model 11.14024032044m";
  }
  static double number()
  {
    return 1100;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );
    out.push_back("g_max_psi_1",p_.g_max_psi_1_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );

    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_0",p_.a_max_psi_0_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_M2",p_.k_sig_M2);
    out.push_back("k_sig_M1",p_.k_sig_M1);


    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }
  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_psi_1_=p.get("g_max_psi_1");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_0_=p.get("a_max_psi_0");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_M2=p.get("k_sig_M2");
    p_.k_sig_M1=p.get("k_sig_M1");



    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model_11(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model_02:public MicrogliaModel
{
  SimplestModel m;

  Model_02* clone()const override { return new Model_02;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;

    double g_max_omega_M1_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_M2;
    double k_sig_M1;



    double k_sig_max_omega_M2;
    double k_sig_max_omega_M1;

    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);

    s.ksig_omega_[1]=p.k_sig_M2;
    s.ksig_omega_[2]=p.k_sig_M1;



    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_[1]=p.k_sig_max_omega_M2;
    s.ksig_max_omega_[2]=p.k_sig_max_omega_M1;

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;


    s.g_max_omega_[1]=p.g_max_omega_M1_;

    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_02(){}
  ~Model_02(){}
  virtual std::string id() const  override
  {
    return "Model_02";
  }
  static double number()
  {
    return 200;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );

    out.push_back("g_max_omega_M1",p_.g_max_omega_M1_ );
    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_M2",p_.k_sig_M2);
    out.push_back("k_sig_M1",p_.k_sig_M1);

    out.push_back("k_sig_max_omega_M2",p_.k_sig_max_omega_M2);
    out.push_back("k_sig_max_omega_M1",p_.k_sig_max_omega_M1);

    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_omega_M1_=p.get("g_max_omega_M1");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_M2=p.get("k_sig_M2");
    p_.k_sig_M1=p.get("k_sig_M1");


    p_.k_sig_max_omega_M2=p.get("k_sig_max_omega_M2");
    p_.k_sig_max_omega_M1=p.get("k_sig_max_omega_M1");

    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model_02(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model_12:public MicrogliaModel
{
  SimplestModel m;

  Model_12* clone()const override { return new Model_12;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;
    double g_max_psi_1_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;

    double g_max_omega_M1_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_psi_0_;

    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_M2;
    double k_sig_M1;



    double k_sig_max_omega_M2;
    double k_sig_max_omega_M1;

    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);

    s.ksig_omega_[1]=p.k_sig_M2;
    s.ksig_omega_[2]=p.k_sig_M1;



    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_[1]=p.k_sig_max_omega_M2;
    s.ksig_max_omega_[2]=p.k_sig_max_omega_M1;

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;

    s.g_max_psi_[1]=p.g_max_psi_1_;
    s.g_max_psi_[4]=p.g_max_psi_2_;
    s.g_max_psi_[5]=p.g_max_psi_3_;
    s.g_max_psi_[6]=p.g_max_psi_4_;
    s.g_max_psi_[7]=p.g_max_psi_5_;

    s.g_max_omega_[1]=p.g_max_omega_M1_;

    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[4]=p.a_max_psi_1_;
    s.a_psi_[5]=p.a_max_psi_2_;
    s.a_psi_[6]=p.a_max_psi_3_;
    s.a_psi_[7]=p.a_max_psi_4_;
    s.a_psi_[8]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_12(){}
  ~Model_12(){}
  virtual std::string id() const  override
  {
    return "Model_12";
  }
  static double number()
  {
    return 1200;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );
    out.push_back("g_max_psi_1",p_.g_max_psi_1_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );

    out.push_back("g_max_omega_M1",p_.g_max_omega_M1_ );
    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_0",p_.a_max_psi_0_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_M2",p_.k_sig_M2);
    out.push_back("k_sig_M1",p_.k_sig_M1);

    out.push_back("k_sig_max_omega_M2",p_.k_sig_max_omega_M2);
    out.push_back("k_sig_max_omega_M1",p_.k_sig_max_omega_M1);

    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }
  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_psi_1_=p.get("g_max_psi_1");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_M1_=p.get("g_max_omega_M1");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_0_=p.get("a_max_psi_0");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_M2=p.get("k_sig_M2");
    p_.k_sig_M1=p.get("k_sig_M1");


    p_.k_sig_max_omega_M2=p.get("k_sig_max_omega_M2");
    p_.k_sig_max_omega_M1=p.get("k_sig_max_omega_M1");

    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model_12(const Parameters& p)
  {
    loadParameters(p);
  }


};

class Model_101:public MicrogliaModel
{
  SimplestModel m;

  Model_101* clone()const override { return new Model_101;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;
    double g_max_psi_1_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_psi_0_;

    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_0;
    double k_sig_1;
    double k_sig_2;
    double k_sig_3;
    double k_sig_4;
    double k_sig_5;



    double k_sig_max_omega_0;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_0;
    double k_sig_max_psi_1;
    double k_sig_max_psi_2;
    double k_sig_max_psi_3;
    double k_sig_max_psi_4;
    double k_sig_max_psi_5;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);



    s.ksig_omega_[3]=p.k_sig_0;
    s.ksig_omega_[4]=p.k_sig_1;
    s.ksig_omega_[5]=p.k_sig_2;
    s.ksig_omega_[6]=p.k_sig_3;
    s.ksig_omega_[7]=p.k_sig_4;
    s.ksig_omega_[8]=p.k_sig_5;

    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_[3]=p.k_sig_max_omega_0;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[7]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[8]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[3]=p.k_sig_max_psi_0;
    s.ksig_max_psi_[4]=p.k_sig_max_psi_1;
    s.ksig_max_psi_[5]=p.k_sig_max_psi_2;
    s.ksig_max_psi_[6]=p.k_sig_max_psi_3;
    s.ksig_max_psi_[7]=p.k_sig_max_psi_4;
    s.ksig_max_psi_[8]=p.k_sig_max_psi_5;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;

    s.g_max_psi_[1]=p.g_max_psi_1_;
    s.g_max_psi_[4]=p.g_max_psi_2_;
    s.g_max_psi_[5]=p.g_max_psi_3_;
    s.g_max_psi_[6]=p.g_max_psi_4_;
    s.g_max_psi_[7]=p.g_max_psi_5_;


    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[4]=p.a_max_psi_1_;
    s.a_psi_[5]=p.a_max_psi_2_;
    s.a_psi_[6]=p.a_max_psi_3_;
    s.a_psi_[7]=p.a_max_psi_4_;
    s.a_psi_[8]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_101(){}
  ~Model_101(){}
  virtual std::string id() const  override
  {
    return "Model_101";
  }
  static double number()
  {
    return 1010;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );
    out.push_back("g_max_psi_1",p_.g_max_psi_1_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );

    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_0",p_.a_max_psi_0_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_0",p_.k_sig_0);
    out.push_back("k_sig_1",p_.k_sig_1);
    out.push_back("k_sig_2",p_.k_sig_2);
    out.push_back("k_sig_3",p_.k_sig_3);
    out.push_back("k_sig_4",p_.k_sig_4);
    out.push_back("k_sig_5",p_.k_sig_5);

    out.push_back("k_sig_max_omega_0",p_.k_sig_max_omega_0);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_0",p_.k_sig_max_psi_0);
    out.push_back("k_sig_max_psi_1",p_.k_sig_max_psi_1);
    out.push_back("k_sig_max_psi_2",p_.k_sig_max_psi_2);
    out.push_back("k_sig_max_psi_3",p_.k_sig_max_psi_3);
    out.push_back("k_sig_max_psi_4",p_.k_sig_max_psi_4);
    out.push_back("k_sig_max_psi_5",p_.k_sig_max_psi_5);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_psi_1_=p.get("g_max_psi_1");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_0_=p.get("a_max_psi_0");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_0=p.get("k_sig_0");
    p_.k_sig_1=p.get("k_sig_1");
    p_.k_sig_2=p.get("k_sig_2");
    p_.k_sig_3=p.get("k_sig_3");
    p_.k_sig_4=p.get("k_sig_4");
    p_.k_sig_5=p.get("k_sig_5");


    p_.k_sig_max_omega_0=p.get("k_sig_max_omega_0");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_0=p.get("k_sig_max_psi_0");
    p_.k_sig_max_psi_1=p.get("k_sig_max_psi_1");
    p_.k_sig_max_psi_2=p.get("k_sig_max_psi_2");
    p_.k_sig_max_psi_3=p.get("k_sig_max_psi_3");
    p_.k_sig_max_psi_4=p.get("k_sig_max_psi_4");
    p_.k_sig_max_psi_5=p.get("k_sig_max_psi_5");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model_101(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model_10:public MicrogliaModel
{
  SimplestModel m;

  Model_10* clone()const override { return new Model_10;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double kcat_psi;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;
    double g_max_psi_1_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;


    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;


    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_psi_0_;

    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double inj_width_;

    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;



    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);


    s.ksig_omega_=std::vector<double>(9,0);




    s.ksig_max_omega_=std::vector<double>(9,0);

    s.ksig_max_psi_=std::vector<double>(9,0);

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;

    s.g_max_psi_[1]=p.g_max_psi_1_;
    s.g_max_psi_[4]=p.g_max_psi_2_;
    s.g_max_psi_[5]=p.g_max_psi_3_;
    s.g_max_psi_[6]=p.g_max_psi_4_;
    s.g_max_psi_[7]=p.g_max_psi_5_;




    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[4]=p.a_max_psi_1_;
    s.a_psi_[5]=p.a_max_psi_2_;
    s.a_psi_[6]=p.a_max_psi_3_;
    s.a_psi_[7]=p.a_max_psi_4_;
    s.a_psi_[8]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (9,0);




    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;


    s.M_=std::vector<double> (9,0);






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_10(){}
  ~Model_10(){}
  virtual std::string id() const  override
  {
    return "Model_10";
  }
  static double number()
  {
    return 1000;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );
    out.push_back("g_max_psi_1",p_.g_max_psi_1_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );

    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_0",p_.a_max_psi_0_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);




    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.kcat_psi=p.get("kcat_psi");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_psi_1_=p.get("g_max_psi_1");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");

    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_0_=p.get("a_max_psi_0");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");



    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model_10(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model_03:public MicrogliaModel
{
  SimplestModel m;

  Model_03* clone()const override { return new Model_03;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;

    double g_max_omega_M1_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_M2;
    double k_sig_M1;
    double k_sig_0;
    double k_sig_1;
    double k_sig_2;
    double k_sig_3;
    double k_sig_4;
    double k_sig_5;



    double k_sig_max_omega_M2;
    double k_sig_max_omega_M1;
    double k_sig_max_omega_0;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);

    s.ksig_omega_[1]=p.k_sig_M2;
    s.ksig_omega_[2]=p.k_sig_M1;


    s.ksig_omega_[3]=p.k_sig_0;
    s.ksig_omega_[4]=p.k_sig_1;
    s.ksig_omega_[5]=p.k_sig_2;
    s.ksig_omega_[6]=p.k_sig_3;
    s.ksig_omega_[7]=p.k_sig_4;
    s.ksig_omega_[8]=p.k_sig_5;

    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_[1]=p.k_sig_max_omega_M2;
    s.ksig_max_omega_[2]=p.k_sig_max_omega_M1;
    s.ksig_max_omega_[3]=p.k_sig_max_omega_0;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[7]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[8]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;


    s.g_max_omega_[1]=p.g_max_omega_M1_;

    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_03(){}
  ~Model_03(){}
  virtual std::string id() const  override
  {
    return "Model_03";
  }
  static double number()
  {
    return 300;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );

    out.push_back("g_max_omega_M1",p_.g_max_omega_M1_ );
    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_M2",p_.k_sig_M2);
    out.push_back("k_sig_M1",p_.k_sig_M1);
    out.push_back("k_sig_0",p_.k_sig_0);
    out.push_back("k_sig_1",p_.k_sig_1);
    out.push_back("k_sig_2",p_.k_sig_2);
    out.push_back("k_sig_3",p_.k_sig_3);
    out.push_back("k_sig_4",p_.k_sig_4);
    out.push_back("k_sig_5",p_.k_sig_5);

    out.push_back("k_sig_max_omega_M2",p_.k_sig_max_omega_M2);
    out.push_back("k_sig_max_omega_M1",p_.k_sig_max_omega_M1);
    out.push_back("k_sig_max_omega_0",p_.k_sig_max_omega_0);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_omega_M1_=p.get("g_max_omega_M1");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_M2=p.get("k_sig_M2");
    p_.k_sig_M1=p.get("k_sig_M1");
    p_.k_sig_0=p.get("k_sig_0");
    p_.k_sig_1=p.get("k_sig_1");
    p_.k_sig_2=p.get("k_sig_2");
    p_.k_sig_3=p.get("k_sig_3");
    p_.k_sig_4=p.get("k_sig_4");
    p_.k_sig_5=p.get("k_sig_5");


    p_.k_sig_max_omega_M2=p.get("k_sig_max_omega_M2");
    p_.k_sig_max_omega_M1=p.get("k_sig_max_omega_M1");
    p_.k_sig_max_omega_0=p.get("k_sig_max_omega_0");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model_03(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model_13:public MicrogliaModel
{
  SimplestModel m;

  Model_13* clone()const override { return new Model_13;}

  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_M1M2_;
    double g_M2M1_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_32_;

    double g_34_;
    double g_43_;

    double g_45_;
    double g_54_;


    double g_max_psi_M1_;
    double g_max_psi_1_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;

    double g_max_omega_M1_;

    double g_max_omega_1_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;

    double N_M2_;
    double N_M1_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;

    double M_M2_;
    double M_M1_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;

    double N_Astr_;
    double N_Neuron_;
    double N_Microglia_;

    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;

    double a_max_psi_Neuron_;
    double a_max_psi_0_;

    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_0_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;


    double k_sig_M2;
    double k_sig_M1;
    double k_sig_0;
    double k_sig_1;
    double k_sig_2;
    double k_sig_3;
    double k_sig_4;
    double k_sig_5;



    double k_sig_max_omega_M2;
    double k_sig_max_omega_M1;
    double k_sig_max_omega_0;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_M2;
    double k_sig_max_psi_M1;
    double k_sig_max_psi_0;
    double k_sig_max_psi_1;
    double k_sig_max_psi_2;
    double k_sig_max_psi_3;
    double k_sig_max_psi_4;
    double k_sig_max_psi_5;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(9,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(9,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(9,0);

    s.ksig_omega_[1]=p.k_sig_M2;
    s.ksig_omega_[2]=p.k_sig_M1;


    s.ksig_omega_[3]=p.k_sig_0;
    s.ksig_omega_[4]=p.k_sig_1;
    s.ksig_omega_[5]=p.k_sig_2;
    s.ksig_omega_[6]=p.k_sig_3;
    s.ksig_omega_[7]=p.k_sig_4;
    s.ksig_omega_[8]=p.k_sig_5;

    s.ksig_max_omega_=std::vector<double>(9,0);
    s.ksig_max_omega_[1]=p.k_sig_max_omega_M2;
    s.ksig_max_omega_[2]=p.k_sig_max_omega_M1;
    s.ksig_max_omega_[3]=p.k_sig_max_omega_0;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[7]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[8]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(9,0);
    s.ksig_max_psi_[1]=p.k_sig_max_psi_M2;
    s.ksig_max_psi_[2]=p.k_sig_max_psi_M1;
    s.ksig_max_psi_[3]=p.k_sig_max_psi_0;
    s.ksig_max_psi_[4]=p.k_sig_max_psi_1;
    s.ksig_max_psi_[5]=p.k_sig_max_psi_2;
    s.ksig_max_psi_[6]=p.k_sig_max_psi_3;
    s.ksig_max_psi_[7]=p.k_sig_max_psi_4;
    s.ksig_max_psi_[8]=p.k_sig_max_psi_5;

    s.g_left_=std::vector<double> (9,0);
    s.g_left_[2]=p.g_M1M2_;

    s.g_left_[4]=p.g_10_;
    s.g_left_[5]=p.g_21_;
    s.g_left_[6]=p.g_32_;
    s.g_left_[7]=p.g_43_;
    s.g_left_[8]=p.g_54_;



    s.g_rigth_=std::vector<double> (9,0);
    s.g_rigth_[1]=p.g_M2M1_;

    s.g_rigth_[3]=p.g_01_;
    s.g_rigth_[4]=p.g_12_;
    s.g_rigth_[5]=p.g_23_;
    s.g_rigth_[6]=p.g_34_;
    s.g_rigth_[7]=p.g_45_;



    s.g_max_omega_=std::vector<double> (9,0);


    s.g_max_psi_=std::vector<double> (9,0);

    s.g_max_psi_[1]=p.g_max_psi_M1_;

    s.g_max_psi_[1]=p.g_max_psi_1_;
    s.g_max_psi_[4]=p.g_max_psi_2_;
    s.g_max_psi_[5]=p.g_max_psi_3_;
    s.g_max_psi_[6]=p.g_max_psi_4_;
    s.g_max_psi_[7]=p.g_max_psi_5_;

    s.g_max_omega_[1]=p.g_max_omega_M1_;

    s.g_max_omega_[3]=p.g_max_omega_1_;
    s.g_max_omega_[4]=p.g_max_omega_2_;
    s.g_max_omega_[5]=p.g_max_omega_3_;
    s.g_max_omega_[6]=p.g_max_omega_4_;
    s.g_max_omega_[7]=p.g_max_omega_5_;


    s.a_=std::vector<double> (9,0);
    s.a_[4]=p.a_1_;
    s.a_[5]=p.a_2_;
    s.a_[6]=p.a_3_;
    s.a_[7]=p.a_4_;
    s.a_[8]=p.a_5_;

    s.a_psi_=std::vector<double> (9,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[4]=p.a_max_psi_1_;
    s.a_psi_[5]=p.a_max_psi_2_;
    s.a_psi_[6]=p.a_max_psi_3_;
    s.a_psi_[7]=p.a_max_psi_4_;
    s.a_psi_[8]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (9,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[4]=p.a_max_omega_1_;
    s.a_omega_[5]=p.a_max_omega_2_;
    s.a_omega_[6]=p.a_max_omega_3_;
    s.a_omega_[7]=p.a_max_omega_4_;
    s.a_omega_[8]=p.a_max_omega_5_;



    s.N_=std::vector<double> (9,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_M2_;
    s.N_[2]=p.N_M1_;
    s.N_[3]=p.N_0_;
    s.N_[4]=p.N_1_;
    s.N_[5]=p.N_2_;
    s.N_[6]=p.N_3_;
    s.N_[7]=p.N_4_;
    s.N_[8]=p.N_5_;


    s.M_=std::vector<double> (9,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_M2_;
    s.M_[2]=p.M_M1_;
    s.M_[3]=p.M_0_;
    s.M_[4]=p.M_1_;
    s.M_[5]=p.M_2_;
    s.M_[6]=p.M_3_;
    s.M_[7]=p.M_4_;
    s.M_[8]=p.M_5_;






    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;

    s.dens_Microglia_=p.N_Microglia_;

    return s;
  }




  myParameters p_;

  // BaseModel interface
public:
  Model_13(){}
  ~Model_13(){}
  virtual std::string id() const  override
  {
    return "Model_13";
  }
  static double number()
  {
    return 1300;
  }

  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_M2M1",p_.g_M2M1_);
    out.push_back("g_M1M2",p_.g_M1M2_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_32",p_.g_32_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_43",p_.g_43_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_54",p_.g_54_ );

    out.push_back("g_max_psi_M1",p_.g_max_psi_M1_ );
    out.push_back("g_max_psi_1",p_.g_max_psi_1_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );

    out.push_back("g_max_omega_M1",p_.g_max_omega_M1_ );
    out.push_back("g_max_omega_1",p_.g_max_omega_1_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );

    out.push_back("M_0",p_.M_0_ );
    out.push_back("N_M2_",p_.N_M2_ );
    out.push_back("N_M1_",p_.N_M1_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );

    out.push_back("M_M2",p_.M_M2_ );
    out.push_back("M_M1",p_.M_M1_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron",p_.N_Neuron_);
    out.push_back("N_Microglia",p_.N_Microglia_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_0",p_.a_max_psi_0_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_0",p_.a_max_omega_0_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("k_sig_M2",p_.k_sig_M2);
    out.push_back("k_sig_M1",p_.k_sig_M1);
    out.push_back("k_sig_0",p_.k_sig_0);
    out.push_back("k_sig_1",p_.k_sig_1);
    out.push_back("k_sig_2",p_.k_sig_2);
    out.push_back("k_sig_3",p_.k_sig_3);
    out.push_back("k_sig_4",p_.k_sig_4);
    out.push_back("k_sig_5",p_.k_sig_5);

    out.push_back("k_sig_max_omega_M2",p_.k_sig_max_omega_M2);
    out.push_back("k_sig_max_omega_M1",p_.k_sig_max_omega_M1);
    out.push_back("k_sig_max_omega_0",p_.k_sig_max_omega_0);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_M2",p_.k_sig_max_psi_M2);
    out.push_back("k_sig_max_psi_M1",p_.k_sig_max_psi_M1);
    out.push_back("k_sig_max_psi_0",p_.k_sig_max_psi_0);
    out.push_back("k_sig_max_psi_1",p_.k_sig_max_psi_1);
    out.push_back("k_sig_max_psi_2",p_.k_sig_max_psi_2);
    out.push_back("k_sig_max_psi_3",p_.k_sig_max_psi_3);
    out.push_back("k_sig_max_psi_4",p_.k_sig_max_psi_4);
    out.push_back("k_sig_max_psi_5",p_.k_sig_max_psi_5);



    out.push_back("inj_width_3",p_.inj_width_3_);
    out.push_back("inj_width_7",p_.inj_width_7_);







    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_M2M1_=p.get("g_M2M1");
    p_.g_M1M2_=p.get("g_M1M2");

    p_.g_01_=p.get("g_01");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");
    p_.g_23_=p.get("g_23");
    p_.g_32_=p.get("g_32");
    p_.g_34_=p.get("g_34");
    p_.g_43_=p.get("g_43");
    p_.g_45_=p.get("g_45");
    p_.g_54_=p.get("g_54");

    p_.g_max_psi_M1_=p.get("g_max_psi_M1");
    p_.g_max_psi_1_=p.get("g_max_psi_1");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_M1_=p.get("g_max_omega_M1");
    p_.g_max_omega_1_=p.get("g_max_omega_1");

    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_M2_=p.get("N_M2");
    p_.N_M1_=p.get("N_M1");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");

    p_.M_M2_=p.get("M_M2");
    p_.M_M1_=p.get("M_M1");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.N_Microglia_=p.get("N_Microglia");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");

    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_0_=p.get("a_max_psi_0");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_0_=p.get("a_max_omega_0");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_M2=p.get("k_sig_M2");
    p_.k_sig_M1=p.get("k_sig_M1");
    p_.k_sig_0=p.get("k_sig_0");
    p_.k_sig_1=p.get("k_sig_1");
    p_.k_sig_2=p.get("k_sig_2");
    p_.k_sig_3=p.get("k_sig_3");
    p_.k_sig_4=p.get("k_sig_4");
    p_.k_sig_5=p.get("k_sig_5");


    p_.k_sig_max_omega_M2=p.get("k_sig_max_omega_M2");
    p_.k_sig_max_omega_M1=p.get("k_sig_max_omega_M1");
    p_.k_sig_max_omega_0=p.get("k_sig_max_omega_0");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_M2=p.get("k_sig_max_psi_M2");
    p_.k_sig_max_psi_M1=p.get("k_sig_max_psi_M1");
    p_.k_sig_max_psi_0=p.get("k_sig_max_psi_0");
    p_.k_sig_max_psi_1=p.get("k_sig_max_psi_1");
    p_.k_sig_max_psi_2=p.get("k_sig_max_psi_2");
    p_.k_sig_max_psi_3=p.get("k_sig_max_psi_3");
    p_.k_sig_max_psi_4=p.get("k_sig_max_psi_4");
    p_.k_sig_max_psi_5=p.get("k_sig_max_psi_5");

    p_.inj_width_3_=p.get("inj_width_3");
    p_.inj_width_7_=p.get("inj_width_7");
  }


  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),
                      toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model_13(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model115:public BaseModel
{
  SimplestModel m;

  Model115* clone()const override { return new Model115;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_34_;
    double g_45_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model115(){}
  ~Model115(){}
  virtual std::string id() const  override
  {
    return "Model 1.15";
  }
  static double number()
  {
    return 1.15;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model115(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model115_22:public BaseModel
{
  SimplestModel m;

  Model115_22* clone()const override { return new Model115_22;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_34_;
    double g_45_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;



    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[2]=p.a_max_1_;
    s.a_omega_[3]=p.a_max_2_;
    s.a_omega_[4]=p.a_max_3_;
    s.a_omega_[5]=p.a_max_4_;
    s.a_omega_[6]=p.a_max_5_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model115_22(){}
  ~Model115_22(){}
  virtual std::string id() const  override
  {
    return "Model 1.15022";
  }
  static double number()
  {
    return 1.15022;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model115_22(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model115_25:public BaseModel
{
  SimplestModel m;

  Model115_25* clone()const override { return new Model115_25;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_34_;
    double g_45_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_psi_Neuron_;
    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[2]=p.a_max_psi_1_;
    s.a_psi_[3]=p.a_max_psi_2_;
    s.a_psi_[4]=p.a_max_psi_3_;
    s.a_psi_[5]=p.a_max_psi_4_;
    s.a_psi_[6]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (7,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[2]=p.a_max_omega_1_;
    s.a_omega_[3]=p.a_max_omega_2_;
    s.a_omega_[4]=p.a_max_omega_3_;
    s.a_omega_[5]=p.a_max_omega_4_;
    s.a_omega_[6]=p.a_max_omega_5_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model115_25(){}
  ~Model115_25(){}
  virtual std::string id() const  override
  {
    return "Model 1.15025";
  }
  static double number()
  {
    return 1.15025;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model115_25(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model121:public BaseModel
{
  SimplestModel m;

  Model121* clone()const override { return new Model121;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_Neuron_;
    double a_max_;
    double a_factor_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;


    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model121(){}
  ~Model121(){}
  virtual std::string id() const  override
  {
    return "Model 1.21";
  }
  static double number()
  {
    return 1.21;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );
    out.push_back("a_factor",p_.a_factor_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.a_factor_=p.get("a_factor");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model121(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model122:public BaseModel
{
  SimplestModel m;

  Model122* clone()const override { return new Model122;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;


    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[2]=p.a_max_1_;
    s.a_omega_[3]=p.a_max_2_;
    s.a_omega_[4]=p.a_max_3_;
    s.a_omega_[5]=p.a_max_4_;
    s.a_omega_[6]=p.a_max_5_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model122(){}
  ~Model122(){}
  virtual std::string id() const  override
  {
    return "Model 1.22";
  }
  static double number()
  {
    return 1.22;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model122(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model123:public BaseModel
{
  SimplestModel m;

  Model123* clone()const override { return new Model123;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_Neuron_;
    double a_max_1_;
    double a_max_2_;
    double a_max_3_;
    double a_max_4_;
    double a_max_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;


    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[2]=p.a_max_1_;
    s.a_psi_[3]=p.a_max_2_;
    s.a_psi_[4]=p.a_max_3_;
    s.a_psi_[5]=p.a_max_4_;
    s.a_psi_[6]=p.a_max_5_;


    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[2]=p.a_max_1_;
    s.a_omega_[3]=p.a_max_2_;
    s.a_omega_[4]=p.a_max_3_;
    s.a_omega_[5]=p.a_max_4_;
    s.a_omega_[6]=p.a_max_5_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model123(){}
  ~Model123(){}
  virtual std::string id() const  override
  {
    return "Model 1.23";
  }
  static double number()
  {
    return 1.23;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max_1",p_.a_max_1_ );
    out.push_back("a_max_2",p_.a_max_2_ );
    out.push_back("a_max_3",p_.a_max_3_ );
    out.push_back("a_max_4",p_.a_max_4_ );
    out.push_back("a_max_5",p_.a_max_5_ );

    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_1_=p.get("a_max_1");
    p_.a_max_2_=p.get("a_max_2");
    p_.a_max_3_=p.get("a_max_3");
    p_.a_max_4_=p.get("a_max_4");
    p_.a_max_5_=p.get("a_max_5");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model123(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model124:public BaseModel
{
  SimplestModel m;

  Model124* clone()const override { return new Model124;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_psi_Neuron_;
    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[2]=p.a_max_psi_1_;
    s.a_psi_[3]=p.a_max_psi_2_;
    s.a_psi_[4]=p.a_max_psi_3_;
    s.a_psi_[5]=p.a_max_psi_4_;
    s.a_psi_[6]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (7,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[2]=p.a_max_omega_1_;
    s.a_omega_[3]=p.a_max_omega_2_;
    s.a_omega_[4]=p.a_max_omega_3_;
    s.a_omega_[5]=p.a_max_omega_4_;
    s.a_omega_[6]=p.a_max_omega_5_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model124(){}
  ~Model124(){}
  virtual std::string id() const  override
  {
    return "Model 1.24";
  }
  static double number()
  {
    return 1.24;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");


    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model124(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model125:public BaseModel
{
  SimplestModel m;

  Model125* clone()const override { return new Model125;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_1_;
    double a_2_;
    double a_3_;
    double a_4_;
    double a_5_;
    double a_max_psi_Neuron_;
    double a_max_psi_1_;
    double a_max_psi_2_;
    double a_max_psi_3_;
    double a_max_psi_4_;
    double a_max_psi_5_;
    double a_max_omega_Neuron_;
    double a_max_omega_1_;
    double a_max_omega_2_;
    double a_max_omega_3_;
    double a_max_omega_4_;
    double a_max_omega_5_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[2]=p.a_1_;
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_3_;
    s.a_[5]=p.a_4_;
    s.a_[6]=p.a_5_;

    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_psi_Neuron_;
    s.a_psi_[2]=p.a_max_psi_1_;
    s.a_psi_[3]=p.a_max_psi_2_;
    s.a_psi_[4]=p.a_max_psi_3_;
    s.a_psi_[5]=p.a_max_psi_4_;
    s.a_psi_[6]=p.a_max_psi_5_;

    s.a_omega_=std::vector<double> (7,0);

    s.a_omega_[0]=p.a_max_omega_Neuron_;
    s.a_omega_[2]=p.a_max_omega_1_;
    s.a_omega_[3]=p.a_max_omega_2_;
    s.a_omega_[4]=p.a_max_omega_3_;
    s.a_omega_[5]=p.a_max_omega_4_;
    s.a_omega_[6]=p.a_max_omega_5_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model125(){}
  ~Model125(){}
  virtual std::string id() const  override
  {
    return "Model 1.25";
  }
  static double number()
  {
    return 1.25;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_1",p_.a_1_ );
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_3",p_.a_3_ );
    out.push_back("a_4",p_.a_4_ );
    out.push_back("a_5",p_.a_5_ );

    out.push_back("a_max_psi_Neuron",p_.a_max_psi_Neuron_ );
    out.push_back("a_max_psi_1",p_.a_max_psi_1_ );
    out.push_back("a_max_psi_2",p_.a_max_psi_2_ );
    out.push_back("a_max_psi_3",p_.a_max_psi_3_ );
    out.push_back("a_max_psi_4",p_.a_max_psi_4_ );
    out.push_back("a_max_psi_5",p_.a_max_psi_5_ );

    out.push_back("a_max_omega_Neuron",p_.a_max_omega_Neuron_ );
    out.push_back("a_max_omega_1",p_.a_max_omega_1_ );
    out.push_back("a_max_omega_2",p_.a_max_omega_2_ );
    out.push_back("a_max_omega_3",p_.a_max_omega_3_ );
    out.push_back("a_max_omega_4",p_.a_max_omega_4_ );
    out.push_back("a_max_omega_5",p_.a_max_omega_5_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_1_=p.get("a_1");
    p_.a_2_=p.get("a_2");
    p_.a_3_=p.get("a_3");
    p_.a_4_=p.get("a_4");
    p_.a_5_=p.get("a_5");
    p_.a_max_psi_Neuron_=p.get("a_max_psi_Neuron");
    p_.a_max_psi_1_=p.get("a_max_psi_1");
    p_.a_max_psi_2_=p.get("a_max_psi_2");
    p_.a_max_psi_3_=p.get("a_max_psi_3");
    p_.a_max_psi_4_=p.get("a_max_psi_4");
    p_.a_max_psi_5_=p.get("a_max_psi_5");

    p_.a_max_omega_Neuron_=p.get("a_max_omega_Neuron");
    p_.a_max_omega_1_=p.get("a_max_omega_1");
    p_.a_max_omega_2_=p.get("a_max_omega_2");
    p_.a_max_omega_3_=p.get("a_max_omega_3");
    p_.a_max_omega_4_=p.get("a_max_omega_4");
    p_.a_max_omega_5_=p.get("a_max_omega_5");


    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model125(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model131:public BaseModel
{
  SimplestModel m;

  Model131* clone()const override { return new Model131;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;


    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_1_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_3_;
    s.N_[5]=p.N_4_;
    s.N_[6]=p.N_5_;


    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_1_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_3_;
    s.M_[5]=p.N_4_;
    s.M_[6]=p.N_5_;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model131(){}
  ~Model131(){}
  virtual std::string id() const  override
  {
    return "Model 1.31";
  }
  static double number()
  {
    return 1.31;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }
  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model131(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model132:public BaseModel
{
  SimplestModel m;

  Model132* clone()const override { return new Model132;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_1_;
    double N_2_;
    double N_3_;
    double N_4_;
    double N_5_;
    double N_N_;
    double M_0_;
    double M_1_;
    double M_2_;
    double M_3_;
    double M_4_;
    double M_5_;
    double M_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;
    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;


    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_1_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_3_;
    s.N_[5]=p.N_4_;
    s.N_[6]=p.N_5_;


    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.M_N_;
    s.M_[1]=p.M_0_;
    s.M_[2]=p.M_1_;
    s.M_[3]=p.M_2_;
    s.M_[4]=p.M_3_;
    s.M_[5]=p.M_4_;
    s.M_[6]=p.M_5_;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model132(){}
  ~Model132(){}
  virtual std::string id() const  override
  {
    return "Model 1.32";
  }
  static double number()
  {
    return 1.32;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_1",p_.N_1_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_3",p_.N_3_ );
    out.push_back("N_4",p_.N_4_ );
    out.push_back("N_5",p_.N_5_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("M_0",p_.M_0_ );
    out.push_back("M_1",p_.M_1_ );
    out.push_back("M_2",p_.M_2_ );
    out.push_back("M_3",p_.M_3_ );
    out.push_back("M_4",p_.M_4_ );
    out.push_back("M_5",p_.M_5_ );
    out.push_back("M_N",p_.M_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);
    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_1_=p.get("N_1");
    p_.N_2_=p.get("N_2");
    p_.N_3_=p.get("N_3");
    p_.N_4_=p.get("N_4");
    p_.N_5_=p.get("N_5");
    p_.N_N_=p.get("N_N");
    p_.M_0_=p.get("M_0");
    p_.M_1_=p.get("M_1");
    p_.M_2_=p.get("M_2");
    p_.M_3_=p.get("M_3");
    p_.M_4_=p.get("M_4");
    p_.M_5_=p.get("M_5");
    p_.M_N_=p.get("M_N");


    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model132(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model141:public BaseModel
{
  SimplestModel m;

  Model141* clone()const override { return new Model141;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_1_;
    double k_sig_2_;
    double k_sig_3_;
    double k_sig_4_;
    double k_sig_5_;
    double k_sig_max_;
    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[2]=p.k_sig_1_;
    s.ksig_omega_[3]=p.k_sig_2_;
    s.ksig_omega_[4]=p.k_sig_3_;
    s.ksig_omega_[5]=p.k_sig_4_;
    s.ksig_omega_[6]=p.k_sig_5_;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model141(){}
  ~Model141(){}
  virtual std::string id() const  override
  {
    return "Model 1.41";
  }
  static double number()
  {
    return 1.41;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig_1",p_.k_sig_1_);
    out.push_back("k_sig_2",p_.k_sig_2_);
    out.push_back("k_sig_3",p_.k_sig_3_);
    out.push_back("k_sig_4",p_.k_sig_4_);
    out.push_back("k_sig_5",p_.k_sig_5_);

    out.push_back("k_sig_max",p_.k_sig_max_);


    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.k_sig_1_=p.get("k_sig_1");
    p_.k_sig_2_=p.get("k_sig_2");
    p_.k_sig_3_=p.get("k_sig_3");
    p_.k_sig_4_=p.get("k_sig_4");
    p_.k_sig_5_=p.get("k_sig_5");

    p_.k_sig_max_=p.get("k_sig_max");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model141(const Parameters& p)
  {
    loadParameters(p);
  }


};










class Model142:public BaseModel
{
  SimplestModel m;

  Model142* clone()const override { return new Model142;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_1;
    double k_sig_max_2;
    double k_sig_max_3;
    double k_sig_max_4;
    double k_sig_max_5;


    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_[2]=p.k_sig_max_1;
    s.ksig_max_omega_[3]=p.k_sig_max_2;
    s.ksig_max_omega_[4]=p.k_sig_max_3;
    s.ksig_max_omega_[5]=p.k_sig_max_4;
    s.ksig_max_omega_[6]=p.k_sig_max_5;

    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[2]=p.k_sig_max_1;
    s.ksig_max_psi_[3]=p.k_sig_max_2;
    s.ksig_max_psi_[4]=p.k_sig_max_3;
    s.ksig_max_psi_[5]=p.k_sig_max_4;
    s.ksig_max_psi_[6]=p.k_sig_max_5;



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;

    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model142(){}
  ~Model142(){}
  virtual std::string id() const  override
  {
    return "Model 1.42";
  }
  static double number()
  {
    return 1.42;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max_1",p_.k_sig_max_1);
    out.push_back("k_sig_max_2",p_.k_sig_max_2);
    out.push_back("k_sig_max_3",p_.k_sig_max_3);
    out.push_back("k_sig_max_4",p_.k_sig_max_4);
    out.push_back("k_sig_max_5",p_.k_sig_max_5);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");

    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_1=p.get("k_sig_max_1");
    p_.k_sig_max_2=p.get("k_sig_max_2");
    p_.k_sig_max_3=p.get("k_sig_max_3");
    p_.k_sig_max_4=p.get("k_sig_max_4");
    p_.k_sig_max_5=p.get("k_sig_max_5");



    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }
  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model142(const Parameters& p)
  {
    loadParameters(p);
  }


};







class Model144:public BaseModel
{
  SimplestModel m;

  Model144* clone()const override { return new Model144;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_omega_1;
    double k_sig_max_omega_2;
    double k_sig_max_omega_3;
    double k_sig_max_omega_4;
    double k_sig_max_omega_5;

    double k_sig_max_psi_1;
    double k_sig_max_psi_2;
    double k_sig_max_psi_3;
    double k_sig_max_psi_4;
    double k_sig_max_psi_5;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_omega_[2]=p.k_sig_max_omega_1;
    s.ksig_max_omega_[3]=p.k_sig_max_omega_2;
    s.ksig_max_omega_[4]=p.k_sig_max_omega_3;
    s.ksig_max_omega_[5]=p.k_sig_max_omega_4;
    s.ksig_max_omega_[6]=p.k_sig_max_omega_5;

    s.ksig_max_psi_=std::vector<double>(7,0);
    s.ksig_max_psi_[2]=p.k_sig_max_psi_1;
    s.ksig_max_psi_[3]=p.k_sig_max_psi_2;
    s.ksig_max_psi_[4]=p.k_sig_max_psi_3;
    s.ksig_max_psi_[5]=p.k_sig_max_psi_4;
    s.ksig_max_psi_[6]=p.k_sig_max_psi_5;



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model144(){}
  ~Model144(){}
  virtual std::string id() const  override
  {
    return "Model 1.44";
  }
  static double number()
  {
    return 1.44;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max_omega_1",p_.k_sig_max_omega_1);
    out.push_back("k_sig_max_omega_2",p_.k_sig_max_omega_2);
    out.push_back("k_sig_max_omega_3",p_.k_sig_max_omega_3);
    out.push_back("k_sig_max_omega_4",p_.k_sig_max_omega_4);
    out.push_back("k_sig_max_omega_5",p_.k_sig_max_omega_5);

    out.push_back("k_sig_max_psi_1",p_.k_sig_max_psi_1);
    out.push_back("k_sig_max_psi_2",p_.k_sig_max_psi_2);
    out.push_back("k_sig_max_psi_3",p_.k_sig_max_psi_3);
    out.push_back("k_sig_max_psi_4",p_.k_sig_max_psi_4);
    out.push_back("k_sig_max_psi_5",p_.k_sig_max_psi_5);


    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_omega_1=p.get("k_sig_max_omega_1");
    p_.k_sig_max_omega_2=p.get("k_sig_max_omega_2");
    p_.k_sig_max_omega_3=p.get("k_sig_max_omega_3");
    p_.k_sig_max_omega_4=p.get("k_sig_max_omega_4");
    p_.k_sig_max_omega_5=p.get("k_sig_max_omega_5");


    p_.k_sig_max_psi_1=p.get("k_sig_max_psi_1");
    p_.k_sig_max_psi_2=p.get("k_sig_max_psi_2");
    p_.k_sig_max_psi_3=p.get("k_sig_max_psi_3");
    p_.k_sig_max_psi_4=p.get("k_sig_max_psi_4");
    p_.k_sig_max_psi_5=p.get("k_sig_max_psi_5");

    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }


  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model144(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model151:public BaseModel
{
  SimplestModel m;

  Model151* clone()const override { return new Model151;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double Kg_2_;
    double Kg_3_;
    double Kg_4_;
    double Kg_5_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_psi_[2]=p.Kg_2_;
    s.Keq_gmax_psi_[3]=p.Kg_3_;
    s.Keq_gmax_psi_[4]=p.Kg_4_;
    s.Keq_gmax_psi_[5]=p.Kg_5_;

    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);
    s.Keq_gmax_omega_[2]=p.Kg_2_;
    s.Keq_gmax_omega_[3]=p.Kg_3_;
    s.Keq_gmax_omega_[4]=p.Kg_4_;
    s.Keq_gmax_omega_[5]=p.Kg_5_;


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model151(){}
  ~Model151(){}
  virtual std::string id() const  override
  {
    return "Model 1.51";
  }
  static double number()
  {
    return 1.51;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("Kg_2",p_.Kg_2_);
    out.push_back("Kg_3",p_.Kg_3_);
    out.push_back("Kg_4",p_.Kg_4_);
    out.push_back("Kg_5",p_.Kg_5_);


    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.Kg_2_=p.get("Kg_2");
    p_.Kg_3_=p.get("Kg_3");
    p_.Kg_4_=p.get("Kg_4");
    p_.Kg_5_=p.get("Kg_5");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model151(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model152:public BaseModel
{
  SimplestModel m;

  Model152* clone()const override { return new Model152;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double Kg_psi_2_;
    double Kg_psi_3_;
    double Kg_psi_4_;
    double Kg_psi_5_;

    double Kg_omega_2_;
    double Kg_omega_3_;
    double Kg_omega_4_;
    double Kg_omega_5_;

    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;

    double k_sig_;
    double k_sig_max_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_psi_[2]=p.Kg_psi_2_;
    s.Keq_gmax_psi_[3]=p.Kg_psi_3_;
    s.Keq_gmax_psi_[4]=p.Kg_psi_4_;
    s.Keq_gmax_psi_[5]=p.Kg_psi_5_;

    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);
    s.Keq_gmax_omega_[2]=p.Kg_omega_2_;
    s.Keq_gmax_omega_[3]=p.Kg_omega_3_;
    s.Keq_gmax_omega_[4]=p.Kg_omega_4_;
    s.Keq_gmax_omega_[5]=p.Kg_omega_5_;


    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);

    s.ksig_max_psi_[3]=p.k_sig_max_;
    s.ksig_max_psi_[4]=p.k_sig_max_*1.5;
    s.ksig_max_psi_[5]=p.k_sig_max_*3;
    s.ksig_max_psi_[6]=p.k_sig_max_*6;

    s.ksig_max_omega_[3]=p.k_sig_max_;
    s.ksig_max_omega_[4]=p.k_sig_max_*1.5;
    s.ksig_max_omega_[5]=p.k_sig_max_*3;
    s.ksig_max_omega_[6]=p.k_sig_max_*6;




    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model152(){}
  ~Model152(){}
  virtual std::string id() const  override
  {
    return "Model 1.52";
  }
  static double number()
  {
    return 1.52;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("Kg_psi_2",p_.Kg_psi_2_);
    out.push_back("Kg_psi_3",p_.Kg_psi_3_);
    out.push_back("Kg_psi_4",p_.Kg_psi_4_);
    out.push_back("Kg_psi_5",p_.Kg_psi_5_);

    out.push_back("Kg_omega_2",p_.Kg_omega_2_);
    out.push_back("Kg_omega_3",p_.Kg_omega_3_);
    out.push_back("Kg_omega_4",p_.Kg_omega_4_);
    out.push_back("Kg_omega_5",p_.Kg_omega_5_);

    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);
    out.push_back("k_sig_max",p_.k_sig_max_);

    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);

    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.Kg_psi_2_=p.get("Kg_psi_2");
    p_.Kg_psi_3_=p.get("Kg_psi_3");
    p_.Kg_psi_4_=p.get("Kg_psi_4");
    p_.Kg_psi_5_=p.get("Kg_psi_5");

    p_.Kg_omega_2_=p.get("Kg_omega_2");
    p_.Kg_omega_3_=p.get("Kg_omega_3");
    p_.Kg_omega_4_=p.get("Kg_omega_4");
    p_.Kg_omega_5_=p.get("Kg_omega_5");

    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_=p.get("DAMP_ratio");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");
    p_.k_sig_=p.get("k_sig");
    p_.k_sig_max_=p.get("k_sig_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model152(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model20:public BaseModel
{
  SimplestModel m;

  Model20* clone()const override { return new Model20;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_psi_;
    double DAMP_ratio_omega_;

    double prot_concentration_;
    double DAMP_MW_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_psi_;
    s.DAMP_omega_ratio_=p.DAMP_ratio_omega_;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model20(){}
  ~Model20(){}
  virtual std::string id() const  override
  {
    return "Model 2.0";
  }
  static double number()
  {
    return 2.0;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio_psi",p_.DAMP_ratio_psi_);

    out.push_back("DAMP_ratio_omega",p_.DAMP_ratio_omega_);

    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_psi_=p.get("DAMP_ratio_psi");
    p_.DAMP_ratio_omega_=p.get("DAMP_ratio_omega");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }



  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model20(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model211:public BaseModel
{
  SimplestModel m;

  Model211* clone()const override { return new Model211;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_12_;
    double g_23_;
    double g_34_;
    double g_45_;

    double g_10_;
    double g_21_;

    double g_max_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_psi_;
    double DAMP_ratio_omega_;
    double prot_concentration_;
    double DAMP_MW_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_psi_;
    s.DAMP_omega_ratio_=p.DAMP_ratio_omega_;
    s.DAMP_omega_ratio_=0;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);





    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;


    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_;
    s.g_max_psi_[3]=p.g_max_;
    s.g_max_psi_[4]=p.g_max_;
    s.g_max_psi_[5]=p.g_max_;

    s.g_max_omega_[2]=p.g_max_;
    s.g_max_omega_[3]=p.g_max_;
    s.g_max_omega_[4]=p.g_max_;
    s.g_max_omega_[5]=p.g_max_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;
    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;


    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model211(){}
  ~Model211(){}
  virtual std::string id() const  override
  {
    return "Model 2.11";
  }
  static double number()
  {
    return 2.11;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max",p_.g_max_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio_psi",p_.DAMP_ratio_psi_);

    out.push_back("DAMP_ratio_omega",p_.DAMP_ratio_omega_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_=p.get("g_max");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_psi_=p.get("DAMP_ratio_psi");
    p_.DAMP_ratio_omega_=p.get("DAMP_ratio_omega");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");



    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }


  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  Model211(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model212:public BaseModel
{
  SimplestModel m;

  Model212* clone()const override { return new Model212;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_psi_;
    double DAMP_ratio_omega_;

    double prot_concentration_;
    double DAMP_MW_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_psi_;
    s.DAMP_omega_ratio_=p.DAMP_ratio_omega_;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model212(){}
  ~Model212(){}
  virtual std::string id() const  override
  {
    return "Model 2.12";
  }
  static double number()
  {
    return 2.12;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio_psi",p_.DAMP_ratio_psi_);

    out.push_back("DAMP_ratio_omega",p_.DAMP_ratio_omega_);

    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_psi_=p.get("DAMP_ratio_psi");
    p_.DAMP_ratio_omega_=p.get("DAMP_ratio_omega");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }


  Model212(const Parameters& p)
  {
    loadParameters(p);
  }



};


class Model213:public BaseModel
{
  SimplestModel m;

  Model213* clone()const override { return new Model213;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_34_;
    double g_45_;
    double g_max_2_;
    double g_max_3_;
    double g_max_4_;
    double g_max_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_psi_;
    double DAMP_ratio_omega_;

    double prot_concentration_;
    double DAMP_MW_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_psi_;
    s.DAMP_omega_ratio_=p.DAMP_ratio_omega_;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_2_;
    s.g_max_psi_[3]=p.g_max_3_;
    s.g_max_psi_[4]=p.g_max_4_;
    s.g_max_psi_[5]=p.g_max_5_;

    s.g_max_omega_[2]=p.g_max_2_;
    s.g_max_omega_[3]=p.g_max_3_;
    s.g_max_omega_[4]=p.g_max_4_;
    s.g_max_omega_[5]=p.g_max_5_;


    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model213(){}
  ~Model213(){}
  virtual std::string id() const  override
  {
    return "Model 2.13";
  }
  static double number()
  {
    return 2.13;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );
    out.push_back("g_max_2",p_.g_max_2_ );
    out.push_back("g_max_3",p_.g_max_3_ );
    out.push_back("g_max_4",p_.g_max_4_ );
    out.push_back("g_max_5",p_.g_max_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio_psi",p_.DAMP_ratio_psi_);

    out.push_back("DAMP_ratio_omega",p_.DAMP_ratio_omega_);

    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_2_=p.get("g_max_2");
    p_.g_max_3_=p.get("g_max_3");
    p_.g_max_4_=p.get("g_max_4");
    p_.g_max_5_=p.get("g_max_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_psi_=p.get("DAMP_ratio_psi");
    p_.DAMP_ratio_omega_=p.get("DAMP_ratio_omega");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  Model213(const Parameters& p)
  {
    loadParameters(p);
  }


};


class Model214:public BaseModel
{
  SimplestModel m;

  Model214* clone()const override { return new Model214;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_psi_;
    double DAMP_ratio_omega_;

    double prot_concentration_;
    double DAMP_MW_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_psi_;
    s.DAMP_omega_ratio_=p.DAMP_ratio_omega_;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_23_;
    s.g_rigth_[5]=p.g_23_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;



    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model214(){}
  ~Model214(){}
  virtual std::string id() const  override
  {
    return "Model 2.14";
  }
  static double number()
  {
    return 2.14;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio_psi",p_.DAMP_ratio_psi_);

    out.push_back("DAMP_ratio_omega",p_.DAMP_ratio_omega_);

    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_psi_=p.get("DAMP_ratio_psi");
    p_.DAMP_ratio_omega_=p.get("DAMP_ratio_omega");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  Model214(const Parameters& p)
  {
    loadParameters(p);
  }


};




class Model215:public BaseModel
{
  SimplestModel m;

  Model215* clone()const override { return new Model215;}


  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_psi_;
    double Keq_omega_;
    double kcat_psi;
    double kcat_omega_;
    double g_01_;
    double g_10_;
    double g_12_;
    double g_21_;

    double g_23_;
    double g_34_;
    double g_45_;
    double g_max_psi_2_;
    double g_max_psi_3_;
    double g_max_psi_4_;
    double g_max_psi_5_;
    double g_max_omega_2_;
    double g_max_omega_3_;
    double g_max_omega_4_;
    double g_max_omega_5_;
    double N_0_;
    double N_2_;
    double N_N_;
    double N_Astr_;
    double N_Neuron_;
    double a_2_;
    double a_factor_;
    double a_max_Neuron_;
    double a_max_;

    double inj_width_;
    double DAMP_ratio_psi_;
    double DAMP_ratio_omega_;

    double prot_concentration_;
    double DAMP_MW_;

    double inj_width_3_;
    double inj_width_7_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.inj_width_=p.inj_width_;
    s.DAMP_psi_ratio_=p.DAMP_ratio_psi_;
    s.DAMP_omega_ratio_=p.DAMP_ratio_omega_;

    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.Keq_gmax_psi_=std::vector<double>(7,p.Keq_psi_);
    s.Keq_gmax_omega_=std::vector<double>(7,p.Keq_omega_);


    s.ksig_omega_=std::vector<double>(7,0);

    s.ksig_max_omega_=std::vector<double>(7,0);
    s.ksig_max_psi_=std::vector<double>(7,0);



    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;
    s.g_left_[3]=p.g_21_;



    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;
    s.g_rigth_[2]=p.g_12_;
    s.g_rigth_[3]=p.g_23_;
    s.g_rigth_[4]=p.g_34_;
    s.g_rigth_[5]=p.g_45_;

    s.g_max_omega_=std::vector<double> (7,0);


    s.g_max_psi_=std::vector<double> (7,0);

    s.g_max_psi_[2]=p.g_max_psi_2_;
    s.g_max_psi_[3]=p.g_max_psi_3_;
    s.g_max_psi_[4]=p.g_max_psi_4_;
    s.g_max_psi_[5]=p.g_max_psi_5_;

    s.g_max_omega_[2]=p.g_max_omega_2_;
    s.g_max_omega_[3]=p.g_max_omega_3_;
    s.g_max_omega_[4]=p.g_max_omega_4_;
    s.g_max_omega_[5]=p.g_max_omega_5_;



    s.a_=std::vector<double> (7,0);
    s.a_[3]=p.a_2_;
    s.a_[4]=p.a_2_*p.a_factor_;
    s.a_[5]=s.a_[4]*p.a_factor_;
    s.a_[6]=s.a_[5]*p.a_factor_;

    s.a_omega_=std::vector<double> (7,0);
    s.a_psi_=std::vector<double> (7,0);

    s.a_psi_[0]=p.a_max_Neuron_;

    s.a_psi_[3]=p.a_max_;
    s.a_psi_[4]=p.a_max_*p.a_factor_;
    s.a_psi_[5]=s.a_psi_[4]*p.a_factor_;
    s.a_psi_[6]=s.a_psi_[5]*p.a_factor_;


    s.a_omega_[0]=p.a_max_Neuron_;

    s.a_omega_[0]=p.a_max_Neuron_;
    s.a_omega_[3]=p.a_max_;
    s.a_omega_[4]=p.a_max_*p.a_factor_;
    s.a_omega_[5]=s.a_omega_[4]*p.a_factor_;
    s.a_omega_[6]=s.a_omega_[5]*p.a_factor_;



    s.N_=std::vector<double> (7,0);

    s.N_[0]=p.N_N_;
    s.N_[1]=p.N_0_;
    s.N_[2]=p.N_0_;
    s.N_[3]=p.N_2_;
    s.N_[4]=p.N_2_*1.5;
    s.N_[5]=p.N_2_*3;
    s.N_[6]=p.N_2_*6;



    s.M_=std::vector<double> (7,0);
    s.M_[0]=p.N_N_;
    s.M_[1]=p.N_0_;
    s.M_[2]=p.N_0_;
    s.M_[3]=p.N_2_;
    s.M_[4]=p.N_2_*1.5;
    s.M_[5]=p.N_2_*3;
    s.M_[6]=p.N_2_*6;





    s.dens_Astr_=p.N_Astr_;



    s.dens_Neur_=p.N_Neuron_;


    return s;
  }


  myParameters p_;

  // BaseModel interface
public:
  Model215(){}
  ~Model215(){}
  virtual std::string id() const  override
  {
    return "Model 2.15";
  }
  static double number()
  {
    return 2.15;
  }
  virtual Parameters getParameters() const override
  {
    Parameters out;
    out.push_back("model",number());
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
    out.push_back("g_12",p_.g_12_ );
    out.push_back("g_21",p_.g_21_ );

    out.push_back("g_23",p_.g_23_ );
    out.push_back("g_34",p_.g_34_ );
    out.push_back("g_45",p_.g_45_ );

    out.push_back("g_max_psi_2",p_.g_max_psi_2_ );
    out.push_back("g_max_psi_3",p_.g_max_psi_3_ );
    out.push_back("g_max_psi_4",p_.g_max_psi_4_ );
    out.push_back("g_max_psi_5",p_.g_max_psi_5_ );
    out.push_back("g_max_omega_2",p_.g_max_omega_2_ );
    out.push_back("g_max_omega_3",p_.g_max_omega_3_ );
    out.push_back("g_max_omega_4",p_.g_max_omega_4_ );
    out.push_back("g_max_omega_5",p_.g_max_omega_5_ );
    out.push_back("N_0",p_.N_0_ );
    out.push_back("N_2",p_.N_2_ );
    out.push_back("N_N",p_.N_N_ );
    out.push_back("N_Astr",p_.N_Astr_);
    out.push_back("N_Neuron_",p_.N_Neuron_);
    out.push_back("a_2",p_.a_2_ );
    out.push_back("a_factor",p_.a_factor_ );

    out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
    out.push_back("a_max",p_.a_max_ );



    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio_psi",p_.DAMP_ratio_psi_);

    out.push_back("DAMP_ratio_omega",p_.DAMP_ratio_omega_);

    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);



    out.push_back("inj_width_3dpl",p_.inj_width_3_);
    out.push_back("inj_width_7dpl2",p_.inj_width_7_);







    return out;

  }
  virtual void loadParameters(const Parameters& p) override
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
    p_.g_12_=p.get("g_12");
    p_.g_21_=p.get("g_21");

    p_.g_23_=p.get("g_23");
    p_.g_34_=p.get("g_34");
    p_.g_45_=p.get("g_45");
    p_.g_max_psi_2_=p.get("g_max_psi_2");
    p_.g_max_psi_3_=p.get("g_max_psi_3");
    p_.g_max_psi_4_=p.get("g_max_psi_4");
    p_.g_max_psi_5_=p.get("g_max_psi_5");
    p_.g_max_omega_2_=p.get("g_max_omega_2");
    p_.g_max_omega_3_=p.get("g_max_omega_3");
    p_.g_max_omega_4_=p.get("g_max_omega_4");
    p_.g_max_omega_5_=p.get("g_max_omega_5");
    p_.N_0_=p.get("N_0");
    p_.N_2_=p.get("N_2");
    p_.N_N_=p.get("N_N");
    p_.a_2_=p.get("a_2");
    p_.DAMP_ratio_psi_=p.get("DAMP_ratio_psi");
    p_.DAMP_ratio_omega_=p.get("DAMP_ratio_omega");
    p_.DAMP_MW_=p.get("DAMP_MW");
    p_.prot_concentration_=p.get("prot_concentration");
    p_.inj_width_=p.get("inj_width");
    p_.N_Astr_=p.get("N_Astr");
    p_.N_Neuron_=p.get("N_Neuron");

    p_.a_factor_=p.get("a_factor");
    p_.a_max_Neuron_=p.get("a_max_Neuron");
    p_.a_max_=p.get("a_max");


    p_.inj_width_3_=p.get("inj_width_3dpl");
    p_.inj_width_7_=p.get("inj_width_7dpl2");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax,double teq) const override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq);

  }

  virtual CortexSimulation run_dt(const Experiment& e,double dx,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,double tequilibrio) const  override
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dx,dts,tequilibrio);
  }

  virtual std::pair<CortexSimulation, std::vector<double>> run(const Experiment& e,double dx,double dtmin0, double dtmin, std::size_t nPoints_per_decade, double dtmax, double teq, double maxlogError,double dtinf) const override
  {
    return m.simulate_Adapted(getParameters(),toModelParameters(this->p_),e,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,teq,maxlogError,dtinf);

  }

  virtual const SimplestModel& myModel()const override{ return m;}
  virtual SimplestModel& myModel()override {return m;}
  virtual SimplestModel::Param toModelParameters()const override{return toModelParameters(this->p_);}


  Model215(const Parameters& p)
  {
    loadParameters(p);
  }


};



#endif // MODELS

