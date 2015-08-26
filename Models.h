#ifndef MODELS
#define MODELS

#include "CortexMeasure.h"
#include "CortexSimulation.h"

#include "BayesIteration.h"

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
    std::vector<double> damp_;
    double Dpsi_;
    double Domega_;
    double epsilon_;

    double kon_psi_;
    double kcat_psi_;
    double kon_omega_;
    double kcat_omega_;

    std::vector<double> ksig_omega_;
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
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;



  };

  SimplestModel() {}
  CortexState nextEuler(const Param &p,const CortexState& c, double dt) const;

  CortexState init(const Param &p, const CortexExperiment &s) const;


  CortexSimulation simulate(const Parameters& par,const Param &p, const CortexExperiment &sp, double dt)const;

  CortexSimulation simulate(const Parameters& par, const Param &p, const Experiment &sp, double dt, double tequilibrio)const;





  void addDamp(CortexState& c,const Param& p)const
  {
    unsigned i=0;
    double damp=p.prot_concentration_*p.DAMP_ratio_/p.DAMP_MW_/1000;
    while (c.x_[i]<p.inj_width_+c.x_[0])
      {
        double inj_size_cell=std::min(p.inj_width_+c.x_[0],c.x_[i]+c.dx_[i])-c.x_[i];

        c.psi_T_[i]+=damp*inj_size_cell/c.dx_[i]/p.epsilon_;
        ++i;
      }
  }





  std::vector<double> dPsi_dt(const Param &p, const CortexState& c) const;


  std::vector<double> Psi_Bound(const Param &p, const CortexState& c) const;


  std::vector<double> dOmega_dt(const Param &p, const CortexState &c)const;

  std::vector<double> Omega_Bound(const Param &p, const CortexState& c) const;




  std::vector<std::vector<double> > dRho_dt(const Param &p, const CortexState &c, bool hasOmega)const;
  CortexState init(const SimplestModel::Param &p, const Experiment &s) const;
};





class BaseModel
{
public:
  virtual std::string id()const=0;

  virtual Parameters getParameters()const=0;

  virtual void loadParameters(const Parameters& p)=0;

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const=0;
  virtual CortexSimulation run(const Experiment& e,double dt,  double tequilibrio) const=0;


  static BaseModel* create(const Parameters& p);

  virtual ~BaseModel(){}

  virtual std::vector<double> getObservedProbFromModel(const std::vector<double>& modelRho)const;

  virtual std::vector<double> getObservedNumberFromData(const std::vector<double>& modelRho)const;

  double likelihood(const Experiment& m, const CortexSimulation& s)const
  {
    unsigned is=0;
    double sumLogLik=0;

    for (unsigned ie=0; ie<m.numMeasures(); ie++)
      {
        const CortexMeasure* cm=m.getMeasure(ie);
        double t=cm->dia()*24*60*60;
        while (s.t_[is]<t
               &&is<s.t_.size())
          ++is;
        std::vector<double> rho_sim;
        std::vector<double> rho_meas;
        for (unsigned ix=0; ix<cm->dx().size(); ++ix)
          {
            rho_sim=getObservedProbFromModel(s.rho_[is][ix]);
            auto ob=cm->meanAstro(ix);
            rho_meas=getObservedNumberFromData(ob);
            sumLogLik+=MultinomialLikelihood(rho_meas,rho_sim);
          }

      }
    return sumLogLik;

  }






private:
  static std::map<double,BaseModel*> models_;
  static std::map<double,BaseModel*> getModels();
};




















class Model00:public BaseModel
{
  SimplestModel m;




  class myParameters
  {
  public:
    double D_;
    double epsilon_;
    double Keq_;
    double kcat_;
    double g_01_;
    double g_10_;
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
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double DAMP_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.damp_=std::vector<double>(1);
    s.inj_width_=p.inj_width_;
    s.DAMP_ratio_=p.DAMP_ratio_;
    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.damp_[0]=p.DAMP_;
    s.Dpsi_=p.D_;
    s.Domega_=0;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_/p.Keq_;
    s.kcat_psi_=p.kcat_;
    s.kon_omega_=0;
    s.kcat_omega_=0;

    s.ksig_omega_=std::vector<double>(7,0);

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;


    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;

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
  virtual std::string id() const
  {
    return "Model 0.0";
  }
  static double number()
  {
    return 0;
  }
  virtual Parameters getParameters() const
  {
    Parameters out;
    out.push_back("model",0.0);
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq",p_.Keq_);
    out.push_back("kcat", p_.kcat_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
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
    out.push_back("DAMP",p_.DAMP_);


    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);





    return out;

  }
  virtual void loadParameters(const Parameters& p)
  {
    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_=p.get("Keq");
    p_.kcat_=p.get("kcat");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
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



  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);

  }

  virtual CortexSimulation run(const Experiment& e,double dt, double teq) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt,teq);

  }





  Model00(const Parameters& p)
  {
    loadParameters(p);
  }


};



class Model10:public BaseModel
{
  SimplestModel m;

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
    double inj_width_;
    double DAMP_ratio_;
    double prot_concentration_;
    double DAMP_MW_;
    double DAMP_;
    double k_sig_;
  };


  SimplestModel::Param toModelParameters(const myParameters& p)const
  {
    SimplestModel::Param s;
    s.damp_=std::vector<double>(1);
    s.inj_width_=p.inj_width_;
    s.DAMP_ratio_=p.DAMP_ratio_;
    s.prot_concentration_=p.prot_concentration_;
    s.DAMP_MW_=p.DAMP_MW_;
    s.damp_[0]=p.DAMP_;
    s.Dpsi_=p.D_;
    s.Domega_=p.D_;
    s.epsilon_=p.epsilon_;

    s.kon_psi_=p.kcat_psi/p.Keq_psi_;
    s.kcat_psi_=p.kcat_psi;
    s.kon_omega_=p.kcat_omega_/p.Keq_omega_;
    s.kcat_omega_=p.kcat_omega_;

    s.ksig_omega_=std::vector<double>(7,0);
    s.ksig_omega_[3]=p.k_sig_;
    s.ksig_omega_[4]=p.k_sig_*1.5;
    s.ksig_omega_[5]=p.k_sig_*3;
    s.ksig_omega_[6]=p.k_sig_*6;

    s.g_left_=std::vector<double> (7,0);
    s.g_left_[2]=p.g_10_;


    s.g_rigth_=std::vector<double> (7,0);
    s.g_rigth_[1]=p.g_01_;

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
  virtual std::string id() const
  {
    return "Model 1.0";
  }
  static double number()
  {
    return 1;
  }
  virtual Parameters getParameters() const
  {
    Parameters out;
    out.push_back("model",1.0);
    out.push_back("D",p_.D_);
    out.push_back("epsilon",p_.epsilon_);
    out.push_back("Keq_psi",p_.Keq_psi_);
    out.push_back("Keq_omega",p_.Keq_omega_);
    out.push_back("kcat_psi", p_.kcat_psi);
    out.push_back("kcat_omega", p_.kcat_omega_);
    out.push_back("g_01",p_.g_01_);
    out.push_back("g_10",p_.g_10_ );
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
    out.push_back("DAMP",p_.DAMP_);
    out.push_back("inj_width",p_.inj_width_);
    out.push_back("DAMP_ratio",p_.DAMP_ratio_);
    out.push_back("prot_concentration",p_.prot_concentration_);
    out.push_back("DAMP_MW",p_.DAMP_MW_);


    out.push_back("k_sig",p_.k_sig_);






    return out;

  }
  virtual void loadParameters(const Parameters& p)
  {

    p_.D_=p.get("D");
    p_.epsilon_=p.get("epsilon");
    p_.Keq_psi_=p.get("Keq_psi");
    p_.Keq_omega_=p.get("Keq_omega");
    p_.kcat_psi=p.get("kcat_psi");
    p_.kcat_omega_=p.get("kcat_omega");
    p_.g_01_=p.get("g_01");
    p_.g_10_=p.get("g_10");
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
    p_.k_sig_=p.get("k_sig");
  }

  virtual CortexSimulation run(const CortexExperiment& e,double dt) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt);


  }
  virtual CortexSimulation run(const Experiment& e,double dt,double teq) const
  {
    return m.simulate(getParameters(),toModelParameters(this->p_),e,dt,teq);

  }


  Model10(const Parameters& p)
  {
    loadParameters(p);
  }


};





#endif // MODELS

