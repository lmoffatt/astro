#include "Models.h"

#include <iostream>


CortexState SimplestModel::nextEuler(const SimplestModel::Param &p, const CortexState &c, double dt)const
{
  auto dPsi=dPsi_dt(p,c);
  bool hasOmega=p.Domega_>0;
  std::vector<double> dOmega;
  if (hasOmega)
    dOmega=dOmega_dt(p,c);
  auto dRho=dRho_dt(p,c,hasOmega);


  CortexState out(c);
  double maxf=0;

  for (double x:dPsi)
    if (std::abs(x)>maxf)
      maxf=std::abs(x);

  if (hasOmega)
    {
      out.omega_T_+=dOmega*dt;
    }
  out.psi_T_+=dPsi*dt;



  out.rho_+=dRho*dt;
  if (hasOmega)
    out.omega_B_=Omega_Bound(p,out);
  out.psi_B_=Psi_Bound(p,out);

  return out;


}

CortexState SimplestModel::init(const SimplestModel::Param &p, const CortexExperiment &s)const
{
  CortexState o(s.x_,s.dx_,s.h_,p.epsilon_,p.N_.size());


  unsigned numX=o.x_.size();
  for (unsigned x=0; x<numX; ++x)
    {
      o.rho_[x][0]=p.dens_Neur_*std::pow(s.h_,2)*s.dx_[x]*1000;  // dens is in cells per liter
      o.rho_[x][1]=p.dens_Astr_*std::pow(s.h_,2)*s.dx_[x]*1000;

    }
  o.psi_B_=Psi_Bound(p,o);
  o.omega_B_=Omega_Bound(p,o);
  return o;
}



std::vector<double> SimplestModel::dPsi_dt(const Param &p,
                                           const CortexState &c)const
{
  unsigned numX=c.rho_.size();

  std::vector<double> o(numX,0);


  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      if (x>0)
        Jn=2.0*p.Dpsi_*(c.psi_T_[x-1]-c.psi_B_[x-1]-c.psi_T_[x]+c.psi_B_[x])/(c.dx_[x-1]+c.dx_[x]);
      if (x+1<numX)
        Jp=2.0*p.Dpsi_*(c.psi_T_[x+1]-c.psi_B_[x+1]-c.psi_T_[x]+c.psi_B_[x])/(c.dx_[x+1]+c.dx_[x]);

      o[x]=(Jn+Jp)/c.dx_[x]-p.kcat_psi_*c.psi_B_[x];
    }
  return o;
}

std::vector<double> SimplestModel::Psi_Bound(const SimplestModel::Param &p, const CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;

  std::vector<double> o(numX,0);

  double molar_section=1.0/(NAv*p.epsilon_*c.h_*c.h_*1000.0);

  double K_psi=p.kcat_psi_/p.kon_psi_;
  for (unsigned x=0; x<numX; ++x)
    {
      double Nt=0;
      for (unsigned k=0; k<numK; ++k)
        Nt+=p.N_[k]*c.rho_[x][k];
      double R_T=Nt*molar_section/c.dx_[x];
      double b=R_T+c.psi_T_[x]+K_psi;
      o[x]=0.5*(b)-0.5*std::sqrt(b*b-4*R_T*c.psi_T_[x]);
    }
  return o;
}

std::vector<double> SimplestModel::dOmega_dt(const Param &p,
                                             const CortexState &c)const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  const double NAv=6.022E23;
  std::vector<double> o(numX,0);
  double molar_section=1.0/(NAv*p.epsilon_*c.h_*c.h_*1000.0);

  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      if (x>0)
        Jn=2.0*p.Domega_*(c.omega_T_[x-1]-c.omega_B_[x-1]-c.omega_T_[x]+c.omega_B_[x])/(c.dx_[x-1]+c.dx_[x]);
      if (x+1<numX)
        Jp=2.0*p.Domega_*(c.omega_T_[x+1]-c.omega_B_[x+1]-c.omega_T_[x]+c.omega_B_[x])/(c.dx_[x+1]+c.dx_[x]);

      double sig=0;
      for (unsigned k=0; k<numK; ++k)
        {
          sig+=p.ksig_omega_[k]*c.rho_[x][k];
        }

      sig*=molar_section/c.dx_[x];
      o[x]=(Jn+Jp)/c.dx_[x]+sig-p.kcat_omega_*c.omega_B_[x];

    }
  return o;

}

std::vector<double> SimplestModel::Omega_Bound(const SimplestModel::Param &p, const CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;

  std::vector<double> o(numX,0);

  double molar_section=1.0/(NAv*p.epsilon_*c.h_*c.h_)*1000.0;

  double K_omega=p.kcat_omega_/p.kon_omega_;
  for (unsigned x=0; x<numX; ++x)
    {
      double Mt=0;
      for (unsigned k=0; k<numK; ++k)
        Mt+=p.M_[k]*c.rho_[x][k];
      double R_T=Mt*molar_section/c.dx_[x];
      double b=R_T+c.omega_T_[x]+K_omega;
      o[x]=0.5*(b)-0.5*std::sqrt(b*b-4*R_T*c.omega_T_[x]);
    }
  return o;
}


std::vector<std::vector<double>> SimplestModel::
dRho_dt(const Param &p,
        const CortexState &c, bool hasOmega)const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();


  std::vector<std::vector<double>> dr(numX,std::vector<double>(numK,0));
  if (hasOmega)
    {
      for (unsigned x=0; x<numX; ++x)
        {
          double psi_F=c.psi_T_[x]-c.psi_B_[x];
          double omega_F=c.omega_T_[x]-c.omega_B_[x];
          for (unsigned k=0; k<numK; ++k)
            {
              double Jr,Jl,Ja;
              if (k+1<numK)
                Jr=p.g_left_[k+1]*c.rho_[x][k+1]
                    -(p.g_rigth_[k]
                      +p.g_max_psi_[k]*p.kon_psi_*psi_F
                      /(p.kcat_psi_+p.kon_psi_*psi_F)
                      +p.g_max_omega_[k]*p.kon_omega_*omega_F
                      /(p.kcat_omega_+p.kon_omega_*omega_F)
                      )*c.rho_[x][k];
              else
                Jr=0;
              if (k>0)
                Jl=-p.g_left_[k]*c.rho_[x][k]
                    +(p.g_rigth_[k-1]
                    +p.g_max_psi_[k-1]*p.kon_psi_*psi_F
                    /(p.kcat_psi_+p.kon_psi_*psi_F)
                    +p.g_max_omega_[k-1]*p.kon_omega_*omega_F
                    /(p.kcat_omega_+p.kon_omega_*omega_F)
                    )*c.rho_[x][k-1];
              else
                Jl=0;
              Ja= (+p.a_[k]
                   +p.a_psi_[k]*p.kon_psi_*psi_F
                   /(p.kcat_psi_+p.kon_psi_*psi_F)
                   +p.a_omega_[k]*p.kon_omega_*omega_F
                   /(p.kcat_omega_+p.kon_omega_*omega_F)
                   )*c.rho_[x][k];


              dr[x][k]=Jl+Jr-Ja;

            }
        }
    }
  else
    {
      for (unsigned x=0; x<numX; ++x)
        {
          for (unsigned k=0; k<numK; ++k)
            {
              double psi_F=c.psi_T_[x]-c.psi_B_[x];
              double Jr,Jl,Ja;
              if (k+1<numK)
                Jr=p.g_left_[k+1]*c.rho_[x][k+1]
                    -(p.g_rigth_[k]
                      +p.g_max_psi_[k]*p.kon_psi_*psi_F
                      /(p.kcat_psi_+p.kon_psi_*psi_F)
                      )*c.rho_[x][k];
              else
                Jr=0;
              if (k>0)
                Jl=-p.g_left_[k]*c.rho_[x][k]
                    +(p.g_rigth_[k-1]
                    +p.g_max_psi_[k-1]*p.kon_psi_*psi_F
                    /(p.kcat_psi_+p.kon_psi_*c.psi_T_[x])
                    )*c.rho_[x][k-1];
              else
                Jl=0;
              Ja= (p.a_[k]
                   +p.a_psi_[k]*p.kon_psi_*psi_F
                   /(p.kcat_psi_+p.kon_psi_*psi_F)
                   )*c.rho_[x][k];

              dr[x][k]=Jl+Jr-Ja;

            }
        }
    }

  return dr;
}





CortexSimulation SimplestModel::simulate(const Parameters& par,
                                         const SimplestModel::Param &p,
                                         const CortexExperiment &sp
                                         ,double dt)const
{


  std::cout<<"starts a Simulation on \n";
  sp.write(std::cout);
  par.write(std::cout);
  std::cout<<"dt ="<<dt<<"\n";


  CortexState c=init(p,sp);

  unsigned numSamples=std::ceil((sp.tsim_+sp.teq_)/sp.sample_time_)+1;


  CortexSimulation s(c,numSamples);
  s.p_=par;
  s.dt_=dt;
  double t=-sp.teq_;
  unsigned i=0;
  s.t_[i]=t;
  s.sdt_[i]=dt;
  s.omega_T_[i]=c.omega_T_;
  s.psi_T_[i]=c.psi_T_;
  s.rho_[i]=c.rho_;

  while (t+dt<0)
    {
      c=nextEuler(p,c,dt);
      t+=dt;
      if (t+sp.teq_>=sp.sample_time_*(i+1))
        {
          ++i;
          s.t_[i]=t;
          s.sdt_[i]=dt;
          s.omega_T_[i]=c.omega_T_;
          s.psi_T_[i]=c.psi_T_;
          s.omega_B_[i]=c.omega_B_;
          s.psi_B_[i]=c.psi_B_;
          s.rho_[i]=c.rho_;

        }
    }

  addDamp(c,p);
  while (t+dt<sp.tsim_)
    {
      c=nextEuler(p,c,dt);
      t+=dt;
      if (t+sp.teq_>=sp.sample_time_*(i+1))
        {
          ++i;
          s.t_[i]=t;
          s.sdt_[i]=dt;
          s.omega_T_[i]=c.omega_T_;
          s.psi_T_[i]=c.psi_T_;
          s.omega_B_[i]=c.omega_B_;
          s.psi_B_[i]=c.psi_B_;
          s.rho_[i]=c.rho_;

          if (i % 50 ==0)
            std::cerr<<"\t sample \t"<<i;
        }

    }
  return s;


}


BaseModel *BaseModel::create(const Parameters &p)
{
  double modelNumeber=p.get("model");
  auto it=models_.find(modelNumeber);
  if (it!=models_.end())
    {
      BaseModel*  out=it->second;
      out->loadParameters(p);
      return out;
    }
  else
    return nullptr;
}

std::map<double, BaseModel *> BaseModel::getModels()
{
  std::map<double, BaseModel *> o;
  o[Model00::number()]=new Model00;
  o[Model10::number()]=new Model10;

  return o;

}
std::map<double, BaseModel *> BaseModel::models_=getModels();

