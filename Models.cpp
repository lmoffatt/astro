#include "Models.h"

#include <iostream>




CortexState SimplestModel::init(const SimplestModel::Param &p, const CortexExperiment &s)const
{

  CortexState o(s.x_,s.dx_,s.h_,p.epsilon_,p.N_.size());
  double rM2,r0;
  if (p.N_.size()==9)
    {
      rM2=p.g_left_[2]/(p.g_rigth_[1]+p.g_left_[2]);
      r0=p.g_left_[4]/(p.g_rigth_[3]+p.g_left_[4]);
    }
  else
    {
      r0=p.g_left_[2]/(p.g_rigth_[1]+p.g_left_[2]);

    }
  unsigned numX=o.x_.size();
  if (p.N_.size()==9)
    {
      for (unsigned x=0; x<numX; ++x)
        {
          o.rho_[x][0]=p.dens_Neur_*std::pow(s.h_,2)*s.dx_[x]*1000;  // dens is in cells per liter
          o.rho_[x][1]=p.dens_Microglia_*std::pow(s.h_,2)*s.dx_[x]*1000*rM2;
          o.rho_[x][2]=p.dens_Microglia_*std::pow(s.h_,2)*s.dx_[x]*1000*(1-rM2);


          o.rho_[x][3]=p.dens_Astr_*std::pow(s.h_,2)*s.dx_[x]*1000*r0;
          o.rho_[x][4]=p.dens_Astr_*std::pow(s.h_,2)*s.dx_[x]*1000*(1-r0);
        }
    }
  else
    {
      for (unsigned x=0; x<numX; ++x)
        {
          o.rho_[x][0]=p.dens_Neur_*std::pow(s.h_,2)*s.dx_[x]*1000;  // dens is in cells per liter


          o.rho_[x][1]=p.dens_Astr_*std::pow(s.h_,2)*s.dx_[x]*1000*r0;
          o.rho_[x][2]=p.dens_Astr_*std::pow(s.h_,2)*s.dx_[x]*1000*(1-r0);
        }
    }
  o.psi_B_=Psi_Bound(p,o);
  o.omega_B_=Omega_Bound(p,o);
  return o;
}




CortexState SimplestModel::init(const SimplestModel::Param &p,
                                const Experiment &s
                                ,double dx)const
{
  std::vector<double> x_grid_in_m=s.x_in_m(dx);

  CortexState o(x_grid_in_m,s.h(),p.epsilon_,p.N_.size());

  double rM2,r0;
  if (p.N_.size()==9)
    {
      rM2=p.g_left_[2]/(p.g_rigth_[1]+p.g_left_[2]);
      r0=p.g_left_[4]/(p.g_rigth_[3]+p.g_left_[4]);
    }
  else
    {
      r0=p.g_left_[2]/(p.g_rigth_[1]+p.g_left_[2]);

    }

  unsigned numX=o.x_.size();
  if (p.N_.size()==9)
    {
      for (unsigned x=0; x<numX; ++x)
        {
          o.rho_[x][0]=p.dens_Neur_*std::pow(s.h(),2)*o.dx_[x]*1000;  // dens is in cells per liter
          o.rho_[x][1]=p.dens_Microglia_*std::pow(s.h(),2)*o.dx_[x]*1000*rM2;
          o.rho_[x][2]=p.dens_Microglia_*std::pow(s.h(),2)*o.dx_[x]*1000*(1-rM2);
          o.rho_[x][3]=p.dens_Astr_*std::pow(s.h(),2)*o.dx_[x]*1000*r0;
          o.rho_[x][4]=p.dens_Astr_*std::pow(s.h(),2)*o.dx_[x]*1000*(1-r0);
        }
    }
  else
    {
      for (unsigned x=0; x<numX; ++x)
        {

          o.rho_[x][0]=p.dens_Neur_*std::pow(s.h(),2)*o.dx_[x]*1000;  // dens is in cells per liter

          o.rho_[x][1]=p.dens_Astr_*std::pow(s.h(),2)*o.dx_[x]*1000*r0;
          o.rho_[x][2]=p.dens_Astr_*std::pow(s.h(),2)*o.dx_[x]*1000*(1-r0);

        }
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
      double Det=std::sqrt(b*b-4*R_T*c.psi_T_[x]);
      o[x]=2*R_T*c.psi_T_[x]/(b+Det);
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

  /// esto de abajo indica que para uso un espesor y ancho c.h_ para el numero rho de celulas
  double molar_section=1.0/(NAv*p.epsilon_*c.h_*c.h_*1000.0);

  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      double Dfn,Dfp;
      if (x>0)
        {
          Dfn=2.0*p.Domega_/(c.dx_[x-1]+c.dx_[x]);
          Jn=(c.omega_T_[x-1]-c.omega_B_[x-1]-c.omega_T_[x]+c.omega_B_[x]);
          Jn=Dfn*Jn;
        }
      if (x+1<numX)
        {
          Dfp=2.0*p.Domega_/(c.dx_[x+1]+c.dx_[x]);
          Jp=(c.omega_T_[x+1]-c.omega_B_[x+1]-c.omega_T_[x]+c.omega_B_[x]);
          Jp=Dfp*Jp;
        }
      double sig=0;
      double psi_F=c.psi_T_[x]-c.psi_B_[x];
      double omega_F=c.omega_T_[x]-c.omega_B_[x];

      if (p.DAMP_omega_ratio_==0)
        {for (unsigned k=0; k<numK; ++k)
            {
              sig+=(p.ksig_omega_[k]
                    +p.ksig_max_omega_[k]*p.kon_omega_*omega_F
                    /(p.kcat_omega_+p.kon_omega_*omega_F)
                    +p.ksig_max_psi_[k]*p.kon_psi_*psi_F
                    /(p.kcat_psi_+p.kon_psi_*psi_F)
                    )*c.rho_[x][k];

            }
        }

      /// esto de abajo indica que para uso un espesor y ancho c.h_ para el numero rho de celulas, el largo de la ///
      ///celda es variable, es dx_[x]
      // por otro lado
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
      double Det=std::sqrt(b*b-4*R_T*c.omega_T_[x]);
      o[x]=2*R_T*c.omega_T_[x]/(b+Det);
    }
  return o;
}













CortexSimulation SimplestModel::simulate(Parameters par,
                                         Param p,
                                         const Experiment &sp
                                         , double dx
                                        , double dtmin,
                                         std::size_t nPoints_per_decade,
                                         double dtmax,
                                         double tequilibrio
                                         )const
{
  /* std::cout<<"starts a Simulation on \n";
  sp.write(std::cout);
  par.write(std::cout);
  std::cout<<"dt ="<<dt<<"\n";
*/

  //  find inj_width information in parameters

  double dtprod=std::pow(10,1.0/nPoints_per_decade);



  CortexState c=init(p,sp,dx);
  CortexState cn=c;
  CortexSimulation s(c,sp.numSimPoints());
  s.p_=par;
  s.dt_=dtmax;
  double t=-tequilibrio;
  unsigned i=0;
  std::vector<std::vector<double>> dRho=c.rho_;

  while ((t+dtmax<=0)&&(c.isValid_))
    {
      nextEuler(cn,dRho,p,c,dtmax);
      c=cn;
      t+=dtmax;
    }
  if (!c.isValid_)
    {  s.isValid_=false;
      return s;
    }
  else
    {
      addDamp(c,p);
      cn=c;
      double dt_run=dtmin;
      while ((i<sp.numSimPoints())&&(c.isValid_))
        {
          t+=dt_run;

          nextEuler(cn,dRho,p,c,dt_run);
          c=cn;
          if (t>=sp.tSimul(i))
            {
              s.t_[i]=t;
              s.sdt_[i]=dt_run;
              s.omega_T_[i]=c.omega_T_;
              s.psi_T_[i]=c.psi_T_;
              s.omega_B_[i]=c.omega_B_;
              s.psi_B_[i]=c.psi_B_;
              s.rho_[i]=c.rho_;
              ++i;

            }

          if (dt_run<dtmax)
            {
              dt_run*=dtprod;
              if (dt_run>dtmax)
                dt_run=dtmax;
            }
          else
            dt_run=dtmax;

        }
      if (!c.isValid_)
        {  s.isValid_=false;
          return s;
        }
      else
        {
          s.isValid_=true;
          return s;
        }
    }
}





BaseModel *BaseModel::create(const Parameters &p)
{
  double modelNumeber=p.model();
  auto it=models_.find(modelNumeber);
  if (it!=models_.end())
    {
      BaseModel*  out=it->second->clone();
      out->loadParameters(p);
      return out;
    }
  else
    return nullptr;
}




std::vector<double> BaseModel::getObservedProbFromModel(const std::vector<double> &modelRho) const
{
  std::vector<double>  v(5);
  double n=modelRho[2]+modelRho[3]+modelRho[4]+modelRho[5]+modelRho[6];
  v[0]=modelRho[2]/n;
  v[1]=modelRho[3]/n;
  v[2]=modelRho[4]/n;
  v[3]=modelRho[5]/n;
  v[4]=modelRho[6]/n;
  return v;
}

std::vector<double> BaseModel::getObservedNumberFromModel(const std::vector<double> &modelRho) const
{
  std::vector<double>  v(5);
  v[0]=modelRho[2];
  v[1]=modelRho[3];
  v[2]=modelRho[4];
  v[3]=modelRho[5];
  v[4]=modelRho[6];
  return v;

}

std::vector<double> BaseModel::getNumberAtInjuryFromModel(const std::vector<double> &modelRho, double f) const
{
  std::vector<double>  v(modelRho);
  return v*f;

}

std::vector<double> BaseModel::getObservedNumberFromData(const std::vector<double> &modelRho) const
{
  std::vector<double>  v(5);
  v[0]=modelRho[0];
  v[1]=modelRho[1];
  v[2]=modelRho[2];
  v[3]=modelRho[3];
  v[4]=modelRho[4];
  return v;
}

std::size_t BaseModel::getNumberOfObservedStates()const
{
  return 5;
}


std::size_t BaseModel::getNumberOfSimulatedStatesAtInjury()const
{
  return 7;
}





std::map<double, BaseModel *> BaseModel::getModels()
{
  std::map<double, BaseModel *> o;
  o[Model00::number()]=new Model00;
  o[Model011::number()]=new Model011;
  o[Model012::number()]=new Model012;

  o[Model012_22::number()]=new Model012_22;

  o[Model012_51::number()]=new Model012_51;

  o[Model013::number()]=new Model013;

  o[Model013_23::number()]=new Model013_23;
  o[Model013_23_31::number()]=new Model013_23_31;

  o[Model013_51::number()]=new Model013_51;



  o[Model021::number()]=new Model021;
  o[Model022::number()]=new Model022;
  o[Model023::number()]=new Model023;
  o[Model031::number()]=new Model031;

  o[Model051::number()]=new Model051;

  o[Model10::number()]=new Model10;
  o[Model111::number()]=new Model111;

  o[Model112::number()]=new Model112;
  o[Model112_22::number()]=new Model112_22;

  o[Model112_22_31::number()]=new Model112_22_31;

  o[Model112_51::number()]=new Model112_51;

  o[Model112_52::number()]=new Model112_52;

  o[Model113::number()]=new Model113;

  o[Model113_42::number()]=new Model113_42;

  o[Model114::number()]=new Model114;

  o[Model114_24_32_44::number()]=new Model114_24_32_44;

  o[Model114_24_44::number()]=new Model114_24_44;


  o[Model115::number()]=new Model115;

  o[Model115_22::number()]=new Model115_22;

  o[Model115_25::number()]=new Model115_25;

  o[Model121::number()]=new Model121;
  o[Model122::number()]=new Model122;
  o[Model123::number()]=new Model123;
  o[Model124::number()]=new Model124;

  o[Model125::number()]=new Model125;


  o[Model131::number()]=new Model131;
  o[Model132::number()]=new Model132;
  o[Model141::number()]=new Model141;
  o[Model142::number()]=new Model142;
  o[Model144::number()]=new Model144;

  o[Model151::number()]=new Model151;

  o[Model152::number()]=new Model152;

  o[Model20::number()]=new Model20;

  o[Model211::number()]=new Model211;
  o[Model212::number()]=new Model212;
  o[Model213::number()]=new Model213;

  o[Model214::number()]=new Model214;
  o[Model215::number()]=new Model215;


  o[Model00m::number()]=new Model00m;
  o[Model013_23_31m::number()]=new Model013_23_31m;

  o[Model100m::number()]=new Model100m;
  o[Model114_24_32_44m::number()]=new Model114_24_32_44m;

  return o;

}
std::map<double, BaseModel *> BaseModel::models_=getModels();

