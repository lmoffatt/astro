#include "Models.h"

#include <iostream>




SimplestModel::SimplestModel() {}

std::pair<CortexState,double>
SimplestModel::nextEuler_Adapt_2
(const CortexState& one,
 const CortexState::Der &d,
 const SimplestModel::Param &p,
 const std::pair<CortexState,double> &c,
 double dt,
 double dtmin,
 double maxlogError,
 std::vector<double>& dts) const
{
  CortexState firstHalf(c.first);
  EulerStep(d,p,firstHalf,dt/2);
  auto d2=dStep(p,firstHalf);
  CortexState two(firstHalf);
  EulerStep(d2,p,two,dt-dt/2);

  double dist=one.maxlogdistance(two);
  if (dist<maxlogError || dt<dtmin*4.0)
    {
      dts.push_back(dt/2);
      dts.push_back(dt-dt/2);

      return {two,std::max(dist,c.second)};
    }
  else
    {
      auto fH=nextEuler_Adapt_2(firstHalf,d,p,c,dt/2,dtmin,maxlogError,dts);
      return nextEuler_Adapt(p,fH,dt-dt/2,dtmin,maxlogError,dts);
    }
}

std::pair<Di<CortexState>,double>
SimplestModel::nextEuler_Adapt_2_i
(CortexState one,
 const CortexState::Der &d,
 const Di<SimplestModel::Param> &p,
 const std::pair<Di<CortexState>, double> &c,
 double dt,
 double dtmin,
 double maxlogError) const
{
  CortexState firstHalf(c.first.first);
  EulerStep(d,p.first,firstHalf,dt/2);
  auto d2=dStep(p.first,firstHalf);
  CortexState two(firstHalf);
  EulerStep(d2,p.first,two,dt-dt/2);
  double dist=one.maxlogdistance(two);
  if (dist<maxlogError || dt<dtmin*4.0)
    {
      auto ones_i=nextEuler_i(p.second,c.first.second,dt);
      return {{one,ones_i},std::max(dist,c.second)};
    }
  else
    {
      auto fH=nextEuler_Adapt_2_i(firstHalf,d,p,c,dt/2,dtmin,maxlogError);
      return nextEuler_Adapt_i(p,fH,dt-dt/2,dtmin,maxlogError);
    }
}



std::pair<CortexState,double> SimplestModel::nextEuler_Adapt
(const SimplestModel::Param &p,
 const std::pair<CortexState,double> &c,
 double dt,
 double dtmin,
 double maxlogError,
 std::vector<double>& dts) const
{
  auto d=dStep(p,c.first);
  CortexState one(c.first);
  EulerStep(d,p,one,dt);
  return nextEuler_Adapt_2(one,d,p,c,dt,dtmin,maxlogError,dts);
}






std::pair<Di<CortexState>,double>
SimplestModel::nextEuler_Adapt_i
(const Di<SimplestModel::Param> &p,
 const std::pair<Di<CortexState>,double> &c,
 double dt,
 double dtmin,
 double maxlogError) const
{
  auto d=dStep(p.first,c.first.first);
  CortexState one(c.first.first);
  EulerStep(d,p.first,one,dt);
  return nextEuler_Adapt_2_i(one,d,p,c,dt,dtmin,maxlogError);

}


void SimplestModel::nextEuler
(const SimplestModel::Param &p,
 CortexState &c,
 double dt) const
{
  auto d=dStep(p,c);
  EulerStep(d,p,c,dt);
}

std::vector<CortexState> SimplestModel::nextEuler_i
(const std::vector<SimplestModel::Param> &p,
 const std::vector<CortexState> &c,
 double dt) const
{
  std::vector<CortexState> out(p.size());
  for (std::size_t i=0; i<p.size(); ++i)
    {
      auto d=dStep(p[i],c[i]);
      out[i]=c[i];
      EulerStep(d,p[i],out[i],dt);
    }
  return out;
}


CortexState::Der SimplestModel::dStep( const SimplestModel::Param &p, const CortexState &c) const
{
  CortexState::Der d;
  d.dPsi=dPsi_dt(p,c);
  d.hasOmega=p.Domega_>0;
  if (d.hasOmega)
    d.dOmega=dOmega_dt(p,c);
  d.dRho=dRho_dt(p,c,d.hasOmega);
  return d;
}



void SimplestModel::EulerStep(const CortexState::Der& d, const SimplestModel::Param &p,  CortexState &c, double dt)const
{
  bool hasOmega=p.Domega_>0;
  if (hasOmega)
    {
      if (!addStep(c.omega_T_,d.dOmega,dt))
        c.isValid_=false;
    }
  if (!addStep(c.psi_T_,d.dPsi,dt))
    c.isValid_=false;

  if (!addMStep(c.rho_,d.dRho,dt))
    c.isValid_=false;
  if (hasOmega)
    Omega_Bound(p,c);
  Psi_Bound(p,c);

}


std::pair<CortexState,std::pair<std::size_t,double>>
SimplestModel::CrankNicholsonStep(const CortexState::Der& d,
                                  const CortexState::dDer& dd,
                                  const SimplestModel::Param &p,
                                  const CortexState &c,
                                  double dt,
                                  double maxlogError,
                                  std::size_t maxloop)const
{

  std::vector<double> d0=d.toVector();
  auto dd1=dd.I_dt(dt);
  std::vector<double> d1=dd1.solve(d0);

  CortexState::Der d2(d);
  d2.fromVector(d1);

  CortexState c1=c;
  EulerStep(d2,p,c1,dt);
  CortexState ch(c);
  CortexState::Der d3=dStep(p,c1);
  EulerStep(d,p,ch,dt/2);
  CortexState c2(ch);
  EulerStep(d3,p,c2,dt-dt/2);
  double dist=c2.maxlogdistance(c1);
  std::size_t iter=0;

  while (dist>maxlogError &&  iter<maxloop)
    {
      std::swap(c1,c2);
      c2=ch;
      d3=dStep(p,c1);
      EulerStep(d3,p,c2,dt-dt/2);
      dist=c2.maxlogdistance(c1);
      ++iter;
    }

  return {c2,{iter,dist}};
}


std::pair<CortexState,std::pair<std::size_t,double>>
SimplestModel::CrankNicholsonStep(const CortexState::Der& d,
                                  const SimplestModel::Param &p,
                                  const CortexState &c,
                                  double dt,
                                  double maxlogError,
                                  std::size_t maxloop)const
{



  CortexState c1=c;
  EulerStep(d,p,c1,dt/2);
  CortexState::Der d02=dStep(p,c1);
  EulerStep(d02,p,c1,dt-dt/2);
  CortexState ch(c);
  CortexState::Der d3=dStep(p,c1);
  EulerStep(d,p,ch,dt/2);
  CortexState c2(ch);
  EulerStep(d3,p,c2,dt-dt/2);
  double dist=c2.maxlogdistance(c1);
  std::size_t iter=0;

  while (dist>maxlogError &&  iter<maxloop)
    {
      std::swap(c1,c2);
      c2=ch;
      d3=dStep(p,c1);
      EulerStep(d3,p,c2,dt-dt/2);
      dist=c2.maxlogdistance(c1);
      ++iter;
    }

  return {c2,{iter,dist}};
}




std::pair<CortexState,double>
SimplestModel::CrankNicholson_Adapt_2
(const CortexState& one,
 const CortexState::Der &d,
 const CortexState::dDer& dd,
 const SimplestModel::Param &p,
 const std::pair<CortexState,double> &c,
 double dt,
 double dtmin,
 double maxlogError,
 double maxlogErrorCN,
 std::size_t maxloop,
 std::pair<std::vector<double>,std::vector< std::size_t>>& dts) const
{
  auto firstHalf=CrankNicholsonStep(d,dd,p,c.first,dt/2,maxlogErrorCN,maxloop);
  auto d2=dStep(p,firstHalf.first);
  auto dd2=ddStep(p,c.first);
  auto two=CrankNicholsonStep(d2,dd2,p,firstHalf.first,dt-dt/2,maxlogErrorCN,maxloop);


  double dist=one.maxlogdistance(two.first);
  if (dist<maxlogError || dt<dtmin*4.0)
    {
      dts.first.push_back(dt/2);
      dts.second.push_back(firstHalf.second.second);
      dts.first.push_back(dt-dt/2);
      dts.second.push_back(two.second.second);
      return {two.first,std::max(dist,c.second)};
    }
  else
    {
      auto fH=CrankNicholson_Adapt_2(firstHalf.first,d,dd,p,c,dt/2,dtmin,maxlogError,maxlogErrorCN,maxloop,dts);
      return nextCrankNicholson_Adapt(p,fH,dt-dt/2,dtmin,maxlogError,maxlogErrorCN,dts,maxloop,true);
    }
}


std::pair<CortexState,double>
SimplestModel::CrankNicholson_Adapt_2
(const CortexState& one,
 const CortexState::Der &d,
 const SimplestModel::Param &p,
 const std::pair<CortexState,double> &c,
 double dt,
 double dtmin,
 double maxlogError,
 double maxlogErrorCN,
 std::size_t maxloop,
 std::pair<std::vector<double>,std::vector< std::size_t>>& dts) const
{
  auto firstHalf=CrankNicholsonStep(d,p,c.first,dt/2,maxlogErrorCN,maxloop);
  auto d2=dStep(p,firstHalf.first);
  auto two=CrankNicholsonStep(d2,p,firstHalf.first,dt-dt/2,maxlogErrorCN,maxloop);


  double dist=one.maxlogdistance(two.first);
  if (dist<maxlogError || dt<dtmin*4.0)
    {
      dts.first.push_back(dt/2);
      dts.second.push_back(firstHalf.second.second);
      dts.first.push_back(dt-dt/2);
      dts.second.push_back(two.second.second);
      return {two.first,std::max(dist,c.second)};
    }
  else
    {
      auto fH=CrankNicholson_Adapt_2(firstHalf.first,d,p,c,dt/2,dtmin,maxlogError,maxlogErrorCN,maxloop,dts);
      return nextCrankNicholson_Adapt(p,fH,dt-dt/2,dtmin,maxlogError,maxlogErrorCN,dts,maxloop,false);
    }
}



std::pair<CortexState,std::pair<double, std::size_t>> SimplestModel::nextCrankNicholson
(const SimplestModel::Param &p,
 const CortexState &c,
 double dt,
 double maxlogErrorCN,
 std::size_t maxloop,
 bool UseDerivative
 ) const
{
  auto d=dStep(p,c);
  if (UseDerivative)
    {
      auto dd=ddStep(p,c);
      return CrankNicholsonStep(d,dd,p,c,dt,maxlogErrorCN,maxloop);
    }
  else
    {
      return CrankNicholsonStep(d,p,c,dt,maxlogErrorCN,maxloop);

    }
}

std::pair<CortexState,double> SimplestModel::nextCrankNicholson_Adapt
(const SimplestModel::Param &p,
 const std::pair<CortexState,double> &c,
 double dt,
 double dtmin,
 double maxlogError,
 double maxlogErrorCN,
 std::pair<std::vector<double>,std::vector< std::size_t>>& dts,
 std::size_t maxloop,
 bool useDerivative
 ) const
{
  auto d=dStep(p,c.first);
  if (useDerivative)
    {
      auto dd=ddStep(p,c.first);
      auto one=CrankNicholsonStep(d,dd,p,c.first,dt,maxlogErrorCN,maxloop);
      return CrankNicholson_Adapt_2(one.first,d,dd,p,c,dt,dtmin,maxlogError,maxlogErrorCN,maxloop,dts);
    }
  else
    {
      auto one=CrankNicholsonStep(d,p,c.first,dt,maxlogErrorCN,maxloop);
      return CrankNicholson_Adapt_2(one.first,d,p,c,dt,dtmin,maxlogError,maxlogErrorCN,maxloop,dts);

    }
}


CortexState SimplestModel::init(const SimplestModel::Param &p, const CortexExperiment &s)const
{

  CortexState o(s.x_,s.dx_,s.h_,p.epsilon_,p.N_.size(),p.Domega_>0);
  double rM2=0,r0=0;
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
  Psi_Bound(p,o);
  if (o.hasOmega_)
    Omega_Bound(p,o);
  return o;
}

CortexSimulation SimplestModel::simulate(const Parameters &par, const SimplestModel::Param &p, const CortexExperiment &sp, double dt) const
{
  std::cout<<"starts a Simulation on \n";
  sp.write(std::cout);
  par.write(std::cout);
  std::cout<<"dt ="<<dt<<"\n";



  CortexState c=init(p,sp);
  CortexState cn=c;
  unsigned numSamples=std::ceil((sp.tsim_+sp.teq_)/sp.sample_time_)+1;


  CortexSimulation s(c,numSamples);
  s.p_=par;
  s.dt_=dt;
  s.h_=sp.h_;
  double t=-sp.teq_;
  unsigned i=0;
  s.t_[i]=t;
  s.maxlogErrt_[i]=dt;
  s.omega_T_[i]=c.omega_T_;
  s.psi_T_[i]=c.psi_T_;
  s.rho_[i]=c.rho_;

  std::vector<std::vector<double>> dRho=c.rho_;
  while (t+dt<0)
    {

      nextEuler(p,c,dt);
      if (!c.isValid_)
        {
          s.isValid_=false;
          return s;
        }
      t+=dt;
      if (t+sp.teq_>=sp.sample_time_*(i+1))
        {
          ++i;
          s.t_[i]=t;
          s.maxlogErrt_[i]=dt;
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
      nextEuler(p,c,dt);
      if (!c.isValid_)
        {
          s.isValid_=false;
          return s;
        }
      t+=dt;
      if (t+sp.teq_>=sp.sample_time_*(i+1))
        {
          ++i;
          s.t_[i]=t;
          s.maxlogErrt_[i]=dt;
          s.omega_T_[i]=c.omega_T_;
          s.psi_T_[i]=c.psi_T_;
          s.omega_B_[i]=c.omega_B_;
          s.psi_B_[i]=c.psi_B_;
          s.rho_[i]=c.rho_;

          if (i % 50 ==0)
            std::cerr<<"\t sample \t"<<i;
        }

    }
  s.isValid_=true;
  return s;


}




CortexState SimplestModel::init(const SimplestModel::Param &p,
                                const Experiment &s
                                ,double dx)const
{
  std::vector<double> x_grid_in_m=s.x_in_m(dx);

  CortexState o(x_grid_in_m,s.h(),p.epsilon_,p.N_.size(),p.Domega_>0);

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


  Psi_Bound(p,o);
  Omega_Bound(p,o);
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

void SimplestModel::Psi_Bound(const SimplestModel::Param &p, CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;


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
      c.psi_B_[x]=2*R_T*c.psi_T_[x]/(b+Det);
    }

}


std::vector<double> SimplestModel::dPsi_Bound_dPsi(const SimplestModel::Param &p, const CortexState &c) const
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
      o[x]=2*R_T/(b+Det)
          *(1
            -c.psi_T_[x]*(1+(c.psi_T_[x]+K_psi-R_T)/Det)
            /(b+Det)
            );
    }
  return o;
}

std::vector<std::vector<double>> SimplestModel::dPsi_Bound_dRho(const SimplestModel::Param &p, const CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;

  std::vector<std::vector<double>> o(numX,std::vector<double>(numK));

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
      for (unsigned k=0; k<numK; ++k)
        {
          double nk=molar_section/c.dx_[x]*p.N_[k];
          o[x][k]=2*nk*c.psi_T_[x]/(b+Det)*(1-R_T*(1+(R_T+K_psi-c.psi_T_[x])/Det)/(b+Det));
        }
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

void SimplestModel::Omega_Bound(const SimplestModel::Param &p,
                                CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;


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
      c.omega_B_[x]=2*R_T*c.omega_T_[x]/(b+Det);
    }
}



std::vector<double> SimplestModel::dOmega_Bound_dOmega(const SimplestModel::Param &p, const CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;

  std::vector<double> o(numX,0);

  double molar_section=1.0/(NAv*p.epsilon_*c.h_*c.h_*1000.0);

  double K_omega=p.kcat_omega_/p.kon_omega_;
  for (unsigned x=0; x<numX; ++x)
    {
      double Nt=0;
      for (unsigned k=0; k<numK; ++k)
        Nt+=p.N_[k]*c.rho_[x][k];
      double R_T=Nt*molar_section/c.dx_[x];
      double b=R_T+c.omega_T_[x]+K_omega;
      double Det=std::sqrt(b*b-4*R_T*c.omega_T_[x]);
      o[x]=2*R_T/(b+Det)*(1-c.omega_T_[x]*
                          (1+(c.omega_T_[x]+K_omega-R_T)/Det)/(b+Det));
    }
  return o;
}

std::vector<std::vector<double>> SimplestModel::domega_Bound_dRho(const SimplestModel::Param &p, const CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;

  std::vector<std::vector<double>> o(numX,std::vector<double>(numK));

  double molar_section=1.0/(NAv*p.epsilon_*c.h_*c.h_*1000.0);

  double K_omega=p.kcat_omega_/p.kon_omega_;
  for (unsigned x=0; x<numX; ++x)
    {
      double Nt=0;
      for (unsigned k=0; k<numK; ++k)
        Nt+=p.N_[k]*c.rho_[x][k];
      double R_T=Nt*molar_section/c.dx_[x];
      double b=R_T+c.omega_T_[x]+K_omega;
      double Det=std::sqrt(b*b-4*R_T*c.omega_T_[x]);
      for (unsigned k=0; k<numK; ++k)
        {
          double nk=molar_section/c.dx_[x]*p.N_[k];
          o[x][k]=2*nk*c.omega_T_[x]/(b+Det)*(1-R_T*(1+(R_T+K_omega-c.omega_T_[x])/Det)/(b+Det));
        }
    }
  return o;
}






CortexState::dDer SimplestModel::ddStep(const SimplestModel::Param& p ,const CortexState &c) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_[0].size();


  CortexState::dDer o(c.psi_B_.size(),c.rho_[0].size(),c.hasOmega_);

  if (c.hasOmega_)
    {
      auto dPsiBdPsi=dPsi_Bound_dPsi(p,c);
      auto dPsiBdRho=dPsi_Bound_dRho(p,c);

      auto dOmegaBdOmega=dOmega_Bound_dOmega(p,c);
      auto dOmegaBdRho=domega_Bound_dRho(p,c);

      for (unsigned x=0; x<numX; ++x)
        {
          double omega_F=c.omega_T_[x]-c.omega_B_[x];
          double psi_F=c.psi_T_[x]-c.psi_B_[x];

          double dLTdomega=0;
          double dLTdpsi=0;
          std::vector<double> dLTdrho(numK);
          for (std::size_t k=0; k<numK; ++k)
            {
              dLTdomega+=p.ksig_max_omega_[k]*p.kon_omega_*p.kcat_omega_/
                  sqr(p.kcat_omega_+p.kon_omega_*(omega_F))
                  *c.rho_[x][k];
              dLTdpsi+=p.ksig_max_psi_[k]*p.kon_psi_*p.kcat_psi_/
                  sqr(p.kcat_psi_+p.kon_psi_*(psi_F))
                  *c.rho_[x][k];
              dLTdrho[k]=p.ksig_omega_[k]
                  +p.ksig_max_omega_[k]*p.kon_omega_*omega_F
                  /(p.kcat_omega_+p.kon_omega_*omega_F)
                  +p.ksig_max_psi_[k]*p.kon_psi_*psi_F
                  /(p.kcat_psi_+p.kon_psi_*psi_F)

                  -p.ksig_max_omega_[k]*p.kon_omega_*p.kcat_omega_/
                  sqr(p.kcat_omega_+p.kon_omega_*(omega_F))
                  *c.rho_[x][k]*dOmegaBdRho[x][k]

                  -p.ksig_max_psi_[k]*p.kon_psi_*p.kcat_psi_/
                  sqr(p.kcat_psi_+p.kon_psi_*(psi_F))
                  *c.rho_[x][k]*dPsiBdRho[x][k];
            }
          dLTdomega*=(1.0-dOmegaBdOmega[x]);
          dLTdpsi*=(1.0-dPsiBdPsi[x]);

          if (x>0)
            {
              o.dpsi_dpsi(x,x-1)=2.0*p.Dpsi_/c.dx_[x]
                  /(c.dx_[x-1]+c.dx_[x])
                  *(1-dPsiBdPsi[x-1]);
              if (c.hasOmega_)
                o.domega_domega(x,x-1)=2.0*p.Domega_/c.dx_[x]
                    /(c.dx_[x-1]+c.dx_[x])
                    *(1-dOmegaBdOmega[x-1]);
              o.domega_dpsi(x,x)=dLTdpsi;
              for (std::size_t k=0; k<numK; ++k)
                {
                  o.dpsi_drho(x,x-1,k)=-2.0*p.Dpsi_/c.dx_[x]
                      /(c.dx_[x-1]+c.dx_[x])
                      *(dPsiBdRho[x-1][k]);
                  o.domega_drho(x,x-1,k)=-2.0*p.Domega_/c.dx_[x]
                      /(c.dx_[x-1]+c.dx_[x])
                      *(dOmegaBdRho[x-1][k]);
                }
            }
          if (x+1<numX)
            {
              o.dpsi_dpsi(x,x+1)=2.0*p.Dpsi_/c.dx_[x]
                  /(c.dx_[x+1]+c.dx_[x])
                  *(1-dPsiBdPsi[x+1]);
              o.domega_domega(x,x+1)=2.0*p.Domega_/c.dx_[x]
                  /(c.dx_[x+1]+c.dx_[x])
                  *(1-dOmegaBdOmega[x+1])
                  +dLTdomega;
              for (std::size_t k=0; k<numK; ++k)
                {
                  o.dpsi_drho(x,x+1,k)=2.0*p.Dpsi_/c.dx_[x]
                      /(c.dx_[x+1]+c.dx_[x])
                      *(-dPsiBdRho[x+1][k]);
                  o.domega_drho(x,x+1,k)=2.0*p.Domega_/c.dx_[x]
                      /(c.dx_[x+1]+c.dx_[x])
                      *(-dOmegaBdRho[x+1][k]);
                }
            }
          o.dpsi_dpsi(x,x)=2.0*p.Dpsi_/c.dx_[x]*
              (1.0/(c.dx_[x-1]+c.dx_[x])+1.0/(c.dx_[x+1]+c.dx_[x]))*
              (1-dPsiBdPsi[x])-
              p.kcat_psi_*dPsiBdPsi[x];

          o.domega_domega(x,x)=2.0*p.Domega_/c.dx_[x]*
              (1.0/(c.dx_[x-1]+c.dx_[x])+1.0/(c.dx_[x+1]+c.dx_[x]))
              *(1-dOmegaBdOmega[x])-
              p.kcat_omega_*dOmegaBdOmega[x]+
              dLTdomega;

          o.domega_dpsi(x,x)=dLTdpsi;

          for (std::size_t k=0; k<numK; ++k)
            {
              o.dpsi_drho(x,x,k)=2.0*p.Dpsi_/c.dx_[x]
                  *(1.0/(c.dx_[x-1]+c.dx_[x])+1.0/(c.dx_[x+1]+c.dx_[x]))
                  *(-dPsiBdRho[x][k])
                  -p.kcat_psi_*dPsiBdRho[x][k];
              o.domega_drho(x,x,k)=2.0*p.Domega_/c.dx_[x]
                  /(c.dx_[x+1]+c.dx_[x])
                  *(-dOmegaBdRho[x][k])
                  -p.kcat_omega_*dOmegaBdRho[x][k]
                  +dLTdrho[k];

            }

          for (unsigned k=0; k<numK; ++k)
            {

              if (k+1<numK)
                {
                  o.drho_drho(x,k,x,k)+= -p.g_rigth_[k]
                      -p.g_max_psi_[k]*psi_F
                      -(p.Keq_gmax_psi_[k]+psi_F)
                      +p.g_max_omega_[k]*omega_F
                      /(p.Keq_gmax_omega_[k]+omega_F)
                      //dpsidrho
                      -p.g_max_psi_[k] *p.Keq_gmax_psi_[k]
                      /sqr(p.Keq_gmax_psi_[k]+psi_F)
                      *(-dPsiBdRho[x][k])*c.rho_[x][k]
                      //domegdrho
                      -p.g_max_omega_[k] *p.Keq_gmax_omega_[k]
                      /sqr(p.Keq_gmax_omega_[k]+omega_F)
                      *(-dOmegaBdRho[x][k])*c.rho_[x][k]
                      ;
                  o.drho_drho(x,k,x,k+1)+=p.g_left_[k+1];
                  o.drho_dpsi(x,k,x)+=
                      -p.g_max_psi_[k] *p.Keq_gmax_psi_[k]
                      /sqr(p.Keq_gmax_psi_[k]+psi_F)
                      *(1-dPsiBdPsi[x])*c.rho_[x][k]
                      ;
                  o.drho_domega(x,k,x)+=
                      -p.g_max_omega_[k] *p.Keq_gmax_omega_[k]
                      /sqr(p.Keq_gmax_omega_[k]+omega_F)
                      *(1-dOmegaBdOmega[x])*c.rho_[x][k]
                      ;

                }
              if (k>0)
                {
                  o.drho_drho(x,k,x,k)+= -p.g_left_[k];
                  o.drho_drho(x,k,x,k-1)+= -p.g_rigth_[k-1]
                      -p.g_max_psi_[k-1]*psi_F
                      -(p.Keq_gmax_psi_[k-1]+psi_F)
                      +p.g_max_omega_[k-1]*omega_F
                      /(p.Keq_gmax_omega_[k-1]+omega_F)
                      //dpsidrho
                      -p.g_max_psi_[k-1] *p.Keq_gmax_psi_[k-1]
                      /sqr(p.Keq_gmax_psi_[k-1]+psi_F)
                      *(-dPsiBdRho[x][k-1])*c.rho_[x][k-1]
                      //domegdrho
                      -p.g_max_omega_[k-1] *p.Keq_gmax_omega_[k-1]
                      /sqr(p.Keq_gmax_omega_[k-1]+omega_F)
                      *(-dOmegaBdRho[x][k-1])*c.rho_[x][k-1]
                      ;
                  o.drho_dpsi(x,k,x)+=
                      -p.g_max_psi_[k-1] *p.Keq_gmax_psi_[k-1]
                      /sqr(p.Keq_gmax_psi_[k-1]+psi_F)
                      *(1-dPsiBdPsi[x])*c.rho_[x][k-1]
                      ;
                  o.drho_domega(x,k,x)+=
                      -p.g_max_omega_[k-1] *p.Keq_gmax_omega_[k-1]
                      /sqr(p.Keq_gmax_omega_[k-1]+omega_F)
                      *(1-dOmegaBdOmega[x])*c.rho_[x][k-1]
                      ;

                }
              o.drho_drho(x,k,x,k)+=-p.a_[k]
                  -p.a_psi_[k]*p.kon_psi_*psi_F
                  /(p.kcat_psi_+p.kon_psi_*psi_F)
                  -p.a_omega_[k]*p.kon_omega_*omega_F
                  /(p.kcat_omega_+p.kon_omega_*omega_F)

                  -p.a_psi_[k]*p.kon_psi_*p.kcat_psi_
                  /sqr(p.kcat_psi_+p.kon_psi_*psi_F)
                  *dPsiBdRho[x][k]*c.rho_[x][k]

                  -p.a_omega_[k]*p.kon_omega_*p.kcat_omega_
                  /sqr(p.kcat_omega_+p.kon_omega_*omega_F)
                  *dOmegaBdRho[x][k]*c.rho_[x][k]
                  ;
              o.drho_dpsi(x,k,x)+=
                  -p.a_psi_[k]*p.kon_psi_*p.kcat_psi_
                  /sqr(p.kcat_psi_+p.kon_psi_*psi_F)
                  *dPsiBdPsi[x]*c.rho_[x][k]
                  ;
              o.drho_domega(x,k,x)+=
                  -p.a_omega_[k]*p.kon_omega_*p.kcat_omega_
                  /sqr(p.kcat_omega_+p.kon_omega_*omega_F)
                  *dOmegaBdOmega[x]*c.rho_[x][k]
                  ;


            }

        }


      return o;
    }
  else
    {
      auto dPsiBdPsi=dPsi_Bound_dPsi(p,c);
      auto dPsiBdRho=dPsi_Bound_dRho(p,c);


      for (unsigned x=0; x<numX; ++x)
        {
          double psi_F=c.psi_T_[x]-c.psi_B_[x];

          double dLTdpsi=0;
          std::vector<double> dLTdrho(numK);
          for (std::size_t k=0; k<numK; ++k)
            {
              dLTdpsi+=p.ksig_max_psi_[k]*p.kon_psi_*p.kcat_psi_/
                  sqr(p.kon_psi_*psi_F)
                  *c.rho_[x][k];
              dLTdrho[k]=p.ksig_omega_[k]
                  +p.ksig_max_psi_[k]*p.kon_psi_*psi_F
                  /(p.kcat_psi_+p.kon_psi_*psi_F)
                  -p.ksig_max_psi_[k]*p.kon_psi_*p.kcat_psi_/
                  sqr(p.kon_psi_*psi_F)
                  *c.rho_[x][k]*dPsiBdRho[x][k];
            }
          dLTdpsi*=(1.0-dPsiBdPsi[x]);

          if (x>0)
            {
              o.dpsi_dpsi(x,x-1)=2.0*p.Dpsi_/c.dx_[x]
                  /(c.dx_[x-1]+c.dx_[x])
                  *(1-dPsiBdPsi[x-1]);
              for (std::size_t k=0; k<numK; ++k)
                {
                  o.dpsi_drho(x,x-1,k)=-2.0*p.Dpsi_/c.dx_[x]
                      /(c.dx_[x-1]+c.dx_[x])
                      *(dPsiBdRho[x-1][k]);
                }
              o.dpsi_dpsi(x,x)+=2.0*p.Dpsi_/c.dx_[x]*
                  (1.0/(c.dx_[x-1]+c.dx_[x]))*
                  (1-dPsiBdPsi[x])-
                  p.kcat_psi_*dPsiBdPsi[x];

              for (std::size_t k=0; k<numK; ++k)
                {
                  o.dpsi_drho(x,x,k)+=2.0*p.Dpsi_/c.dx_[x]
                      *(1.0/(c.dx_[x-1]+c.dx_[x]))
                      *(-dPsiBdRho[x][k])
                      -p.kcat_psi_*dPsiBdRho[x][k];

                }

            }
          if (x+1<numX)
            {
              o.dpsi_dpsi(x,x+1)=2.0*p.Dpsi_/c.dx_[x]
                  /(c.dx_[x+1]+c.dx_[x])
                  *(1-dPsiBdPsi[x+1]);
              for (std::size_t k=0; k<numK; ++k)
                {
                  o.dpsi_drho(x,x+1,k)=2.0*p.Dpsi_/c.dx_[x]
                      /(c.dx_[x+1]+c.dx_[x])
                      *(-dPsiBdRho[x+1][k]);
                }
              o.dpsi_dpsi(x,x)+=2.0*p.Dpsi_/c.dx_[x]*
                  (+1.0/(c.dx_[x+1]+c.dx_[x]))*
                  (1-dPsiBdPsi[x]);
              for (std::size_t k=0; k<numK; ++k)
                {
                  o.dpsi_drho(x,x,k)+=2.0*p.Dpsi_/c.dx_[x]
                      *(1.0/(c.dx_[x+1]+c.dx_[x]))
                      *(-dPsiBdRho[x][k]);

                }

            }
          o.dpsi_dpsi(x,x)+=- p.kcat_psi_*dPsiBdPsi[x];

          for (std::size_t k=0; k<numK; ++k)
            {
              o.dpsi_drho(x,x,k)+=-p.kcat_psi_*dPsiBdRho[x][k];

            }


          for (unsigned k=0; k<numK; ++k)
            {

              if (k+1<numK)
                {
                  o.drho_drho(x,k,x,k)+= -p.g_rigth_[k]
                      -p.g_max_psi_[k]*psi_F
                      -(p.Keq_gmax_psi_[k]+psi_F)
                      //dpsidrho
                      -p.g_max_psi_[k] *p.Keq_gmax_psi_[k]
                      /sqr(p.Keq_gmax_psi_[k]+psi_F)
                      *(-dPsiBdRho[x][k])*c.rho_[x][k]
                      //domegdrho
                      ;
                  o.drho_drho(x,k,x,k+1)+=p.g_left_[k+1];
                  o.drho_dpsi(x,k,x)+=
                      -p.g_max_psi_[k] *p.Keq_gmax_psi_[k]
                      /sqr(p.Keq_gmax_psi_[k]+psi_F)
                      *(1-dPsiBdPsi[x])*c.rho_[x][k]
                      ;

                }
              if (k>0)
                {
                  o.drho_drho(x,k,x,k)+= -p.g_left_[k];
                  o.drho_drho(x,k,x,k-1)+= -p.g_rigth_[k-1]
                      -p.g_max_psi_[k-1]*psi_F
                      -(p.Keq_gmax_psi_[k-1]+psi_F)
                      //dpsidrho
                      -p.g_max_psi_[k-1] *p.Keq_gmax_psi_[k-1]
                      /sqr(p.Keq_gmax_psi_[k-1]+psi_F)
                      *(-dPsiBdRho[x][k-1])*c.rho_[x][k-1]
                      //domegdrho
                      ;
                  o.drho_dpsi(x,k,x)+=
                      -p.g_max_psi_[k-1] *p.Keq_gmax_psi_[k-1]
                      /sqr(p.Keq_gmax_psi_[k-1]+psi_F)
                      *(1-dPsiBdPsi[x])*c.rho_[x][k-1]
                      ;

                }
              o.drho_drho(x,k,x,k)+=-p.a_[k]
                  -p.a_psi_[k]*p.kon_psi_*psi_F
                  /(p.kcat_psi_+p.kon_psi_*psi_F)

                  -p.a_psi_[k]*p.kon_psi_*p.kcat_psi_
                  /sqr(p.kcat_psi_+p.kon_psi_*psi_F)
                  *dPsiBdRho[x][k]*c.rho_[x][k]

                  ;
              o.drho_dpsi(x,k,x)+=
                  -p.a_psi_[k]*p.kon_psi_*p.kcat_psi_
                  /sqr(p.kcat_psi_+p.kon_psi_*psi_F)
                  *dPsiBdPsi[x]*c.rho_[x][k]
                  ;


            }

        }


      return o;

    }
}
std::vector<std::vector<double> >
SimplestModel::dRho_dt(const SimplestModel::Param &p, const CortexState &c, bool hasOmega) const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  std::vector<std::vector<double> > drho_(numX,std::vector<double>(numK,0.0));


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
                      +p.g_max_psi_[k]*psi_F
                      /(p.Keq_gmax_psi_[k]+psi_F)
                      +p.g_max_omega_[k]*omega_F
                      /(p.Keq_gmax_omega_[k]+omega_F)
                      )*c.rho_[x][k];
              else
                Jr=0;
              if (k>0)
                Jl=-p.g_left_[k]*c.rho_[x][k]
                    +(p.g_rigth_[k-1]
                    +p.g_max_psi_[k-1]*psi_F
                    /(p.Keq_gmax_psi_[k-1]+psi_F)
                    +p.g_max_omega_[k-1]*omega_F
                    /(p.Keq_gmax_omega_[k-1]+omega_F)
                    )*c.rho_[x][k-1];
              else
                Jl=0;
              Ja= (+p.a_[k]
                   +p.a_psi_[k]*p.kon_psi_*psi_F
                   /(p.kcat_psi_+p.kon_psi_*psi_F)
                   +p.a_omega_[k]*p.kon_omega_*omega_F
                   /(p.kcat_omega_+p.kon_omega_*omega_F)
                   )*c.rho_[x][k];


              drho_[x][k]=Jl+Jr-Ja;

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
                      +p.g_max_psi_[k]*psi_F
                      /(p.Keq_gmax_psi_[k]+psi_F)
                      )*c.rho_[x][k];
              else
                Jr=0;
              if (k>0)
                Jl=-p.g_left_[k]*c.rho_[x][k]
                    +(p.g_rigth_[k-1]
                    +p.g_max_psi_[k-1]*psi_F
                    /(p.Keq_gmax_psi_[k-1]+psi_F)
                    )*c.rho_[x][k-1];
              else
                Jl=0;
              Ja= (p.a_[k]
                   +p.a_psi_[k]*p.kon_psi_*psi_F
                   /(p.kcat_psi_+p.kon_psi_*psi_F)
                   )*c.rho_[x][k];

              drho_[x][k]=Jl+Jr-Ja;

            }
        }
    }

  return drho_;
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
  CortexSimulation s(c,sp.numSimPoints());
  s.p_=par;
  s.dt_=dtmax;
  s.h_=sp.h();
  double t=-tequilibrio;
  unsigned i=0;
  std::vector<std::vector<double>> dRho=c.rho_;
  double dt_run=dtmin;

  while ((t+dt_run<=0)&&(c.isValid_))
    {
      nextEuler(p,c,dtmax);
      t+=dt_run;
      if (dt_run<-t)
        {
          dt_run*=dtprod;
          if (dt_run>-t)
            dt_run=-t;
        }
    }
  if (!c.isValid_)
    {  s.isValid_=false;
      return s;
    }
  else
    {
      addDamp(c,p);
      dt_run=dtmin;
      while ((i<sp.numSimPoints())&&(c.isValid_))
        {
          t+=dt_run;

          nextEuler(p,c,dt_run);
          if (t>=sp.tSimul(i))
            {
              s.t_[i]=t;
              s.maxlogErrt_[i]=dt_run;
              s.psi_T_[i]=c.psi_T_;
              s.psi_B_[i]=c.psi_B_;
              if (c.hasOmega_)
                {
                  s.omega_T_[i]=c.omega_T_;
                  s.omega_B_[i]=c.omega_B_;
                }
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

std::pair<CortexSimulation, std::vector<double>> SimplestModel::simulate_Adapted(
    Parameters par,
    SimplestModel::Param p,
    const Experiment &sp,
    double dx,
    double dtmin,
    std::size_t nPoints_per_decade,
    double dtmax,
    double tequilibrio,
    double maxLogError,double dtinfratio) const

{
  double dtprod=std::pow(10,1.0/nPoints_per_decade);
  std::pair<CortexState,double> c{init(p,sp,dx),0};
  CortexSimulation s(c.first,sp.numSimPoints());
  s.p_=par;
  s.dt_=dtmax;
  s.h_=sp.h();
  double t=-tequilibrio;
  unsigned i=0;
  std::size_t n=std::ceil(tequilibrio/dtmax);
  double dtrun=tequilibrio/n;
  std::vector<double> dts;
  while (t+dtrun<=0)
    {
      c=nextEuler_Adapt(p,c,dtrun,dtinfratio*dtrun,maxLogError,dts);
      t+=dtrun;
    }
  if (!std::isfinite(c.second))
    {
      s.isValid_=false;
      return {s,dts};
    }
  else
    {
      double dt_run=dtmin;
      addDamp(c.first,p);
      while ((c.first.isValid_)&&(std::isfinite(c.second)))
        {
          if (t>=sp.tSimul(i)&&(std::isfinite(c.second)))
            {
              s.t_[i]=t;
              s.maxlogErrt_[i]=c.second;
              s.omega_T_[i]=c.first.omega_T_;
              s.psi_T_[i]=c.first.psi_T_;
              s.omega_B_[i]=c.first.omega_B_;
              s.psi_B_[i]=c.first.psi_B_;
              s.rho_[i]=c.first.rho_;
              c.second=0;
              ++i;
              if(i>=sp.numSimPoints())
                {
                  if (!c.first.isValid_||(!std::isfinite(c.second)))
                    {
                      s.isValid_=false;
                      return {s,dts};
                    }
                  else
                    {
                      s.isValid_=true;
                      return {s,dts};
                    }
                }


            }
          if (t+dt_run>sp.tSimul(i))
            dt_run=sp.tSimul(i)-t;
          else if (dt_run<dtmax)
            {
              dt_run*=dtprod;
              if (dt_run>dtmax)
                dt_run=dtmax;
            }
          else
            dt_run=dtmax;
          t+=dt_run;
          c=nextEuler_Adapt(p,c,dt_run,dtinfratio*dt_run,maxLogError,dts);

        }
      s.isValid_=false;
      return {s,dts};
    }
}


std::pair<CortexSimulation, std::pair<std::vector<double>,std::vector< std::size_t>>> SimplestModel::simulate_CrankNicholson_Adapted(
    Parameters par,
    SimplestModel::Param p,
    const Experiment &sp,
    double dx,
    double dtmin,
    std::size_t nPoints_per_decade,
    double dtmax,
    double tequilibrio,
    double maxLogError,
    double maxlogErrorCN,
    double dtinfratio,
    std::size_t maxloop,
    bool UseDerivative) const

{
  double dtprod=std::pow(10,1.0/nPoints_per_decade);
  std::pair<CortexState,double> c{init(p,sp,dx),0};
  CortexSimulation s(c.first,sp.numSimPoints());
  s.p_=par;
  s.dt_=dtmax;
  s.h_=sp.h();
  double t=-tequilibrio;
  unsigned i=0;
  std::size_t n=std::ceil(tequilibrio/dtmax);
  double dtrun=tequilibrio/n;
  std::pair<std::vector<double>,std::vector< std::size_t>>  dts;
  while (t+dtrun<=0)
    {
      c=nextCrankNicholson_Adapt(p,c,dtrun,dtinfratio*dtrun,maxLogError,maxlogErrorCN,dts,maxloop,UseDerivative);
      t+=dtrun;
    }
  if (!std::isfinite(c.second))
    {
      s.isValid_=false;
      return {s,dts};
    }
  else
    {
      double dt_run=dtmin;
      addDamp(c.first,p);
      while ((c.first.isValid_)&&(std::isfinite(c.second)))
        {
          if (t>=sp.tSimul(i)&&(std::isfinite(c.second)))
            {
              s.t_[i]=t;
              s.maxlogErrt_[i]=c.second;
              s.omega_T_[i]=c.first.omega_T_;
              s.psi_T_[i]=c.first.psi_T_;
              s.omega_B_[i]=c.first.omega_B_;
              s.psi_B_[i]=c.first.psi_B_;
              s.rho_[i]=c.first.rho_;
              c.second=0;
              ++i;
              if(i>=sp.numSimPoints())
                {
                  if (!c.first.isValid_||(!std::isfinite(c.second)))
                    {
                      s.isValid_=false;
                      return {s,dts};
                    }
                  else
                    {
                      s.isValid_=true;
                      return {s,dts};
                    }
                }


            }
          if (t+dt_run>sp.tSimul(i))
            dt_run=sp.tSimul(i)-t;
          else if (dt_run<dtmax)
            {
              dt_run*=dtprod;
              if (dt_run>dtmax)
                dt_run=dtmax;
            }
          else
            dt_run=dtmax;
          t+=dt_run;
          c=nextCrankNicholson_Adapt(p,c,dtrun,dtinfratio*dt_run,maxLogError,maxlogErrorCN,dts,maxloop,UseDerivative);

        }
      s.isValid_=false;
      return {s,dts};
    }
}

CortexSimulation SimplestModel::simulate_CrankNicholson(
    Parameters par,
    SimplestModel::Param p,
    const Experiment &sp,
    double dx,
    double dtmin,
    std::size_t nPoints_per_decade,
    double dtmax,
    double tequilibrio,
    double maxlogErrorCN,
    std::size_t maxloop,
    bool UseDerivative) const

{
  double dtprod=std::pow(10,1.0/nPoints_per_decade);
  std::pair<CortexState,std::pair<std::size_t,double>> c{init(p,sp,dx),{0,0}};
  CortexSimulation s(c.first,sp.numSimPoints());
  s.p_=par;
  s.dt_=dtmax;
  s.h_=sp.h();
  double t=-tequilibrio;
  unsigned i=0;
  std::size_t n=std::ceil(tequilibrio/dtmax);
  double dtrun=tequilibrio/n;
  while (t+dtrun<=0)
    {
      c=nextCrankNicholson(p,c.first,dtrun,maxlogErrorCN,maxloop,UseDerivative);
      t+=dtrun;
    }
  if (!std::isfinite(c.second.second))
    {
      s.isValid_=false;
      return s;
    }
  else
    {
      double dt_run=dtmin;
      addDamp(c.first,p);
      while ((c.first.isValid_)&&(std::isfinite(c.second.second)))
        {
          if (t>=sp.tSimul(i)&&(std::isfinite(c.second.second)))
            {
              s.t_[i]=t;
              s.maxlogErrt_[i]=c.second.second;
              s.omega_T_[i]=c.first.omega_T_;
              s.psi_T_[i]=c.first.psi_T_;
              s.omega_B_[i]=c.first.omega_B_;
              s.psi_B_[i]=c.first.psi_B_;
              s.rho_[i]=c.first.rho_;
              c.second.second=0;
              ++i;
              if(i>=sp.numSimPoints())
                {
                  if (!c.first.isValid_||(!std::isfinite(c.second.second)))
                    {
                      s.isValid_=false;
                      return s;
                    }
                  else
                    {
                      s.isValid_=true;
                      return s;
                    }
                }


            }
          if (t+dt_run>sp.tSimul(i))
            dt_run=sp.tSimul(i)-t;
          else if (dt_run<dtmax)
            {
              dt_run*=dtprod;
              if (dt_run>dtmax)
                dt_run=dtmax;
            }
          else
            dt_run=dtmax;
          t+=dt_run;
          c=nextCrankNicholson(p,c.first,dtrun,maxlogErrorCN,maxloop,UseDerivative);

        }
      s.isValid_=false;
      return s;
    }
}



CortexSimulation SimplestModel::simulate(Parameters par,
                                         Param p,
                                         const Experiment &sp
                                         , double dx
                                         , const std::pair<std::vector<double>,std::vector<std::size_t>>& dts,
                                         double tequilibrio
                                         )const
{
  CortexState c{init(p,sp,dx)};
  CortexSimulation s(c,sp.numSimPoints());
  s.p_=par;
  s.dt_=dts.first[0];
  s.h_=sp.h();
  double t=-tequilibrio;
  unsigned i=0;

  std::size_t its=0;
  double dt_run=dts.first[0];
  while ((t<0)&&(c.isValid_))
    {
      dt_run=dts.first[its];
      nextEuler(p,c,dt_run);
      t+=dt_run;
      ++its;

    }
  if (!c.isValid_)
    {
      s.isValid_=false;
      return s;
    }
  else
    {
      addDamp(c,p);
      while (c.isValid_)
        {
          if (t+dt_run*std::sqrt(std::numeric_limits<double>::epsilon())>=sp.tSimul(i)&&(c.isValid_))
            {
              s.t_[i]=sp.tSimul(i);
              s.maxlogErrt_[i]=dt_run;
              s.omega_T_[i]=c.omega_T_;
              s.psi_T_[i]=c.psi_T_;
              s.omega_B_[i]=c.omega_B_;
              s.psi_B_[i]=c.psi_B_;
              s.rho_[i]=c.rho_;
              ++i;
              if(i>=sp.numSimPoints())
                {
                  s.isValid_=c.isValid_;
                  return s;
                }



            }
          dt_run=dts.first[its];
          nextEuler(p,c,dt_run);
          t+=dt_run;
          ++its;

        }
      s.isValid_=false;
      return s;
    }
}






CortexSimulation SimplestModel::simulate_CrankNicholson_dt(Parameters par,
                                                           Param p,
                                                           const Experiment &sp
                                                           , double dx
                                                           , const std::vector<std::pair<double,std::size_t>>& dts,
                                                           double tequilibrio,
                                                           bool UseDerivative
                                                           )const
{
  std::pair<CortexState,std::pair<std::size_t,double>>  c{init(p,sp,dx),{0,0}};
  CortexSimulation s(c.first,sp.numSimPoints());
  s.p_=par;
  s.dt_=dts[0].first;
  s.h_=sp.h();
  double t=-tequilibrio;
  unsigned i=0;

  std::size_t its=0;
  double dt_run=dts[0].first;
  while ((t<0)&&(c.first.isValid_))
    {
      dt_run=dts[its].first;
      c=nextCrankNicholson(p,c.first,dt_run,0,dts[0].second,UseDerivative);
      t+=dt_run;
      ++its;

    }
  if (!c.first.isValid_)
    {
      s.isValid_=false;
      return s;
    }
  else
    {
      addDamp(c.first,p);
      while (c.first.isValid_)
        {
          if (t+dt_run*std::sqrt(std::numeric_limits<double>::epsilon())>=sp.tSimul(i)&&(c.first.isValid_))
            {
              s.t_[i]=sp.tSimul(i);
              s.maxlogErrt_[i]=dt_run;
              s.omega_T_[i]=c.first.omega_T_;
              s.psi_T_[i]=c.first.psi_T_;
              s.omega_B_[i]=c.first.omega_B_;
              s.psi_B_[i]=c.first.psi_B_;
              s.rho_[i]=c.first.rho_;
              ++i;
              if(i>=sp.numSimPoints())
                {
                  s.isValid_=c.first.isValid_;
                  return s;
                }



            }
          dt_run=dts[its].first;
          c=nextCrankNicholson(p,c.first,dt_run,0,dts[0].second,UseDerivative);
          t+=dt_run;
          ++its;

        }
      s.isValid_=false;
      return s;
    }
}









Di<CortexSimulation> SimplestModel::simulate_Adapted_i(const Di<Parameters>& par, const Di<SimplestModel::Param>& p, const Experiment &sp, double dx, double desiredError, double dt0, std::size_t nPoints_per_decade,
                                                       double dtmax,double dtmin,
                                                       double tequilibrio) const

{
  double dtprod=std::pow(10,1.0/nPoints_per_decade);


  std::pair<Di<CortexState>,double> c{p.second.size(),0};

  Di<CortexSimulation> s(p.second.size());
  c.first.first=init(p.first,sp,dx);
  for (std::size_t i=0; i<p.second.size(); ++i)
    c.first.second[i]=init(p.second[i],sp,dx);

  s.first=CortexSimulation(c.first.first,sp.numSimPoints());
  s.first.p_=par.first;
  s.first.dt_=dtmax;
  s.first.h_=sp.h();
  for (std::size_t i=0; i<p.second.size(); ++i)
    {
      s.second[i]=CortexSimulation(c.first.second[i],sp.numSimPoints());
      s.second[i].p_=par.second[i];
      s.second[i].dt_=dtmax;
      s.second[i].h_=sp.h();

    }
  double t=-tequilibrio;
  unsigned i=0;

  while ((t+dtmax<=0)&&(std::isfinite(c.second)))
    {
      c=nextEuler_Adapt_i(p,c,dtmax,dtmin,desiredError);
      t+=dtmax;
    }
  if (!std::isfinite(c.second))
    {
      s.first.isValid_=false;
      return s;
    }
  else
    {
      double dt_run=dt0;
      addDamp(c.first.first,p.first);
      for (std::size_t ipar=0; ipar<p.second.size(); ++ipar)
        addDamp(c.first.second[ipar],p.second[ipar]);

      while ((i<sp.numSimPoints())&&(c.first.first.isValid_)&&(std::isfinite(c.second)))
        {
          t+=dt0;
          c=nextEuler_Adapt_i(p,c,dt0,dtmin,desiredError);
          if (t>=sp.tSimul(i)&&(std::isfinite(c.second)))
            {
              s.first.t_[i]=t;
              s.first.maxlogErrt_[i]=c.second;
              s.first.omega_T_[i]=c.first.first.omega_T_;
              s.first.psi_T_[i]=c.first.first.psi_T_;
              s.first.omega_B_[i]=c.first.first.omega_B_;
              s.first.psi_B_[i]=c.first.first.psi_B_;
              s.first.rho_[i]=c.first.first.rho_;
              for (std::size_t ipar=0; ipar<p.second.size(); ++ipar)
                {
                  s.second[ipar].t_[i]=t;
                  s.second[ipar].maxlogErrt_[i]=c.second;
                  s.second[ipar].omega_T_[i]=c.first.second[ipar].omega_T_;
                  s.second[ipar].psi_T_[i]=c.first.second[ipar].psi_T_;
                  s.second[ipar].omega_B_[i]=c.first.second[ipar].omega_B_;
                  s.second[ipar].psi_B_[i]=c.first.second[ipar].psi_B_;
                  s.second[ipar].rho_[i]=c.first.second[ipar].rho_;

                }
              c.second=0;

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
      if (!c.first.first.isValid_||(!std::isfinite(c.second)))
        {
          s.first.isValid_=false;
          for (std::size_t ipar=0; ipar<p.second.size(); ++ipar)
            s.second[ipar].isValid_=false;
          return s;
        }
      else
        {
          s.first.isValid_=true;
          for (std::size_t ipar=0; ipar<p.second.size(); ++ipar)
            s.second[ipar].isValid_=true;
          return s;
        }
    }
}



void SimplestModel::addDamp(CortexState &c, const SimplestModel::Param &p) const
{
  unsigned i=0;
  double damp=p.prot_concentration_*p.DAMP_psi_ratio_/p.DAMP_MW_/1000;
  while (c.x_[i]<p.inj_width_+c.x_[0])
    {
      double inj_size_cell=std::min(p.inj_width_+c.x_[0],c.x_[i]+c.dx_[i])-c.x_[i];

      c.psi_T_[i]+=damp*inj_size_cell/c.dx_[i]/p.epsilon_;
      ++i;
    }
  i=0;
  if (p.DAMP_omega_ratio_>0)
    {
      double damp2=p.prot_concentration_*p.DAMP_omega_ratio_/p.DAMP_MW_/1000;
      while (c.x_[i]<p.inj_width_+c.x_[0])
        {
          double inj_size_cell=std::min(p.inj_width_+c.x_[0],c.x_[i]+c.dx_[i])-c.x_[i];

          c.omega_T_[i]+=damp2*inj_size_cell/c.dx_[i]/p.epsilon_;
          ++i;
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

std::vector<std::string> BaseModel::getModelStateLabels() const
{
  return {"NotAstrocyte","type0","typeI","typeII","typeIII","typeIV","typeV"};
}

std::vector<std::string> BaseModel::getObservedStateLabels() const
{
  return {"typeI","typeII","typeIII","typeIV","typeV"};

}

std::vector<std::string> BaseModel::getApoptoticStatesAtInjuryLabels() const
{
  return  getModelStateLabels();
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

  o[Model_01::number()]=new Model_01;
  o[Model_02::number()]=new Model_02;
  o[Model_03::number()]=new Model_03;
  o[Model_10::number()]=new Model_10;
  o[Model_101::number()]=new Model_101;
  o[Model_11::number()]=new Model_11;

  o[Model_12::number()]=new Model_12;

  o[Model_13::number()]=new Model_13;


  return o;

}
std::map<double, BaseModel *> BaseModel::models_=getModels();

