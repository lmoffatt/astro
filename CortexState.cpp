#include "CortexState.h"
#include"read.h"
#include <string>
#include <sstream>
#include <list>
#include <iostream>
#include <set>
#include <algorithm>
void TissuePhoto::read(std::string& line, std::istream &s)
{
  std::string name;
  std::stringstream ss(line);
  ss>>name;
  if (name=="foto")
    {
      ss>>num;
      std::stringstream sid(id);
      sid<<"f"<<num;
      id=sid.str();
      safeGetline(s,line);
      while (true)
        {
          ss.clear();
          ss.str(line);
          name.clear();
          ss.clear();
          ss>>name;
          while (name.empty()&&safeGetline(s,line))
            {

              ss.str(line);
              ss.clear();
              ss>>name;
            }
          if (line.find("celulas")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              double x,y,prob;
              unsigned typo;
              while (ss2>>x>>y>>typo>>prob)
                {
                  std::string idAstro;
                  std::stringstream sid(idAstro);
                  sid<<this->id<<"a"<<astr_.size();

                  astr_.push_back(Astrocyte(sid.str(),x,y,typo,prob));
                  safeGetline(s,line);
                  ss2.str(line);
                  ss2.clear();
                }
            }
          else if (line.find("limite foto")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2(line);
              double x,y;
              std::vector<position> v;
              while (ss2>>x>>y)
                {
                  v.push_back(position(x,y));
                  safeGetline(s,line);
                  ss2.str(line);
                  ss2.clear();
                }
              lf_=LimiteFoto(v);

            }
          else if (line.find("limite tejido")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2(line);
              double x,y;
              std::vector<position> v;
              while (ss2>>x>>y)
                {
                  v.push_back(position(x,y));
                  safeGetline(s,line);
                  ss2.str(line);
                  ss2.clear();
                }
              lt_=LimiteTejido(v);

            }
          else if (line.find("limite lesion")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2(line);
              double x,y;
              std::vector<position> v;
              while (ss2>>x>>y)
                {
                  v.push_back(position(x,y));
                  safeGetline(s,line);
                  ss2.str(line);
                  ss2.clear();
                }
              ll_=LimiteLesion(v);

            }

          else if (line.find("vaso")!=line.npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2(line);
              double x,y;
              std::vector<position> v;
              while (ss2>>x>>y)
                {
                  v.push_back(position(x,y));
                  safeGetline(s,line);
                  ss2.str(line);
                  ss2.clear();
                }
              vasos_.push_back(LimiteVaso(v));

            }

          else if (line.find("pin")!=line.npos)
            {
              std::stringstream ss(line);
              std::string pinName; unsigned pinNum;
              ss>>pinName>>pinNum;

              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2(line);
              double x,y;
              std::vector<position> v;
              while (ss2>>x>>y)
                {
                  v.push_back(position(x,y));
                  safeGetline(s,line);

                  ss2.str(line);
                  ss2.clear();
                }
              pines_[pinNum]=Pin(v);

            }

          else break;

        }
    }

}

void TissuePhoto::write(std::ostream &s)
{
  s<<"foto"<<"\t"<<num<<"\n";

  for (auto e:pines_)
    {
      s<<"pin"<<"\t"<<e.first<<"\n";
      s<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      Pin pin=e.second;
      for (auto e:pin.limits())
        s<<e.x<<"\t"<<e.y<<"\n";
      s<<"\n";
    }

  s<<"celulas"<<"\n";
  s<<astr_.front().getHeader()<<"\n";
  for (Astrocyte a:astr_)
    {
      a.write(s);
    }
  s<<"\n";

  for (LimiteVaso v:vasos_)
    {
      s<<"vaso"<<"\n";
      s<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:v.limits())
        s<<e.x<<"\t"<<e.y<<"\n";
      s<<"\n";
    }


  if (!ll_.limits().empty())
    {
      s<<"limite lesion"<<"\n";
      s<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:ll_.limits())
        s<<e.x<<"\t"<<e.y<<"\n";
      s<<"\n";
    }

  if (!lt_.limits().empty())
    {
      s<<"limite tejido"<<"\n";
      s<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:lt_.limits())
        s<<e.x<<"\t"<<e.y<<"\n";
      s<<"\n";
    }
  if (!lf_.limits().empty())
    {
      s<<"limite foto"<<"\n";
      s<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:lf_.limits())
        s<<e.x<<"\t"<<e.y<<"\n";
      s<<"\n";
    }

}

bool TissuePhoto::align( TissuePhoto &other)
{
  // find common pins

  std::vector<unsigned> commonPins;

  for (auto mypin:pines_)
    {
      if (other.pines_.find(mypin.first)!=other.pines_.end())
        commonPins.push_back(mypin.first);

    }

  if (commonPins.empty())
    return false;
  else
    {// calculate the displacement for all found pins
      std::vector<double> dxs;
      std::vector<double> dys;
      std::cerr<<"\nfoto"<<other.num<<"\n";
      double sdx=0,sdy=0;
      for (auto aPin:commonPins)
        {
          double dx,dy;
          Pin mine, his;
          mine=pines_[aPin];
          his=other.pines_[aPin];
          dx=mine.limits().front().x-his.limits().front().x;
          std::cerr<<"pin"<<aPin<<"\t";
          std::cerr<<"dx="<<dx<<"\t";
          dxs.push_back(dx);
          dy=mine.limits().front().y-his.limits().front().y;
          std::cerr<<"dy="<<dy<<"\n";


          dys.push_back(dy);
          sdx+=dx;
          sdy+=dy;
        }

      // now calculate the mean sin of the rotation angle for all found pins
      // together with the expansion factor
      // and the measurement error
      double dx=sdx/dxs.size();
      double dy=sdy/dxs.size();

      // now correct the position of the other photo taking this as good

      other.correctPosition(dx,dy);

      // finally add the corrected pins that do not appear in this to the list of this

      for (auto hisPin:other.pines_)
        {
          if (pines_.find(hisPin.first)==pines_.end())
            pines_.insert(hisPin);
        }
      return true;
    }
}

void TissuePhoto::correctPosition(double dx, double dy)
{
  this->ll_.correctPosition(dx,dy);
  this->lt_.correctPosition(dx,dy);
  this->lf_.correctPosition(dx,dy);
  for (auto& pin:pines_)
    {
      pin.second.correctPosition(dx,dy);
    }
  for (auto& v:vasos_)
    {
      v.correctPosition(dx,dy);
    }

  for (auto& a:astr_)
    {
      a.correctPosition(dx,dy);
    }
}








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
    out.omega_+=dOmega*dt;
  out.psi_+=dPsi*dt;

  out.rho_+=dRho*dt;

  return out;


}

CortexState SimplestModel::init(const SimplestModel::Param &p, const CortexExperiment &s)const
{
  CortexState o(s.x_,s.dx_,s.h_,p.N_.size());


  unsigned numX=o.x_.size();
  for (unsigned x=0; x<numX; ++x)
    {
      o.rho_[x][0]=p.dens_Neur_*std::pow(s.h_,2)*s.dx_[x]*1000;  // dens is in cells per liter
      o.rho_[x][1]=p.dens_Astr_*std::pow(s.h_,2)*s.dx_[x]*1000;

    }
  return o;
}

CortexState SimplestModel::injectDamp(const CortexState &c, double damp)const
{
  CortexState o(c);
  o.psi_[0]+=damp;  // TODO: correct by dx
  return o;
}

std::vector<double> SimplestModel::dPsi_dt(const Param &p,
                                           const CortexState &c)const
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();

  //Avogadros number
  const double NAv=6.022E23;

  std::vector<double> o(numX,0);

  double f=p.kcat_psi_/(NAv*p.epsilon_*c.h_*c.h_)*p.kon_psi_/1000.0;

  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      if (x>0)
        Jn=2.0*p.Dpsi_*(c.psi_[x-1]-c.psi_[x])/(c.dx_[x-1]+c.dx_[x]);
      if (x+1<numX)
        Jp=2.0*p.Dpsi_*(c.psi_[x+1]-c.psi_[x])/(c.dx_[x+1]+c.dx_[x]);

      double d=p.kcat_psi_+p.kon_psi_*c.psi_[x];
      double Nt=0;
      for (unsigned k=0; k<numK; ++k)
        {
          Nt+=p.N_[k]*c.rho_[x][k];
        }

      double a=f/c.dx_[x]/d*Nt;
 //    if ((c.psi_[x]>0)&&a>1)
 //     std::cerr<<"psi ="<<c.psi_[x]<<"a="<<a<<"  dx="<<c.dx_[x]<<"\t";


      o[x]=(Jn+Jp)/c.dx_[x]-a*c.psi_[x];
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
  double f=p.kcat_omega_/(NAv*p.epsilon_*c.h_*c.h_)*p.kon_omega_/1000;

  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      if (x>0)
        Jn=2.0*p.Domega_*(c.omega_[x-1]-c.omega_[x])/(c.dx_[x-1]+c.dx_[x]);
      if (x+1<numX)
        Jp=2.0*p.Domega_*(c.omega_[x+1]-c.omega_[x])/(c.dx_[x+1]+c.dx_[x]);

      double d=p.kcat_omega_+p.kon_omega_*c.omega_[x];
      double sig=0;
      double Mt=0;
      for (unsigned k=0; k<numK; ++k)
        {
          Mt+=p.M_[k]*c.rho_[x][k];
          sig+=p.ksig_omega_[k]*c.rho_[x][k];
        }
      double a=f/c.dx_[x]*c.omega_[x]/d*Mt;

      if (c.omega_[x]>0)
              a=a+0;
       sig/=(NAv*p.epsilon_*c.dx_[x]*c.h_*c.h_)/1000;
      o[x]=(Jn+Jp)/c.dx_[x]+sig-a;
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
          for (unsigned k=0; k<numK; ++k)
            {
              double Jr,Jl,Ja;
              if (k+1<numK)
                Jr=p.g_left_[k+1]*c.rho_[x][k+1]
                    -(p.g_rigth_[k]
                      +p.g_max_psi_[k]*p.kon_psi_*c.psi_[x]
                      /(p.kcat_psi_+p.kon_psi_*c.psi_[x])
                      +p.g_max_omega_[k]*p.kon_omega_*c.omega_[x]
                      /(p.kcat_omega_+p.kon_omega_*c.omega_[x])
                      )*c.rho_[x][k];
              else
                Jr=0;
              if (k>0)
                Jl=-p.g_left_[k]*c.rho_[x][k]
                    +(p.g_rigth_[k-1]
                    +p.g_max_psi_[k-1]*p.kon_psi_*c.psi_[x]
                    /(p.kcat_psi_+p.kon_psi_*c.psi_[x])
                    +p.g_max_omega_[k-1]*p.kon_omega_*c.omega_[x]
                    /(p.kcat_omega_+p.kon_omega_*c.omega_[x])
                    )*c.rho_[x][k-1];
              else
                Jl=0;
              Ja= (+p.a_[k]
                   +p.a_psi_[k]*p.kon_psi_*c.psi_[x]
                   /(p.kcat_psi_+p.kon_psi_*c.psi_[x])
                   +p.a_omega_[k]*p.kon_omega_*c.omega_[x]
                   /(p.kcat_omega_+p.kon_omega_*c.omega_[x])
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
              double Jr,Jl,Ja;
              if (k+1<numK)
                Jr=p.g_left_[k+1]*c.rho_[x][k+1]
                    -(p.g_rigth_[k]
                      +p.g_max_psi_[k]*p.kon_psi_*c.psi_[x]
                      /(p.kcat_psi_+p.kon_psi_*c.psi_[x])
                      )*c.rho_[x][k];
              else
                Jr=0;
              if (k>0)
                Jl=-p.g_left_[k]*c.rho_[x][k]
                    +(p.g_rigth_[k-1]
                    +p.g_max_psi_[k-1]*p.kon_psi_*c.psi_[x]
                    /(p.kcat_psi_+p.kon_psi_*c.psi_[x])
                    )*c.rho_[x][k-1];
              else
                Jl=0;
              Ja= (p.a_[k]
                   +p.a_psi_[k]*p.kon_psi_*c.psi_[x]
                   /(p.kcat_psi_+p.kon_psi_*c.psi_[x])
                   )*c.rho_[x][k];

              dr[x][k]=Jl+Jr-Ja;

            }
        }
    }

  return dr;
}





CortexSimulation SimplestModel::simulate(const SimplestModel::Param &p,
                                         const CortexExperiment &sp
                                         ,double dt)const
{


  std::cout<<"starts a Simulation \n";
  std::cout<<"dt ="<<dt<<"\n";


  CortexState c=init(p,sp);

  unsigned numSamples=sp.tsim_/sp.sample_time_;


  CortexSimulation s(c,numSamples);

  double t=0;
  unsigned i=0;
  s.t_[i]=t;
  s.sdt_[i]=dt;
  s.omega_[i]=c.omega_;
  s.psi_[i]=c.psi_;
  s.rho_[i]=c.rho_;

  while (t<sp.teq_)
    {
      c=nextEuler(p,c,dt);
      t+=dt;
      if (t>sp.sample_time_*(i+1))
        {
          ++i;
          s.t_[i]=t;
          s.sdt_[i]=dt;
          s.omega_[i]=c.omega_;
          s.psi_[i]=c.psi_;
          s.rho_[i]=c.rho_;

        }
    }

  c.addDamp(p.damp_);
  while (t<sp.tsim_)
    {
      c=nextEuler(p,c,dt);
      t+=dt;
      if (t>sp.sample_time_*(i+1))
        {
          ++i;
          s.t_[i]=t;
          s.sdt_[i]=dt;
          s.omega_[i]=c.omega_;
          s.psi_[i]=c.psi_;
          s.rho_[i]=c.rho_;
          if (i % 50 ==0)
            std::cerr<<"\t sample \t"<<i;
        }
    }
  return s;


}





void TissueSection::align()
{
  // use first photo as the master

  std::list<unsigned> remaining;
  for (auto it:fotos)
    remaining.push_back(it.first);
  unsigned master=remaining.front();
  remaining.pop_front();
  TissuePhoto& f=fotos[master];

  while (!remaining.empty())
    {
      for (auto it=remaining.begin(); it!=remaining.end(); ++it)
        {
          if (f.align(fotos[*it]))
            {
              it=remaining.erase(it);
            }
        }
    }

}



void TissueSection::merge()
{
  TissuePhoto foto;

  for (auto it:fotos)
    {
      foto.include(it.second);
    }


  foto.id="Merged";
  foto.num=0;
  fotos.clear();
  fotos[0]=foto;



}

void TissueSection::distances()
{
  for (auto& f:fotos)
    f.second.calculate_distances();
}



class Intervaling
{
public:
  virtual unsigned num() const=0;

  virtual std::vector<double> limits()const=0;

  unsigned npos=std::numeric_limits<unsigned>::max();
  virtual unsigned getIndex(double x)const=0;
  ~Intervaling(){}
};


class pointDefined: public Intervaling
{


  // Intervaling interface
public:
  virtual unsigned num() const
  {
    return ticks_.size();
  }
  virtual std::vector<double> limits() const
  {
    std::vector<double> out(ticks_.size());
    unsigned i=0;
    for (auto it=ticks_.begin();it!=ticks_.end(); ++it)
      {
        out[i]=it->first;
        ++i;
      }
    return out;

  }
  virtual unsigned getIndex(double x) const
  {
    auto it=ticks_.lower_bound(x);
    if (it!=ticks_.end())
      return it->second;
    else
      return npos;
  }

  pointDefined(std::vector<double> c)
  {
    std::sort(c.begin(),c.end());
    for (unsigned i=0; i<c.size(); ++i)
      ticks_[c[i]]=i;




  }

private:
  std::map<double, unsigned> ticks_;

};

CortexMeasure *TissuePhoto::measure(std::string id,std::vector<double> x)
{
  pointDefined p(x);

  std::vector<std::vector<double>> numx(p.num(),std::vector<double>(Astrocyte::numTypes(),0));


  std::vector<std::vector<double>> var(p.num(),std::vector<double>(Astrocyte::numTypes(),0));


  for (Astrocyte a:astr_)
    {
      switch (a.type()){
        case 1:
        case 2:
          numx[p.getIndex(a.distance_to_lession())][a.type()-1]++;
          break;
        case 3:
          numx[p.getIndex(a.distance_to_lession())][2]+=a.prob();
          numx[p.getIndex(a.distance_to_lession())][0]+=1.0-a.prob();
          // we only sum the variance for the landing state;
          //we calculate the whole covariance matrix at the end
          var[p.getIndex(a.distance_to_lession())][2]+=a.prob()*(1.0-a.prob());
          break;
        case 4:
        case 5:
        case 6:
        case 7:
          numx[p.getIndex(a.distance_to_lession())][a.type()-1]+=a.prob();
          numx[p.getIndex(a.distance_to_lession())][a.type()-2]+=1.0-a.prob();
          // we only sum the variance for the landing state;
          //we calculate the whole covariance matrix at the end
          var[p.getIndex(a.distance_to_lession())][a.type()-1]+=a.prob()*(1.0-a.prob());
          break;
        default:
          break;

        }
    }
  std::vector<std::vector<std::vector<double>>>
      covar(p.num()
            ,std::vector<std::vector<double>>(Astrocyte::numTypes(),
                                              std::vector<double>(Astrocyte::numTypes(),0)));
  for (unsigned ix=0; ix<p.num(); ++ix)
    for (unsigned i=0; i<Astrocyte::numTypes(); ++i)
      {
        switch (i){
          case 0:
          case 1:
            break;
          case 2:
            covar[ix][2][2]+=var[ix][2];
            covar[ix][0][0]+=var[ix][2];
            covar[ix][0][2]-=var[ix][2];
            covar[ix][2][0]-=var[ix][2];
            break;
          case 3:
          case 4:
          case 5:
          case 6:
            covar[ix][i][i]+=var[ix][i];
            covar[ix][i-1][i-1]+=var[ix][i];
            covar[ix][i][i-1]-=var[ix][i];
            covar[ix][i-1][i]-=var[ix][i];
            break;
          default:
            break;
          }


      }




  CortexMeasure* m=new CortexMeasure(id,p.limits(),numx,covar);
  return m;


}









void tissueElement::correctPosition(double dx, double dy)
{
  for (auto& pos:limits_)
    {
      pos.x+=dx;
      pos.y+=dy;
    }
}


std::ostream &Astrocyte::write(std::ostream &s)
{
  s<<id()<<"\t";
  if (distance_to_lession_<std::numeric_limits<double>::infinity())
    s<<distance_to_lession_<<"\t"<<distance_to_tissue_<<"\t";
  s<<pos().x<<"\t"<<pos().y<<"\t"<<type()<<"\t"<<prob()<<"\n";
  return s;
}

std::string Astrocyte::getHeader()
{
  std::string ss;
  std::stringstream s(ss);
  s<<"ID"<<"\t";
  if (distance_to_lession_<std::numeric_limits<double>::infinity())
    s<<"d_les (um)"<<"\t"<<"d_tej (um)"<<"\t";
  s<<"X (um)"<<"\t"<<"Y (um)"<<"\t"<<"TIPO"<<"\t"<<"PB"<<"\n";

  return s.str();
}

void Astrocyte::calculateDistances(const TissuePhoto &f)
{
  distance_to_lession_=f.ll_.distance(*this);
  distance_to_tissue_=f.lt_.distance(*this);
}





void CortexExperiment::read(std::string &line, std::istream &s)
{

  std::string name;
  std::stringstream ss(line);
  ss>>name;
  if (name=="experiment")
    {
      safeGetline(s,line);
      while (true)
        {
          ss.str(line);
          name.clear();
          ss.clear();
          ss>>name;
          while (name.empty()&&safeGetline(s,line))
            {
              ss.str(line);
              ss.clear();
              ss>>name;
            }
          if (line.find("nodes")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              double x0,step,end,step2=0,end2;
              char colon, op0,op1,colon2,colon3;
              ss2>>x0>>colon>>op0>>step>>colon2>>end;
              x0=x0*1e-6;
              step=step*1e-6;
              end=end*1e-6;
              double x=x0;
              double dx0=step;
              double dxp=dx0;
              double dx;
              x_.clear();
              dx_.clear();
              while (x<end)
                {
                  x_.push_back(x);
                  dx_.push_back(step);
                  x+=step;
                }
              x_.push_back(x);
              dx_.push_back(step);

              while (ss2>>colon3)
                {
                  ss2>>op1>>step2>>colon3>>end2;
                  while (x<end2*1e-6)
                    {
                      if (op1=='*')
                        {

                          dx0*=step2;
                          x+=dx0;
                          dx=(dx0+dxp)/2;
                        }
                      else
                        {

                          dx=step2*1e-6;
                          x+=dx;
                        }

                      x_.push_back(x);
                      dx_.push_back(dx);
                    }


                }


            }
          else if (line.find("sample time")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              ss2>>this->sample_time_;
            }
          else if (line.find("thickness")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              ss2>>this->h_;
              this->h_*=1e-6;
            }
          else if (line.find("equilibrium time")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              ss2>>this->teq_;
            }
          else if (line.find("total time")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              ss2>>this->tsim_;
            }
          else break;
          line.clear();
        }
    }
}




void Parameters::read(std::string &line, std::istream &s)
{

  std::string name;
  std::stringstream ss(line);
  ss>>name;
  if (name=="parameters")
    {
      safeGetline(s,line);
      while (true)
        {
          ss.str(line);
          name.clear();
          ss.clear();
          ss>>name;
          while (name.empty()&&safeGetline(s,line))
            {
              ss.str(line);
              ss.clear();
              ss>>name;
            }
          if (name.empty())
            break;
          else
            {
              double val;
              ss>>val;
              this->push_back(name,val,ss.str());
              line.clear();
            }
        }
    }
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


std::ostream &CortexSimulation::write(std::ostream &s)
{
  s<<"simulation\n";
  s<<id_<<"\n";


  s<<"numSamples"<<"\t"<<t_.size()<<"\t";
  s<<"numNodes"<<"\t"<<x_.size()<<"\t";
  s<<"numStates"<<"\t"<<rho_.front().front().size()<<"\n"<<"\n";

  s<<"x"<<"\n";
  for (unsigned i=0; i<x_.size(); ++i)
    {
      s<<x_[i]<<"\t";
    }
  s<<"\n";

  s<<"dx"<<"\n";
  for (unsigned i=0; i<dx_.size(); ++i)
    {
      s<<dx_[i]<<"\t";
    }
  s<<"\n";



  s<<"psi"<<"\n";
  for (unsigned i=0; i<t_.size(); ++i)
    {
      s<<t_[i]<<"\t";
      for (unsigned j=0; j<psi_[i].size(); ++j)
        {
          s<<psi_[i][j]<<"\t";
        }
      s<<"\n";
    }
  s<<"\n";

  s<<"omega"<<"\n";
  for (unsigned i=0; i<t_.size(); ++i)
    {
      s<<t_[i]<<"\t";
      for (unsigned j=0; j<omega_[i].size(); ++j)
        {
          s<<omega_[i][j]<<"\t";
        }
      s<<"\n";
    }
  s<<"\n";

  for (unsigned k=0;k<rho_.front().front().size(); ++k)
    {
      s<<"rho"<<"\t"<<k<<"\n";
      for (unsigned i=0; i<t_.size(); ++i)
        {
          s<<t_[i]<<"\t";
          for (unsigned j=0; j<rho_[i].size(); ++j)
            {
              s<<rho_[i][j][k]<<"\t";
            }
          s<<"\n";
        }
      s<<"\n";
    }




return s;
}

std::ostream &CortexSimulation::write(std::ostream &s, const std::string &var, const std::string &par)
{
   std::vector<double> xs;
  std::vector<double> ts;
  std::vector<double> ks;

  unsigned numK=rho_.front().front().size();


  if (par.empty())
    {
      xs=x_;
      ts=t_;

    }
  else
    {

      std::stringstream ss(par);
      double number;
      std::string var2;
      char equal;
      while (ss>>var2>>equal>>number)
        {
          if (var2=="t")
            {
              ts.push_back(number);
            }
          else if (var2=="x")
            {
              xs.push_back(number*1e-6);
            }
          else if (var2=="k")
            {
              ks.push_back(number);
            }

        }


    }

  if (var=="psi")
    {
      if (!xs.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<xs.size(); ++j)
            {
              double x=xs[j];
              while ((x_[jr]<x)&&jr<x_.size())
                ++jr;
              js.push_back(jr);
            }

          s<<"psi vs time at different positions"<<"\n";

          s<<"time"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<x_[js[jjs]]<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<t_.size(); ++i)
            {
              s<<t_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<psi_[i][js[jjs]]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
      if (!ts.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<ts.size(); ++j)
            {
              double t=ts[j];
              while ((t_[jr]<t)&&jr<t_.size())
                ++jr;
              js.push_back(jr);
            }
          s<<"psi vs position at different times"<<"\n";

          s<<"position"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<t_[js[jjs]]<<"\t";
            }
          s<<"\n";

          for (unsigned i=0; i<x_.size(); ++i)
            {
              s<<x_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<psi_[js[jjs]][i]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }

    }
  else if (var=="omega")
    {
      if (!xs.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<xs.size(); ++j)
            {
              double x=xs[j];
              while ((x_[jr]<x)&&jr<x_.size())
                ++jr;
              js.push_back(jr);
            }

          s<<"omega vs time at different positions"<<"\n";

          s<<"time"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<x_[js[jjs]]<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<t_.size(); ++i)
            {
              s<<t_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<omega_[i][js[jjs]]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
      if (!ts.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<ts.size(); ++j)
            {
              double t=ts[j];
              while ((t_[jr]<t)&&jr<t_.size())
                ++jr;
              js.push_back(jr);
            }

          s<<"omega vs position at different times"<<"\n";

          s<<"position"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<t_[js[jjs]]<<"\t";
            }
          s<<"\n";

          for (unsigned i=0; i<x_.size(); ++i)
            {
              s<<x_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<omega_[js[jjs]][i]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }

    }
  else if (var=="rho")
    {
      if (!xs.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<xs.size(); ++j)
            {
              double x=xs[j];
              while ((x_[jr]<x)&&jr<x_.size())
                ++jr;
              js.push_back(jr);
            }

          if (ks.empty())
            for (unsigned k=0;k<numK; ++k)
              ks.push_back(k);


          s<<"rho vs time at different states and positions"<<"\n";

          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {
              s<<"time"<<"\t";
              for (unsigned ik=0; ik<ks.size(); ++ik)
                {
                  unsigned k=ks[ik];

                  s<<"k="<<k<<",x="<<x_[js[jjs]]<<"\t";
                }
              s<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<t_.size(); ++i)
            {
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {
                  s<<t_[i]<<"\t";
                  for (unsigned ik=0; ik<ks.size(); ++ik)
                    {
                      unsigned k=ks[ik];
                      s<<rho_[i][js[jjs]][k]<<"\t";
                    }
                  s<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
      if (!ts.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<ts.size(); ++j)
            {
              double t=ts[j];
              while ((t_[jr]<t)&&jr<t_.size())
                ++jr;
              js.push_back(jr);
            }

          if (ks.empty())
            for (unsigned k=0;k<numK; ++k)
              ks.push_back(k);


          s<<"rho vs position at different states and times"<<"\n";

          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {
              s<<"position"<<"\t";
              for (unsigned ik=0; ik<ks.size(); ++ik)
                {
                  unsigned k=ks[ik];

                  s<<"k="<<k<<",t="<<t_[js[jjs]]<<"\t";
                }
              s<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<x_.size(); ++i)
            {
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {
                  s<<x_[i]<<"\t";
                  for (unsigned ik=0; ik<ks.size(); ++ik)
                    {
                      unsigned k=ks[ik];
                      s<<rho_[js[jjs]][i][k]<<"\t";
                    }
                  s<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
    }

  return s;
}



void CortexSimulation::read(std::string& line, std::istream &s)
{
  std::string name;
  std::stringstream ss(line);
  ss>>name;
  if (name=="simulation")
    {
      name.clear();
      while (name.empty()&&safeGetline(s,line))
        {
          ss.str(line);
          ss.clear();
          ss>>name;
        }
      id_=name;
      name.clear();
      unsigned numSamples;
      unsigned numNodes;
      unsigned numStates;

      while (true)
        {
          while (name.empty()&&safeGetline(s,line))
            {
              ss.str(line);
              ss.clear();
              ss>>name;
            }
          if (name.find("numSamples")!=name.npos)
            {
              if  (ss>>numSamples>>name>>numNodes>>name>>numStates)
                {
                  *this=CortexSimulation(id_,numSamples,numNodes,numStates);
                }
              name.clear();

            }
          else if (name=="x")
            {
              name.clear();
              safeGetline(s,line);
              ss.clear();
              ss.str(line);
              unsigned i=0;
              double myx;
              while (ss>>myx)
                {
                  x_[i]=myx;
                  ++i;
                }

            }
          else if (name.find("dx")!=name.npos)
            {
              name.clear();
              safeGetline(s,line);
              ss.clear();
              ss.str(line);
              unsigned i=0;
              while (ss>>dx_[i])
                {
                  ++i;
                }
            }
          else if (name.find("psi")!=name.npos)
            {
              name.clear();
              unsigned i=0;
              while (true)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while (ss>>psi_[i][j])
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else if (name.find("omega")!=name.npos)
            {
              name.clear();
              unsigned i=0;
              while (true)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while (ss>>omega_[i][j])
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else if (name.find("rho")!=name.npos)
            {
              name.clear();
              unsigned k;
              ss>>k;
              unsigned i=0;

              while (i<numSamples)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while ((ss>>rho_[i][j][k])&&j<numNodes)
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else break;
        }



    }
}
