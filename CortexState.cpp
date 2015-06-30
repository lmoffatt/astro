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








CortexState SimplestModel::nextEuler(const SimplestModel::Param &p, const CortexState &c, double dt)
{
  auto dPsi=dPsi_dt(p,c);
  auto dOmega=dOmega_dt(p,c);
  auto dRhoOmegaLigand=dRhoOmegaLigand_dt(p,c);
  auto dRhoOmegaState=dRhoOmegaState_dt(p,c);
  auto dRhoPsiLigand=dRhoPsiLigand_dt(p,c);
  auto dRhoPsiState=dRhoPsiState_dt(p,c);


  CortexState out(c);
  out.omega_+=dOmega*dt;
  out.psi_+=dPsi*dt;
  out.rho_i_+=(dRhoPsiLigand+dRhoPsiState)*dt;
  out.rho_j_+=(dRhoOmegaLigand+dRhoOmegaState)*dt;

  out.rho_=sum(out.rho_i_);

  return out;


}

CortexState SimplestModel::init(const SimplestModel::Param &p,const SimplestModel::SimParam &s)
{
  CortexState o(s.x_,s.dx_,p.N_,s.n_,s.dn_,p.M_,s.m_,s.dm_);


  unsigned numX=o.x_.size();
  for (unsigned x=0; x<numX; ++x)
    {
      o.rho_[x][0]=p.nAstr_[x];
      o.rho_i_[x][0][0]=p.nAstr_[x];
      o.rho_j_[x][0][0]=p.nAstr_[x];
    }
  return o;
}

CortexState SimplestModel::injectDamp(const CortexState &c, double damp)
{
     CortexState o(c);
     o.psi_[0]+=damp;  // TODO: correct by dx
     return o;
}

std::vector<double> SimplestModel::dPsi_dt(const Param &p,
                                           const CortexState &c)
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();
  unsigned numI=c.rho_i_.front().front().size();

  
  std::vector<double> o(numX,0);
  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      if (x>0)
        Jn=2.0*p.Dpsi_*(c.psi_[x-1]-c.psi_[x])/(c.dx_[x-1]+c.dx_[x]);
      if (x+1<numX)
        Jp=2.0*p.Dpsi_*(c.psi_[x+1]-c.psi_[x])/(c.dx_[x+1]+c.dx_[x]);

      double a=0;

      for (unsigned k=0; k<numK; ++k)
        {
          double nRecep=0;
          for (unsigned i=0; i<numI; ++i)
            {
              nRecep+=(c.N_[k]-c.n_[i])*c.rho_i_[x][k][i];
            }
          a+=p.kon_psi_[k]*nRecep;
        }
      o[x]=(Jn+Jp)/c.dx_[x]-c.psi_[x]*a;
    }
  return o;
}

std::vector<double> SimplestModel::dOmega_dt(const Param &p,
                                             const CortexState &c)
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();
  unsigned numI=c.rho_j_.front().front().size();


  std::vector<double> o(numX,0);
  for (unsigned x=0; x<numX; ++x)
    {
      double Jn=0;
      double Jp=0;
      if (x>0)
        Jn=2.0*p.Domega_*(c.omega_[x-1]-c.omega_[x])/(c.dx_[x-1]+c.dx_[x]);
      if (x+1<numX)
        Jp=2.0*p.Domega_*(c.omega_[x+1]-c.omega_[x])/(c.dx_[x+1]+c.dx_[x]);

      double a=0;
      double sig=0;
      for (unsigned k=0; k<numK; ++k)
        {
          double nRecep=0;
          for (unsigned i=0; i<numI; ++i)
            {
              nRecep+=(c.N_[k]-c.n_[i])*c.rho_j_[x][k][i];
            }
          a+=p.kon_omega_[k]*nRecep;
          sig+=p.ksig_omega_[k]*c.rho_[x][k];
        }
      o[x]=(Jn+Jp)/c.dx_[x]+sig-c.omega_[x]*a;
    }
  return o;
}


std::vector<std::vector<std::vector<double>>> SimplestModel::
dRhoOmegaLigand_dt(const Param &p,
                   const CortexState &c)
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();
  unsigned numJ=c.rho_j_.front().front().size();


  std::vector<std::vector<std::vector<double>>> dr(numX,std::vector<std::vector<double>>(numK,std::vector<double>(numJ,0)));
  for (unsigned x=0; x<numX; ++x)
    {
      for (unsigned k=0; k<numK; ++k)
        {
          for (unsigned j=0; j<numJ; ++j)
            {
              double Jpos,Jneg;
              if (j+1<numJ)
                Jpos=p.kcat_omega_[k]*c.m_[j+1]/c.dm_[j+1]*c.rho_j_[x][k][j+1]
                    -p.kon_omega_[k]*c.omega_[x]*(c.M_[k]-c.m_[j])/c.dm_[j]
                    *c.rho_j_[x][k][j];
              else
                Jpos=0;
              if (j>0)
                Jneg=-p.kcat_omega_[k]*c.m_[j]/c.dm_[j]*c.rho_j_[x][k][j]
                    -p.kon_omega_[k]*c.omega_[x]*(c.M_[k]-c.m_[j-1])/c.dm_[j-1]
                    *c.rho_j_[x][k][j-1];
              else
                Jneg=0;
              dr[x][k][j]=Jpos+Jneg;

            }
        }
    }
  return dr;
}

std::vector<std::vector<std::vector<double>>> SimplestModel::
dRhoOmegaState_dt(const Param &p,
                  const CortexState &c)
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();
  unsigned numJ=c.rho_j_.front().front().size();
  unsigned numI=c.rho_i_.front().front().size();


  std::vector<std::vector<std::vector<double>>> dr(
        numX,std::vector<std::vector<double>>(numK,std::vector<double>(numJ,0)));
  for (unsigned x=0; x<numX; ++x)
    {
      for (unsigned k=0; k<numK; ++k)
        {
          double gkk1_n=0;
          double gk_1k_n=0;
          double ak_n=0;

          for (unsigned i=0; i<numI; ++i)
            {
              gkk1_n+=g_psi(p,k,c.n_[i])*c.rho_i_[x][k][i];
              ak_n+=a_psi(p,k,c.n_[i])*c.rho_i_[x][k][i];
            }
          gkk1_n/=c.rho_[x][k];
          ak_n/=c.rho_[x][k];

          if (k>0)
            {
              for (unsigned i=0; i<numI; ++i)
                {
                  gk_1k_n+=g_psi(p,k-1,c.n_[i])*c.rho_i_[x][k-1][i];
                }
              gk_1k_n/=c.rho_[x][k-1];
            }



          for (unsigned j=0; j<numJ; ++j)
            {
              double Jpos,Jneg;

              if (k+1<numK)
                Jpos=p.g_rev_[k+1]*c.rho_j_[x][k+1][j]
                    -gkk1_n*c.rho_j_[x][k][j]
                    -g_omega(p,k,c.m_[j])*c.rho_j_[x][k][j];
              else
                Jpos=0;
              if (k>0)
                Jneg=-p.g_rev_[k]*c.rho_j_[x][k][j]
                    +gk_1k_n*c.rho_j_[x][k-1][j]
                    +g_omega(p,k-1,c.m_[j])*c.rho_j_[x][k-1][j];
              else Jneg=0;
              double JA=-ak_n*c.rho_j_[x][k][j]
                  -a_omega(p,k,c.m_[j])*c.rho_j_[x][k][j];

              dr[x][k][j]=Jpos+Jneg+JA;

            }
        }
    }
  return dr;
}



std::vector<std::vector<std::vector<double>>> SimplestModel::
dRhoPsiLigand_dt(const Param &p,
                 const CortexState &c)
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();
  unsigned numI=c.rho_i_.front().front().size();


  std::vector<std::vector<std::vector<double>>>
      dr(numX,std::vector<std::vector<double>>(numK,std::vector<double>(numI)));
  for (unsigned x=0; x<numX; ++x)
    {
      for (unsigned k=0; k<numK; ++k)
        {
          for (unsigned i=0; i<numI; ++i)
            {
              double Jpos,Jneg;
              if (i+1<numI)
                Jpos=p.kcat_psi_[k]*c.n_[i+1]/c.dn_[i+1]*c.rho_i_[x][k][i+1]
                    -p.kon_psi_[k]*c.psi_[x]*(c.N_[k]-c.n_[i])/c.dn_[i]
                    *c.rho_i_[x][k][i];
              else Jpos=0;
              if (i>0)
                Jneg=-p.kcat_psi_[k]*c.n_[i]/c.dn_[i]*c.rho_i_[x][k][i]
                    -p.kon_psi_[k]*c.psi_[x]*(c.N_[k]-c.n_[i-1])/c.dn_[i-1]
                    *c.rho_i_[x][k][i-1];
              else
                Jneg=0;
              dr[x][k][i]=Jpos+Jneg;

            }
        }
    }
  return dr;
}

std::vector<std::vector<std::vector<double>>> SimplestModel::
dRhoPsiState_dt(const Param &p,
                const CortexState &c)
{
  unsigned numX=c.rho_.size();
  unsigned numK=c.rho_.front().size();
  unsigned numJ=c.rho_j_.front().front().size();
  unsigned numI=c.rho_i_.front().front().size();


  std::vector<std::vector<std::vector<double>>> dr(
        numX,std::vector<std::vector<double>>(numK,std::vector<double>(numI,0)));
  for (unsigned x=0; x<numX; ++x)
    {
      for (unsigned k=0; k<numK; ++k)
        {
          double gkk1_n=0;
          double gk_1k_n=0;
          double ak_n=0;
          if (k>0)
            {
              for (unsigned j=0; j<numJ; ++j)
                {
                  gk_1k_n+=g_omega(p,k-1,c.n_[j])*c.rho_j_[x][k-1][j];
                }
              gk_1k_n/=c.rho_[x][k-1];
            }
          for (unsigned j=0; j<numJ; ++j)
            {
              gkk1_n+=g_omega(p,k,c.n_[j])*c.rho_j_[x][k][j];
              ak_n+=a_omega(p,k,c.n_[j])*c.rho_j_[x][k][j];
            }
          gkk1_n/=c.rho_[x][k];
          ak_n/=c.rho_[x][k];


          for (unsigned i=0; i<numI; ++i)
            {
              double Jpos,Jneg;
              if (k+1<numK)
                Jpos=p.g_rev_[k+1]*c.rho_i_[x][k+1][i]
                    -gkk1_n*c.rho_i_[x][k][i]
                    -g_omega(p,k,c.m_[i])*c.rho_i_[x][k][i];
              else Jpos=0;
              if (k>0)
                Jneg=-p.g_rev_[k]*c.rho_i_[x][k][i]
                    +gk_1k_n*c.rho_i_[x][k-1][i]
                    +g_omega(p,k-1,c.m_[i])*c.rho_i_[x][k-1][i];
              else Jneg=0;
              double JA=-ak_n*c.rho_i_[x][k][i]
                  -a_omega(p,k,c.m_[i])*c.rho_i_[x][k][i];

              dr[x][k][i]=Jpos+Jneg+JA;

            }
        }
    }
  return dr;
}









double SimplestModel::g_omega(const Param &p,unsigned k,double m)
{
  return p.g_max_omega_[k]*std::pow(m,p.h_omega_[k])
      /(std::pow(m,p.h_omega_[k])+std::pow(p.n_50_omega_[k],p.h_omega_[k]));
}

double SimplestModel::a_omega(const Param &p,unsigned k,double m)
{
  return p.a_max_omega_[k]*std::pow(m,p.ha_omega_[k])
      /(std::pow(m,p.ha_omega_[k])+std::pow(p.na_50_omega_[k],p.ha_omega_[k]));
}


double SimplestModel::g_psi(const Param &p,unsigned k,double m)
{
  return p.g_max_psi_[k]*std::pow(m,p.h_psi_[k])
      /(std::pow(m,p.h_psi_[k])+std::pow(p.n_50_psi_[k],p.h_psi_[k]));
}

double SimplestModel::a_psi(const Param &p,unsigned k,double m)
{
  return p.a_max_psi_[k]*std::pow(m,p.ha_psi_[k])
      /(std::pow(m,p.ha_psi_[k])+std::pow(p.na_50_psi_[k],p.ha_psi_[k]));
}

void SimplestModel::run(const SimplestModel::SimParam &s, const SimplestModel::Param &p)
{
  CortexState c=init(p,s);

  double t=0;
  double dt=s.dt_;
  while (t<s.teq_)
    {
      c=nextEuler(p,c,dt);
      t+=dt;
    }

  c.addDamp(s.);
  while (t<s.tsim_)
    {
      c=nextEuler(p,c,dt);
      t+=dt;
    }


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
