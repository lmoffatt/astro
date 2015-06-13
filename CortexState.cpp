#include "CortexState.h"
#include"read.h"
#include <string>
#include <sstream>
#include <list>
#include <iostream>

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
              line;
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
  s<<"X (um)"<<"\t"<<"Y (um)"<<"TIPO"<<"\t"<<"PB"<<"\n";
  for (Astrocyte a:astr_)
    {
      s<<a.pos().x<<"\t"<<a.pos().y<<"\t"<<a.type()<<"\t"<<a.prob()<<"\n";
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
      std::cerr<<"foto"<<other.num<<"\n";
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







CortexState SimplestModel::next(CortexState &c)
{
  CortexState n(c);

  unsigned nboxes=c.nboxes;


  double ap=0;
  for (unsigned is=0; is<c.nstates; ++is)
    ap+=a[is]*c.pAstro_[0][is]*c.cAstro[0];


  double jr=(c.damp_[1]-c.damp_[0])/(c.dx_[0]+c.dx_[1]);
  double jl=0;
  n.damp_[0]=c.damp_[0]
      +Dd_*2.0*jr/c.dx_[0]
      -Kd_*c.damp_[0]-ap;
  for (unsigned ix=1; ix<c.damp_.size()-1; ++ix)
    {
      ap=0;
      for (unsigned is=0; is<c.nstates; ++is)
        ap+=a[is]*c.pAstro_[ix][is]*c.cAstro[ix];
      jr=(c.damp_[ix+1]-c.damp_[ix])/(c.dx_[ix]+c.dx_[ix+1]);
      jl=(c.damp_[ix-1]-c.damp_[ix])/(c.dx_[ix]+c.dx_[ix-1]);

      n.damp_[ix]=c.damp_[ix]
          +Dd_*2.0*(jr+jl)/c.dx_[0]
          -Kd_*c.damp_[ix]-ap;

    }
  ap=0;
  for (unsigned is=0; is<c.nstates; ++is)
    ap+=a[is]*c.pAstro_[nboxes][is]*c.cAstro[nboxes];


  jl=(c.damp_[nboxes-1]-c.damp_[nboxes])/(c.dx_[nboxes]+c.dx_[nboxes-1]);
  n.damp_[nboxes]=c.damp_[nboxes]
      +Dd_*2.0*jl/c.dx_[nboxes]
      -(Kd_+ap)*c.damp_[nboxes];


  // mediators

  double vpb=0;
  for (unsigned is=0; is<c.nstates; ++is)
    vpb+=(v[is]-c.med_[0]*b[is])*c.pAstro_[0][is]*c.cAstro[0];;


  jr=(c.med_[1]-c.med_[0])/(c.dx_[0]+c.dx_[1]);
  jl=0;
  n.med_[0]=c.med_[0]
      +Dm_*2.0*jr/c.dx_[0]
      -Km_*c.med_[0]+vpb;
  for (unsigned ix=1; ix<nboxes-1; ++ix)
    {
      vpb=0;
      for (unsigned is=0; is<c.nstates; ++is)
        vpb+=(v[is]-c.med_[ix]*b[is])*c.pAstro_[ix][is]*c.cAstro[ix];;
      jr=(c.med_[ix+1]-c.med_[ix])/(c.dx_[ix]+c.dx_[ix+1]);
      jl=(c.med_[ix-1]-c.med_[ix])/(c.dx_[ix]+c.dx_[ix-1]);

      n.med_[ix]=c.med_[ix]
          +Dm_*2.0*(jr+jl)/c.dx_[ix]
          -Km_*c.med_[ix]+vpb;

    }
  vpb=0;
  for (unsigned is=0; is<c.nstates; ++is)
    vpb+=(v[is]-c.med_[nboxes]*b[is])*c.pAstro_[nboxes][is]*c.cAstro[nboxes];;


  jl=(c.med_[nboxes-1]-c.med_[nboxes])/(c.dx_[nboxes-1]+c.dx_[nboxes]);
  n.med_[nboxes]=c.med_[nboxes]
      +Dm_*2.0*jl/c.dx_[nboxes]
      -Km_*c.med_[nboxes]+vpb;


  std::vector<double> g(0,c.nstates-1);

  for (unsigned ix=0; ix<c.nboxes; ++ix)
    {
      for (unsigned is=0; is<c.nstates-1;++is )
        {
          g[is]=a[is]*q[is]*std::pow(c.damp_[ix],nhd[is])
              /(std::pow(c.damp_[ix],nhd[is])+std::pow(Ed[is],nhd[is]))
              +b[is]*r[is]*std::pow(c.med_[ix],nhm[is])
              /(std::pow(c.med_[ix],nhm[is])+std::pow(Em[is],nhm[is]));
        }
      n.pAstro_[ix][0]=c.pAstro_[ix][1]*fneg[0]-c.pAstro_[ix][0]*(fpos[0]+g[0]);

      for (unsigned is=1; is<c.nstates-1;++is )
        {
          n.pAstro_[ix][is]=+c.pAstro_[ix][is-1]*(fpos[is-1]+g[is-1])
              +c.pAstro_[ix][is+1]*fneg[is]
              -c.pAstro_[ix][is]*(fpos[is]+g[is]+fneg[is-1]);
        }

      n.pAstro_[ix][c.nstates]=+c.pAstro_[ix][c.nstates-1]*(fpos[c.nstates-1]+g[c.nstates-1])
          -c.pAstro_[ix][c.nstates]*(fneg[c.nstates-1]);

    }

  return n;

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


void tissueElement::correctPosition(double dx, double dy)
{
  for (auto& pos:limits_)
    {
      pos.x+=dx;
      pos.y+=dy;
    }
}
