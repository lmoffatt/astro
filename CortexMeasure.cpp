#include "CommandManager.h"
#include "CortexMeasure.h"
#include <string>
#include <list>
#include <iostream>
#include <sstream>
#include <set>
#include <algorithm>
#include <limits>

#include "Parameters.h"

bool TissuePhoto::read(std::string& line, std::istream &s, std::ostream &logs)
{
  std::string name;
  std::stringstream ss(line);
  ss>>name;
  bool out=false;
  num=0;
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

                  astr_.push_back(Astrocyte(sid.str(),rata,dia_,x,y,typo,prob));
                  safeGetline(s,line);
                  ss2.str(line);
                  ss2.clear();
                }
            }
          else if (line.find("celula_")!=std::string::npos)
            {
              safeGetline(s,line);
              safeGetline(s,line);
              safeGetline(s,line);
              std::stringstream ss2;
              ss2.str(line);
              std::string code;
              double d_les,d_tej,d_vas;
              double x,y,prob;
              unsigned typo;
              while (ss2>>code>>d_les>>d_tej>>d_vas>>x>>y>>typo>>prob)
                {
                  std::string idAstro;
                  std::stringstream sid(idAstro);
                  sid<<this->id<<"a"<<astr_.size();

                  astr_.push_back(Astrocyte(sid.str(),rata,dia_,x,y,typo,prob));
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
              limiteFoto_.push_back(LimiteFoto(v));

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
          else if (line.find("limite de lesion")!=std::string::npos)
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
              limiteFoto_.push_back(LimiteFoto(v));

            }
          else if (line.find("limite superior")!=std::string::npos)
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
              ls_=LimiteSuperior(v);

            }

          else if (line.find("limite inferior")!=std::string::npos)
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
              li_=LimiteInferior(v);
              limiteFoto_.back().insert(v.begin(),v.end());

            }
          else if (line.find("limite posterior")!=std::string::npos)
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
              lp_=LimitePosterior(v);
              limiteFoto_.back().insert(v.begin(),v.end());

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
          else if (line.find("ancho lesion")!=line.npos)
            {
              std::stringstream ss(line);
              std::string ancho,lesion;
              ss>>ancho>>lesion>>injury_Width_;
              safeGetline(s,line);

            }
          else if (line.find("area de lesion")!=line.npos)
            {
              std::stringstream ss(line);
              std::string ancho,de,lesion;
              ss>>ancho>>de>>lesion>>injury_Area_;
              safeGetline(s,line);

            }
          else break;

        }
      out=true;
    }
  else
    if (name=="rata" )
      {
        ss>>rata;
        std::stringstream sid(id);
        sid<<"rata"<<rata;
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
                while (ss2>>x>>y>>typo)
                  {
                    prob=1;
                    ss2>>prob;
                    std::string idAstro;
                    std::stringstream sid(idAstro);
                    sid<<this->id<<"a"<<astr_.size();

                    astr_.push_back(Astrocyte(sid.str(),rata,dia_,x,y,typo,prob));
                    safeGetline(s,line);
                    ss2.str(line);
                    ss2.clear();
                  }
              }
            else if (line.find("celula_")!=std::string::npos)
              {
                safeGetline(s,line);
                safeGetline(s,line);
                safeGetline(s,line);
                std::stringstream ss2;
                ss2.str(line);
                std::string code;
                double d_les,d_tej,d_vas;
                double x,y,prob;
                unsigned typo;
                while (ss2>>code>>d_les>>d_tej>>d_vas>>x>>y>>typo>>prob)
                  {
                    std::string idAstro;
                    std::stringstream sid(idAstro);
                    sid<<this->id<<"a"<<astr_.size();

                    astr_.push_back(Astrocyte(sid.str(),rata,dia_,x,y,typo,prob));
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
                limiteFoto_.push_back(LimiteFoto(v));

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
            else if (line.find("limite de lesion")!=std::string::npos)
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
                limiteFoto_.push_back(LimiteFoto(v));

              }
            else if (line.find("limite superior")!=std::string::npos)
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
                ls_=LimiteSuperior(v);
                limiteFoto_.back().insert(v.begin(),v.end());

              }

            else if ((line.find("limite inferior")!=std::string::npos)
                     || (line.find("limite inferrior")!=std::string::npos))
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
                li_=LimiteInferior(v);
                limiteFoto_.back().insert(v.begin(),v.end());

              }
            else if (line.find("limite posterior")!=std::string::npos)
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
                lp_=LimitePosterior(v);
                limiteFoto_.back().insert(v.begin(),v.end());

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
            else if (line.find("ancho lesion")!=line.npos)
              {
                std::stringstream ss(line);
                std::string ancho,lesion;
                ss>>ancho>>lesion>>injury_Width_;
                safeGetline(s,line);

              }
            else if (line.find("area de lesion")!=line.npos)
              {
                std::stringstream ss(line);
                std::string ancho,de,lesion;
                ss>>ancho>>de>>lesion>>injury_Area_;
                safeGetline(s,line);

              }
            else
              break;

          }
        out=true;
      }
  if (this->limiteFoto_.empty())
    {   if (!(ll_.limits().empty()||ls_.limits().empty()||
              lp_.limits().empty()||li_.limits().empty()))
        limiteFoto_.push_back(LimiteFoto(ll_,ls_,lp_,li_));
      else
        logs<<"limite de la foto no definido!!!";

    }


  calculate_distances();
  return out;
}

void TissuePhoto::write(std::ostream &s, bool headers)
{
  if (headers)
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

  if (headers)
    s<<"celulas"<<"\n";
  if (headers)
    {
      s<<astr_.front().getHeader()<<"\n";
      headers=false;
    }
  for (Astrocyte a:astr_)
    {
      a.write(s);
    }
  if (headers)
    s<<"\n";

  for (LimiteVaso v:vasos_)
    {
      if (headers)
        s<<"vaso"<<"\n";
      if (headers)
        s<<"\tCat\t\tdia\trata\t"<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:v.limits())
        s<<"\tLimite vaso\t\t"<<dia_<<"\t"<<rata<<"\t"<<e.x<<"\t"<<e.y<<"\n";
      if (headers)
        s<<"\n";
    }


  if (!ll_.limits().empty())
    {
      if (headers)
        s<<"limite lesion"<<"\n";
      if (headers)
        s<<"\tCat\t\tdia\trata\t"<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:ll_.limits())
        s<<"\tLimite lesio\t\t"<<dia_<<"\t"<<rata<<"\t"<<e.x<<"\t"<<e.y<<"\n";
      if (headers)
        s<<"\n";
    }

  if (!lt_.limits().empty())
    {
      if (headers)
        s<<"limite tejido"<<"\n";
      if (headers)
        s<<"\tCat\t\tdia\trata\t"<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:lt_.limits())
        s<<"\tLimite tejido\t\t"<<dia_<<"\t"<<rata<<"\t"<<e.x<<"\t"<<e.y<<"\n";
      if (headers)
        s<<"\n";
    }
  if (!li_.limits().empty())
    {
      if (headers)
        s<<"limite inferior"<<"\n";
      if (headers)
        s<<"\tCat\t\tdia\trata\t"<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:li_.limits())
        s<<"\tLimite inferior\t\t"<<dia_<<"\t"<<rata<<"\t"<<e.x<<"\t"<<e.y<<"\n";
      if (headers)
        s<<"\n";
    }
  if (!ls_.limits().empty())
    {
      if (headers)
        s<<"limite superior"<<"\n";
      if (headers)
        s<<"\tCat\t\tdia\trata\t"<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:lt_.limits())
        s<<"\tLimite superior\t\t"<<dia_<<"\t"<<rata<<"\t"<<e.x<<"\t"<<e.y<<"\n";
      if (headers)
        s<<"\n";
    }
  for (LimiteFoto v:limiteFoto_)
    {
      if (headers)
        s<<"limite foto"<<"\n";
      if (headers)
        s<<"\tCat\t\tdia\trata\t"<<"X (um)"<<"\t"<<"Y (um)"<<"\n";
      for (auto e:v.limits())
        s<<"\tLimite foto\t\t"<<dia_<<"\t"<<rata<<"\t"<<e.x<<"\t"<<e.y<<"\n";
      if (headers)
        s<<"\n";
    }

}

bool TissuePhoto::align( TissuePhoto &other, std::ostream& logs)
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
      // logs<<"\nfoto"<<other.num<<"\n";
      double sdx=0,sdy=0;
      for (auto aPin:commonPins)
        {
          double dx,dy;
          Pin mine, his;
          mine=pines_[aPin];
          his=other.pines_[aPin];
          dx=mine.limits().front().x-his.limits().front().x;
          logs<<"pin"<<aPin<<"\t";
          logs<<"dx="<<dx<<"\t";
          dxs.push_back(dx);
          dy=mine.limits().front().y-his.limits().front().y;
          logs<<"dy="<<dy<<"\n";


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

  for (auto& v:limiteFoto_)
    {
      v.correctPosition(dx,dy);
    }


  for (auto& a:astr_)
    {
      a.correctPosition(dx,dy);
    }
}

double TissuePhoto::distance_to_neareast_Vaso(const position &p)const
{
  double d=std::numeric_limits<double>::infinity();
  for (const LimiteVaso& v:vasos_)
    {
      double dn=v.distance(p);
      if (dn<d)
        d=dn;
    }
  return d;
}

void TissuePhoto::updateMinMax()
{
  xmin_=+std::numeric_limits<double>::infinity();
  xmax_=-std::numeric_limits<double>::infinity();
  ymin_=+std::numeric_limits<double>::infinity();
  ymax_=-std::numeric_limits<double>::infinity();

  if (limiteFoto_.empty())
    {
      limiteFoto_.back().insert(ll_.limits().begin(),ll_.limits().end());


    }


  for (const LimiteFoto& lf:limiteFoto_ )
    {
      if (lf.xmin_<xmin_)
        xmin_=lf.xmin_;
      if (lf.xmax_>xmax_)
        xmax_=lf.xmax_;
      if (lf.ymin_<ymin_)
        ymin_=lf.ymin_;
      if (lf.ymax_>ymax_)
        ymax_=lf.ymax_;
    }
  if (true){
      const LimiteInferior& lf=li_;
      if (lf.xmin_<xmin_)
        xmin_=lf.xmin_;
      if (lf.xmax_>xmax_)
        xmax_=lf.xmax_;
      if (lf.ymin_<ymin_)
        ymin_=lf.ymin_;
      if (lf.ymax_>ymax_)
        ymax_=lf.ymax_;

    }

  if (true){
      const LimiteSuperior& lf=ls_;
      if (lf.xmin_<xmin_)
        xmin_=lf.xmin_;
      if (lf.xmax_>xmax_)
        xmax_=lf.xmax_;
      if (lf.ymin_<ymin_)
        ymin_=lf.ymin_;
      if (lf.ymax_>ymax_)
        ymax_=lf.ymax_;

    }
  if (true){
      const LimitePosterior& lf=lp_;
      if (lf.xmin_<xmin_)
        xmin_=lf.xmin_;
      if (lf.xmax_>xmax_)
        xmax_=lf.xmax_;
      if (lf.ymin_<ymin_)
        ymin_=lf.ymin_;
      if (lf.ymax_>ymax_)
        ymax_=lf.ymax_;

    }
  if (true){
      const LimiteLesion& lf=ll_;
      if (lf.xmin_<xmin_)
        xmin_=lf.xmin_;
      if (lf.xmax_>xmax_)
        xmax_=lf.xmax_;
      if (lf.ymin_<ymin_)
        ymin_=lf.ymin_;
      if (lf.ymax_>ymax_)
        ymax_=lf.ymax_;

    }

}











void TissueSection::align(std::ostream& logs)
{
  // use first photo as the master

  std::list<unsigned> remaining;
  for (auto it:fotos)
    remaining.push_back(it.first);
  unsigned master=remaining.front();
  remaining.pop_front();
  TissuePhoto* f=fotos[master];

  while (!remaining.empty())
    {
      for (auto it=remaining.begin(); it!=remaining.end(); ++it)
        {
          if (f->align(*fotos[*it],logs))
            {
              it=remaining.erase(it);
            }
        }
    }

}



void TissueSection::merge()
{
  TissuePhoto* foto=new TissuePhoto;

  for (auto it:fotos)
    {
      foto->include(*it.second);
    }
  foto->updateMinMax();


  foto->id="Merged";
  foto->num=0;
  fotos.clear();
  fotos[0]=foto;



}

void TissueSection::distances()
{
  for (auto& f:fotos)
    f.second->calculate_distances();
}

std::vector<CortexMeasure> TissueSection::measure(std::mt19937_64 &mt, std::vector<double> x, double minimalDistanceTissue, double minimalDistanceVaso, std::size_t maxpoints)
{
  std::vector<CortexMeasure> o;
  for (auto f:fotos)
    o.push_back(f.second->measure(mt,id(),dia_,x,minimalDistanceTissue,minimalDistanceVaso,maxpoints));
  return o;
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
    if (ticks_.size()>0)
      return ticks_.size()-1;
    else return 0;
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
      {
        unsigned out=it->second;
        return out-1;
      }
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

CortexMeasure TissuePhoto::measure(std::mt19937_64& mt,std::string id,double dia,std::vector<double> x,double minimal_distance_to_tissue,double minimal_distance_to_vaso
                                   , std::size_t maxpoints)
{


  pointDefined p(x);

  std::vector<double> area(p.num(),0);
  updateMinMax();
  double dArea=(xmax_-xmin_)*(ymax_-ymin_)/maxpoints;
  xud_=std::uniform_real_distribution<double>(xmin_,xmax_);
  yud_=std::uniform_real_distribution<double>(ymin_,ymax_);

  for (std::size_t i=0; i<maxpoints; i++)
    {
      position pos;
      pos.x=xud_(mt);
      pos.y=yud_(mt);

      if (IsInside(pos))
        {
          double distance_to_lession=ll_.distance(pos);
          double distance_to_superior=ls_.distance(pos);
          double distance_to_inferior=li_.distance(pos);
          double distance_to_posterior=lp_.distance(pos);

          double distance_to_vaso=distance_to_neareast_Vaso(pos);
          if ((distance_to_superior>minimal_distance_to_tissue)
              &&(distance_to_inferior>minimal_distance_to_tissue)
              &&(distance_to_posterior>minimal_distance_to_tissue)
              &&(distance_to_vaso>minimal_distance_to_vaso))
            {
              auto i=p.getIndex(distance_to_lession);
              if (i!=p.npos)
                area[i]+=dArea;
            }

        }
    }

  double injlength=0;

  for (std::size_t i=0; i<ll_.limits().size()-1; ++i)
    {
      position p=ll_.limits()[i];
      injlength+=p.distance(ll_.limits()[i+1]);
    }


  std::vector<std::vector<double>> numx(p.num(),std::vector<double>(Astrocyte::numTypes(),0));



  std::vector<std::vector<double>> var(p.num(),std::vector<double>(Astrocyte::numTypes(),0));


  for (Astrocyte a:astr_)
    {
      if ((a.distance_to_tissue()>minimal_distance_to_tissue)
          &&(a.distance_to_vaso()>minimal_distance_to_vaso))
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


  if (injury_Width_==0)
    injury_Width_=injury_Area_/injlength;


  CortexMeasure m(id,dia,
                  100e-6
                  ,rata
                  ,minimal_distance_to_tissue
                  ,minimal_distance_to_vaso
                  ,injury_Width_
                  ,injlength
                  ,p.limits()
                  ,area
                  ,numx
                  ,covar);
  return m;


}









void tissueElement::correctPosition(double dx, double dy)
{
  for (auto& pos:limits_)
    {
      pos.x+=dx;
      pos.y+=dy;
    }
  xmin_+=dx;
  xmax_+=dx;
  ymin_+=dy;
  ymax_+=dy;
}

bool lineIsAtRigth(const position& p, const position& beginLine, const position& endLine)
{
  double xlineCrossPoint=beginLine.x
      +(p.y-beginLine.y)/(endLine.y-beginLine.y)
      *(endLine.x-beginLine.x);
  return (xlineCrossPoint>p.x);

}



bool tissueElement::isInside(const position &p)const
{
  if ((p.x<xmin_)||(p.x>xmax_)||(p.y<ymin_)||(p.y>ymax_))
    return false;
  std::size_t number_of_vertical_crosses=0;
  std::size_t i=0;
  const position* b=&limits()[i];
  ++i;
  const position* e=&limits()[i];
  bool lastAbove=b->y>p.y;
  bool nowAbove=e->y>p.y;
  while (true)
    {
      while (lastAbove==nowAbove)
        {
          b=e;
          lastAbove=nowAbove;
          if (i<limits().size()-1)
            ++i;
          else
            i=0;
          e=&limits()[i];
          nowAbove=e->y>p.y;
          if (i==0)
            break;
        }
      if (lastAbove!=nowAbove)
        {
          double xlineCrossPoint=limits_[i-1].x
              +(p.y-limits_[i-1].y)/(limits()[i].y-limits()[i-1].y)
              *(limits()[i].x-limits()[i-1].x);
          if (xlineCrossPoint>p.x)
            number_of_vertical_crosses++;
          lastAbove=nowAbove;
        }
      if (i==0)
        break;

    }
  return number_of_vertical_crosses % 2==1;

}


std::ostream &Astrocyte::write(std::ostream &s)
{
  s<<"\tAstrocito\t"<<id()<<"\t"<<dia_<<"\t"<<rata_<<"\t";
  s<<pos().x<<"\t"<<pos().y<<"\t"<<type()<<"\t"<<prob()<<"\t";
  if (distance_to_tissue_<std::numeric_limits<double>::infinity())
    s<<distance_to_lession_<<"\t"<<distance_to_tissue_<<"\t"
    <<distance_to_vaso_<<"\t";
  else  if (distance_to_superior_<std::numeric_limits<double>::infinity())
    s<<distance_to_lession_<<"\t"<<distance_to_superior_<<"\t"
    <<distance_to_inferior_<<"\t"<<distance_to_posterior_<<"\t"
    <<distance_to_vaso_<<"\t";
  s<<"\n";
  return s;
}

std::string Astrocyte::getHeader()
{
  std::string ss;
  std::stringstream s(ss);
  s<<"Cat"<<"\t"<<"ID"<<"\t"<<"dia"<<"\t"<<"rata"<<"\t";
  s<<"X (um)"<<"\t"<<"Y (um)"<<"\t"<<"TIPO"<<"\t"<<"PB"<<"\t";
  if (distance_to_tissue_<std::numeric_limits<double>::infinity())
    s<<"d_les (um)"<<"\t"<<"d_tej (um)"<<"\t"<<"d_vaso (um)"<<"\t";
  else if (distance_to_superior_<std::numeric_limits<double>::infinity())
    s<<"d_les (um)"<<"\t"<<"d_sup (um)"<<"\t"<<"d_inf (um)"<<"\t"<<"d_pos (um)"<<"\t"
    <<"d_vaso (um)"<<"\t";

  s<<"\n";
  return s.str();
}

void Astrocyte::calculateDistances(const TissuePhoto &f)
{
  distance_to_lession_=f.ll_.distance(*this);
  distance_to_tissue_=f.lt_.distance(*this);
  distance_to_vaso_=f.distance_to_neareast_Vaso(p_);
  distance_to_superior_=f.ls_.distance(*this);
  distance_to_inferior_=f.li_.distance(*this);
  distance_to_posterior_=f.lp_.distance(*this);

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














const std::vector<double> &Experiment::xpos() const
{
  return m_[0].xpos();
}

std::vector<double> Experiment::x_in_m(double dx,  double sinkLength) const
{

  std::vector<double> o;

 double xpos=0;

  double xmax=m_[0].xpos().back()*1e-6;
  dx=dx*1e-6;
  while (xpos<xmax)
    {
      o.push_back(xpos);
      xpos+=dx;
    }

  double f=2; //increase factor of dt to fill the sink
  double xend=xmax+sinkLength;
  while (xpos<xend)
    {
      o.push_back(xpos);
      dx*=f;
      xpos+=dx;
    }
  return o;
}

double Experiment::h() const
{
  return m_[0].h();
}

std::size_t Experiment::numMeasures() const
{
  return m_.size();
}

std::size_t Experiment::numSimPoints() const
{
  if (isMeasure_)
    return tMeasures_.size();
else
  return tSimulates_.size();
}


double Experiment::tsim() const
{
  if (isMeasure_)
    return tmeas_;
  else
  return tsim_;
}

double Experiment::tMeas(unsigned i) const
{
  return tMeasures_[i];
}


double Experiment::tSimul(unsigned i) const
{
  if (isMeasure_)
    return tMeasures_[i];
else
  return tSimulates_[i];
}


Experiment::Experiment(std::string ide, std::vector<CortexMeasure> mv,const  std::vector<double> &tsimul):
  isMeasure_(true),m_(mv),tMeasures_(mv.size()),tsim_(0),tSimulates_()
{
  setId(ide);
  std::size_t isimul=0;
  for (unsigned i=0; i<m_.size(); ++i)
    {
      tMeasures_[i]=m_[i].dia()*60*60*24;
      if (tMeasures_[i]>tsim_)
        tsim_=tMeasures_[i];
      while (isimul<tsimul.size()&&tsimul[isimul]<tMeasures_[i])
        {
          tSimulates_.push_back(tsimul[isimul]);
          ++isimul;
        }
      tSimulates_.push_back(tMeasures_[i]);
    }


}

const CortexMeasure *Experiment::getMeasure(unsigned i) const
{
  return &m_[i];
}


void updateMinMaxPos(double newx, double newy, double &xmin, double &xmax, double &ymin, double &ymax){
  if (newx<xmin)
    xmin=newx;
  if (newx>xmax)
    xmax=newx;
  if (newy<ymin)
    ymin=newy;
  if (newy>ymax)
    ymax=newy;
}
