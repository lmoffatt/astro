#ifndef CORTEXSTATE
#define CORTEXSTATE
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <fstream>
#include <limits>
#include <sstream>

#include "LevenbergMarquardt.h"
#include "BaseClass.h"


#include <random>
struct position
{
  position(double xpos,double ypos):
    x(xpos),y(ypos){}

  position():x(),y(){}
  double x;
  double y;

  double distance(const position& other)const
  {
    return std::sqrt(std::pow(x-other.x,2)+std::pow(y-other.y,2));
  }

};

void updateMinMaxPos(double newx
                     , double newy,
                     double& xmin
                     , double& xmax
                     , double& ymin
                     , double& ymax);

class elementBase
{
public:
  virtual position pos()const=0;
  virtual std::string id()const=0;
  virtual ~elementBase(){}
};

class TissuePhoto;

class Astrocyte: public elementBase
{
public:
  virtual position pos() const{
    return p_;
  }
  virtual std::string id() const
  {
    return id_;
  }

  std::ostream & write(std::ostream &s);
  std::string getHeader();
  void correctPosition(double dx,double dy)
  {
    p_.x+=dx;
    p_.y+=dy;
  }

  void calculateDistances(const TissuePhoto& f);

  unsigned type()const
  {
    return type_;
  }

  double prob()const
  {
    return prob_;
  }


  static unsigned numTypes()
  {
    return 7;
  }

  double distance_to_lession()const
  {
    return distance_to_lession_;
  }

  double distance_to_superior()const
  {
    return distance_to_superior_;
  }
  double distance_to_inferior()const
  {
    return distance_to_inferior_;
  }
  double distance_to_posterior()const
  {
    return distance_to_posterior_;
  }

  double distance_to_tissue()const
  {
    return distance_to_tissue_;
  }
  double distance_to_vaso()const
  {
    return distance_to_vaso_;
  }


  Astrocyte(std::string name,std::size_t rata, double dia, double xpos,double ypos,unsigned type,double prob):
    id_(name)
  , rata_(rata)
  ,dia_(dia)
  ,p_(xpos,ypos)
  ,type_(type)
  ,prob_(prob)
  ,distance_to_lession_(std::numeric_limits<double>::infinity())
  ,distance_to_tissue_(std::numeric_limits<double>::infinity())
  ,distance_to_vaso_(std::numeric_limits<double>::infinity())
  ,distance_to_superior_(std::numeric_limits<double>::infinity())
  ,distance_to_inferior_(std::numeric_limits<double>::infinity())
  ,distance_to_posterior_(std::numeric_limits<double>::infinity())

  {}
private:


  std::string id_;
  std::size_t rata_;
  double dia_;
  position p_;
  unsigned type_;
  double prob_;  // of current type (1-prob) is of previous type,
  double distance_to_lession_;
  double distance_to_tissue_;
  double distance_to_vaso_;
  double distance_to_superior_;
  double distance_to_inferior_;
  double distance_to_posterior_;
};




class tissueElement: public elementBase
{
public:
  virtual position pos() const
  {
    return p;
  }



  virtual std::vector<position>& limits()
  {
    return limits_;
  }

  virtual void insert(std::vector<position>::const_iterator be, std::vector<position>::const_iterator en )
  {
    limits_.insert(limits_.end(),be,en);
    updateMinMax();
  }

  virtual void insert(const std::vector<position>& v )
  {
    double hh=limits_.begin()->distance(*v.begin());
    double hb=limits_.begin()->distance(v.back());
    double bh=limits_.back().distance(*v.begin());
    double bb=limits_.back().distance(v.back());

    std::vector<double> d={hh,hb,bh,bb};
    unsigned imin=0;
    unsigned i;
    //double min=d[imin];
    for (i=1;i<4; ++i)
      if (d[i]<d[imin])
        {
          imin=i;
        }
    switch (imin) {
      case 0:
        limits_.insert(limits_.begin(),v.rbegin(),v.rend());
        break;
      case 1:
        limits_.insert(limits_.begin(),v.begin(),v.end());
        break;
      case 2:
        limits_.insert(limits_.end(),v.rbegin(),v.rend());
        break;
      case 3:
        limits_.insert(limits_.end(),v.begin(),v.end());
        break;
      default:
        break;
      }

    updateMinMax();
  }




  virtual const std::vector<position>& limits()const
  {
    return limits_;
  }


  virtual std::string id() const
  {
    return name_;
  }

  void correctPosition(double dx,double dy);


  double distance(const Astrocyte& a)const
  {
    double d=std::numeric_limits<double>::infinity();
    for (std::size_t i=0; i<limits().size(); ++i)
      {
        d=std::min(a.pos().distance(limits()[i]),d);
      }
    return d;
  }

  double distance(const position& a)const
  {
    double d=std::numeric_limits<double>::infinity();
    for (std::size_t i=0; i<limits().size(); ++i)
      {
        d=std::min(a.distance(limits()[i]),d);
      }
    return d;
  }

  tissueElement(std::string name,
                std::vector<position> limit):
    name_(name)
  ,limits_(limit){
    updateMinMax();
  }




  bool isInside(const position& p)const;


  tissueElement(){}
  double xmax_,xmin_,ymax_,ymin_;

protected:
  void updateMinMax()
  {
    if (!limits_.empty())
      {
        xmax_=limits_[0].x;
        xmin_=limits_[0].x;
        ymax_=limits_[0].y;
        ymin_=limits_[0].y;
        for (const position& pos:limits_)
          updateMinMaxPos(pos.x,pos.y,xmin_,xmax_,ymin_,ymax_);
      }
  }

private:
  std::string name_;
  position p;
  std::vector<position> limits_;




  // elementBase interface
};

class LimiteLesion: public tissueElement
{public:
  LimiteLesion(std::vector<position> lim):
    tissueElement("Limite_Lesion",lim){}

  LimiteLesion(){}

};

class LimiteSuperior: public tissueElement
{public:
  LimiteSuperior(std::vector<position> lim):
    tissueElement("Limite_superior",lim){}

  LimiteSuperior(){}

};

class LimiteInferior: public tissueElement
{public:
  LimiteInferior(std::vector<position> lim):
    tissueElement("Limite_inferior",lim){}

  LimiteInferior(){}

};


class LimitePosterior: public tissueElement
{public:
  LimitePosterior(std::vector<position> lim):
    tissueElement("Limite_posterior",lim){}

  LimitePosterior(){}

};


class LimiteTejido: public tissueElement
{
public:
  LimiteTejido(std::vector<position> lim):
    tissueElement("Limite_Tejido",lim){}

  LimiteTejido(){}

};

class LimiteFoto: public tissueElement
{
public:
  LimiteFoto(std::vector<position> lim):
    tissueElement("Limite_Foto",lim){}

  LimiteFoto(const LimiteLesion& le,const LimiteSuperior& su,
             const LimitePosterior& po,const LimiteInferior& in):
    tissueElement("Limite_Foto",{le.limits()})
  {
    insert(su.limits());
    insert(po.limits());
    insert(in.limits());
  }
  LimiteFoto(){}

};



class Pin:public tissueElement
{
public:
  Pin(std::vector<position> lim):
    tissueElement("Pin",lim){}
  Pin():tissueElement(){}



};

class LimiteVaso: public tissueElement
{
public:
  LimiteVaso(std::vector<position> lim):
    tissueElement("Limite_Vaso",lim){}
  LimiteVaso(){}


};


class CortexMeasure: public BaseObject
{
public:

  // BaseClass interface
public:
  static std::string ClassName(){return "CortexMeasure";}
  virtual std::string myClass() const override
  {
    return ClassName();
  }

  // BaseObject interface
public:
  virtual CortexMeasure *create() const override
  {
    return new CortexMeasure;
  }
  CortexMeasure(){}
  double dia()const
  {
    return dia_;
  }

  double h()const
  {
    return h_;
  }

  double inj_Width() const
  {
    return injWidth_;
  }

  double inj_Area()const
  {
    return inj_Width()*injLength_;
  }

  double inj_Length()const
  {
    return injLength_;
  }

  const std::vector<double>& xpos()const
  {
    return x_;
  }
  double xpos(unsigned idx)const
  {
    return x_[idx];
  }

  const std::vector<double>& numAstro()const
  {
    return numAstro_;
  }

  const std::vector<double>& areaAstro()const
  {
    return measAreaAstro_;
  }

  double numAstro(unsigned idx)const
  {
    return numAstro_[idx];
  }


  bool add(const CortexMeasure& other)
  {
     if ((other.dia_==dia_)&&
         (other.h_==h_)&&
         (other.x_==x_))
       {
         injWidth_+=other.injWidth_*other.rata_.size();
         injWidth_/=1.0+other.rata_.size();
         injLength_+=other.injLength_;

         rata_.insert(rata_.begin(),other.rata_.begin(),other.rata_.end());
         measAreaAstro_+=other.measAreaAstro_;
         meanAstro_+=other.meanAstro_;
         covAstro_+=other.covAstro_;
         for (unsigned i=0; i<meanAstro_.size(); ++i)
           {
             for (unsigned j=0; j<meanAstro_[i].size();++j)
               {
                 numAstro_[i]+=meanAstro_[i][j];
               }
             densAstro_[i]=numAstro_[i]/measAreaAstro_[i];

           }

         return true;
       }
     else {
         return false;
       }



  }

  const std::vector<std::vector<double>>& meanAstro()const
  {
    return meanAstro_;
  }
  const std::vector<double>& meanAstro(unsigned idx)const
  {
    return meanAstro_[idx];
  }
  double meanAstro(unsigned idx,unsigned type)const
  {
    return meanAstro_[idx][type];
  }

  double pAstro(unsigned idx,unsigned type)const
  {
    return meanAstro_[idx][type]/numAstro_[idx];
  }

  double covAstro(unsigned idx, unsigned typei, unsigned typej)const
  {
    return covAstro_[idx][typei][typej];
  }
  double pcovAstro(unsigned idx, unsigned typei, unsigned typej)const
  {
    return covAstro_[idx][typei][typej]/numAstro_[idx];
  }

  std::vector<std::vector<double> > covAstro(unsigned idx)const
  {
    return covAstro_[idx];
  }



  std::ostream& print(std::ostream& s)const
  {
    s<<id()<<"\t"<<dia()<<"\n";
    s<<"limit"<<"\t"<<"tipo 1"<<"\t"<<"tipo 2"<<"\t"<<"tipo 3"<<"\t"<<"tipo 4"<<"\t";
    s<<"tipo 5"<<"\t"<<"tipo 6"<<"\t";
    s<<"se 1"<<"\t"<<"se 2"<<"\t"<<"se 3"<<"\t"<<"se 4"<<"\t"<<"se 5";
    s<<"\t"<<"se 6"<<"\t"<<"se 7"<<"\t";
    s<<"cov 13"<<"\t"<<"cov 34"<<"\t"<<"se 45"<<"\t"<<"se 56"<<"\n";

    for (unsigned ix=0; ix<x_.size(); ++ix)
      {
        s<<x_[ix]<<"\t";
        for (int i=0; i<6; ++i)
          s<<meanAstro(ix,i)<<"\t";
        for (int i=0; i<6; ++i)
          s<<sqrt(covAstro(ix,i,i))<<"\t";
        s<<covAstro(ix,0,2)<<"\t"<<covAstro(ix,2,3)<<"\t"<<covAstro(ix,3,4)<<"\t";
        s<<covAstro(ix,4,5)<<"\n";


      }
    s<<id()<<" ratios \n";
    s<<"limit"<<"\t"<<"numAstro"<<"\t";
    s<<"tipo 1"<<"\t"<<"tipo 2"<<"\t"<<"tipo 3"<<"\t"<<"tipo 4"<<"\t";
    s<<"tipo 5"<<"\t"<<"tipo 6"<<"\t";
    s<<"se 1"<<"\t"<<"se 2"<<"\t"<<"se 3"<<"\t"<<"se 4"<<"\t"<<"se 5";
    s<<"\t"<<"se 6"<<"\t"<<"se 7"<<"\t";
    s<<"cov 13"<<"\t"<<"cov 34"<<"\t"<<"se 45"<<"\t"<<"se 56"<<"\n";

    for (unsigned ix=0; ix<x_.size(); ++ix)
      {

        s<<x_[ix]<<"\t"<<numAstro(ix)<<"\t";
        for (int i=0; i<6; ++i)
          s<<pAstro(ix,i)<<"\t";
        for (int i=0; i<6; ++i)
          s<<sqrt(pcovAstro(ix,i,i))<<"\t";
        s<<pcovAstro(ix,0,2)<<"\t"<<pcovAstro(ix,2,3)<<"\t"<<pcovAstro(ix,3,4)<<"\t";
        s<<pcovAstro(ix,4,5)<<"\n";


      }

    return s;
  }






  CortexMeasure(std::string id,
                double dia,
                double h,
                std::size_t rata,
                double minimalTissueDistance
                ,double minimalVasoDistance
                ,double injWidth
                ,double injLength
                ,std::vector<double> dx
                ,std::vector<double> measAreaAstro
                ,std::vector<std::vector<double>> meanAstro
                ,std::vector<std::vector<std::vector<double>>> covAstro)
    :dia_(dia),h_(h),rata_(1,rata)
    ,minimalTissueDistance_(minimalTissueDistance)
    ,minimalVasoDistance_(minimalVasoDistance)
    ,injWidth_(injWidth)
    ,injLength_(injLength)
    ,x_(dx)
    ,measAreaAstro_(measAreaAstro)
    ,numAstro_(std::vector<double>(meanAstro.size(),0))
    ,meanAstro_(meanAstro),densAstro_(numAstro_),
      covAstro_(covAstro){

    setId(id);
    for (unsigned i=0; i<meanAstro_.size(); ++i)
      {
        for (unsigned j=0; j<meanAstro_[i].size();++j)
          {
            numAstro_[i]+=meanAstro_[i][j];
          }
        densAstro_[i]=numAstro_[i]/measAreaAstro[i];

      }

  }






protected:
  void update() override
  {}



private:
  double dia_;
  double h_; // in um
  std::vector<std::size_t> rata_;
  double minimalTissueDistance_;
  double minimalVasoDistance_;
  double injWidth_;
  double injLength_;
  std::vector<double> x_;
  std::vector<double> measAreaAstro_;
  std::vector<double> numAstro_;
  std::vector<std::vector<double>> meanAstro_;
  std::vector<double> densAstro_;

  std::vector<std::vector<std::vector<double>>> covAstro_;

  // BaseObject interface
public:
  virtual void clear()override
  {
    x_.clear();
    numAstro_.clear();
    meanAstro_.clear();
    measAreaAstro_.clear();
    densAstro_.clear();
    covAstro_.clear();
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"dia",dia_);
    writeField(s,"tissue_width",h_);
    writeField(s,"inj_width",injWidth_);
    writeField(s,"inj_length",injLength_);
    writeField(s,"minimal_tissue_distance",minimalTissueDistance_);
    writeField(s,"minimal_vaso_distance",minimalVasoDistance_);

    writeField(s,"x_pos",x_);
    writeField(s,"area_covered_by_Astrocytes",measAreaAstro_);
    writeField(s,"total_number_of_Astrocytes",numAstro_);

    writeField(s,"number_of_Astrocytes_of_each_type",meanAstro_);
    writeField(s,"density_of_Astrocytes_of_each_type",densAstro_);
    writeField(s,"covariance_of_number_of_Astrocytes_of_each_type",covAstro_);

    return s;
  }

  void writeDataFrameHeader(std::ostream& os)const
  {
    os<<"id"<<"\t";
    os<<"ratas"<<"\t";
    os<<"dia"<<"\t";
    os<<"h"<<"\t";
    os<<"minimalTissueDistance"<<"\t";
    os<<"minimalVasoDistance"<<"\t";
    os<<"injWidth"<<"\t";
    os<<"x_start"<<"\t";
    os<<"x_end"<<"\t";
    os<<"Area"<<"\t";
    os<<"astroType"<<"\t";
    os<<"numAstrocites"<<"\t";
    os<<"densityAstrocites"<<"\t";
    os<<"pAstrocites";

  }
  void writeDataFrame(std::ostream& os, std::size_t idx, std::size_t type) const
  {

    os<<id()<<"\t";
    for (std::size_t ir=0; ir<rata_.size(); ++ir)
       os<<rata_[ir]<<"_";
    os<<"\t";
    os<<dia_<<"\t";
    os<<h_<<"\t";
    os<<minimalTissueDistance_<<"\t";
    os<<minimalVasoDistance_<<"\t";
    os<<injWidth_<<"\t";
    os<<x_[idx]*1e-6<<"\t";
    os<<x_[idx+1]*1e-6<<"\t";
    os<<measAreaAstro_[idx]<<"\t";
    os<<type<<"\t";
    //double simulated_liters=(x_[idx+1]-x_[idx])*1e-6*h()*h()*1000;
    double measured_liters=measAreaAstro_[idx]*1e-12*h()*1000;
    os<<meanAstro_[idx][type]<<"\t";
    os<<meanAstro_[idx][type]/measured_liters<<"\t";
    os<<meanAstro_[idx][type]/numAstro_[idx];

  }


  void writeDataFrame(std::ostream& os) const
  {
    writeDataFrameHeader(os);
    os<<"\n";
    for (std::size_t idx=0; idx<meanAstro_.size(); ++idx)
      for (std::size_t type=0; type<meanAstro_[0].size(); ++type)
        {
          writeDataFrame(os,idx,type);
          os<<"\n";
        }
  }


  bool readBody(std::string& line,std::istream &s, std::ostream& logs) override
  {
    if (!readField(line,s,"dia",dia_,logs))
      {
        logs<<"dia expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"tissue_width",h_,logs))
      {
        logs<<"tissue_width expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"inj_width",injWidth_,logs))
      {
        logs<<"inj_width expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"inj_length",injLength_,logs))
      {
        logs<<"inj_length expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"minimal_tissue_distance",minimalTissueDistance_,logs))
      {
        logs<<"minimal_tissue_distance expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"minimal_vaso_distance",minimalVasoDistance_,logs))
      {
        logs<<"minimal_vaso_distance expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"x_pos",x_,logs))
      {
        logs<<"x_pos expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"area_covered_by_Astrocytes",measAreaAstro_,logs))
      {
        logs<<"area_covered_by_Astrocytes expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"total_number_of_Astrocytes",numAstro_,logs))
      {
        logs<<"total_number_of_Astrocytes expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"number_of_Astrocytes_of_each_type",meanAstro_,logs))
      {
        logs<<"number_of_Astrocytes_of_each_type expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"density_of_Astrocytes_of_each_type",densAstro_,logs))
      {
        logs<<"density_of_Astrocytes_of_each_type expected; found: "<<line<<std::endl;
        return false;
      }
    else if (!readField(line,s,"covariance_of_number_of_Astrocytes_of_each_type",covAstro_,logs))
      {
        logs<<"covariance_of_number_of_Astrocytes_of_each_type expected; found: "<<line<<std::endl;
        return false;
      }
    else  return true;

  }
};








class Experiment: public BaseObject
{
  // BaseClass interface
public:
  static std::string ClassName()
  {
    return "Experiment";
  }

  virtual std::string myClass() const override
  {
    return ClassName();
  }

  // BaseObject interface
public:
  virtual Experiment *create() const override
  {
    return new Experiment;
  }


  const std::vector<double>& xpos() const;
  std::vector<double> x_in_m(double dx, double sinkLength=10e-2) const;

  double h() const;

  std::size_t numMeasures()const;

  double tsim()const;

  double tMeas(unsigned i)const;

  void insert_tSim(double t);

 void push_back(CortexMeasure c);


  Experiment(std::string ide, std::vector<CortexMeasure> mv,const std::vector<double>& tsimul);

  const CortexMeasure* getMeasure(unsigned i)const;

  Experiment(){}


  // BaseObject interface
public:
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"simulation_time",tsim_);
    writeField(s,"time_of_Measures",tMeasures_);
    writeField(s,"time_of_Simulated_Points",tSimulates_);

    writeField(s,"Measures",m_);
    return s;
  }

  void clear()override
  {
    m_.clear();
    tMeasures_.clear();
    tSimulates_.clear();
  }

  virtual bool readBody(std::string& line,std::istream &s, std::ostream& logs) override
  {
    if (!readField(line,s,"simulation_time",tsim_,logs))
      return false;
    else if (!readField(line,s,"time_of_Measures",tMeasures_,logs))
      return false;
    else if (!readField(line,s,"time_of_Simulated_Points",tSimulates_,logs))
      return false;

    else if(!readField(line,s,"Measures",m_,logs))
      return false;
    else
      return true;
  }

  double tSimul(unsigned i) const;
  std::size_t numSimPoints() const;




  void writeDataFrame(std::ostream& os) const
  {
    m_[0].writeDataFrameHeader(os);
    os<<"\n";
    for (auto& m:m_)
      {
        for (std::size_t idx=0; idx<m.meanAstro().size(); ++idx)
          for (std::size_t type=0; type<m.meanAstro()[0].size(); ++type)
            {
              m.writeDataFrame(os,idx,type);
              os<<"\n";
            }
      }
  }



protected:
  void update()override{}




private:
  std::vector<CortexMeasure> m_;
  std::vector<double> tMeasures_;
  double tsim_;
  std::vector<double> tSimulates_;




};



class CortexExperiment
{
public:
  // grilla espacial
  std::vector<double> x_;

  std::vector<double>dx_;

  double h_;


  double sample_time_;
  double teq_;
  double tsim_;

  void read(std::string& line,std::istream& s);


  void write(std::ostream& s)const
  {
    s<<"cortex experiment"<<"\n";

  }


  
};












class TissuePhoto
{
public:

  bool read(std::string& line, std::istream& s, std::ostream& logs);


  void write(std::ostream& s, bool headers=false);


  /// it makes a copy of other where the position
  /// of the elements is based on the coordinate s
  /// system of this
  bool align(TissuePhoto& other, std::ostream &logs);

  void correctPosition(double dx, double dy);


  void include(TissuePhoto& other)
  {
    astr_.insert(astr_.end(),other.astr_.begin(),other.astr_.end());
    lt_.insert(other.lt_.limits().begin(),other.lt_.limits().end());
    ll_.insert(other.ll_.limits().begin(),other.ll_.limits().end());

    limiteFoto_.insert(
          limiteFoto_.end(),other.limiteFoto_.begin(),other.limiteFoto_.end());

    vasos_.insert(vasos_.end(),other.vasos_.begin(),other.vasos_.end());
    pines_.insert(other.pines_.begin(),other.pines_.end());

    if (other.injury_Width_>0)
      injury_Width_=other.injury_Width_;
  }


  void calculate_distances()
  {
    for (Astrocyte& a:astr_)
      {
        a.calculateDistances(*this);
      }

  }

  double distance_to_neareast_Vaso(const position &p) const;



  void updateMinMax();


  unsigned num;
  std::size_t rata;

  std::string id;
  double dia_;
  std::vector<Astrocyte> astr_;

  LimiteTejido lt_;
  LimiteLesion ll_;
  LimiteSuperior ls_;
  LimiteInferior li_;
  LimitePosterior lp_;

  std::vector<LimiteFoto> limiteFoto_;
  double xmin_,xmax_,ymin_,ymax_;

  std::uniform_real_distribution<double> xud_, yud_;

  double injury_Width_=0;
  double injury_Area_=0;


  std::vector<LimiteVaso> vasos_;

  std::map<unsigned,Pin> pines_;

  bool IsInside(const position& p)const
  {
    for (const LimiteVaso& v:vasos_)
      if (v.isInside(p))
        return false;

    for (const LimiteFoto& f:limiteFoto_)
      if (f.isInside(p))
        return true;
    return false;
  }


  CortexMeasure measure(std::mt19937_64 &mt, std::string id, double dia, std::vector<double> x, double minimal_distance_to_tissue=0, double minimal_distance_to_vaso=0, std::size_t maxpoints=1E4);
};





class TissueSection : public BaseObject
{
public:
  TissueSection (std::string name,double dia):
    id_(name),
    dia_(dia)

  {}
  TissueSection()=default;

  ~TissueSection(){}



  void write(std::ostream& s);


  std::string id()const{return id_;}

  std::map<std::size_t,TissuePhoto*> fotos;


  void align(std::ostream &logs);

  void merge();

  void distances();


  std::vector<CortexMeasure> measure(std::mt19937_64& mt,std::vector<double> x,double minimalDistanceTissue,double minimalDistanceVaso,std::size_t maxpoints);

  template<class CM>
  void measure(CM* cm,std::mt19937_64& mt,std::vector<double> x,double minimalDistanceTissue,double minimalDistanceVaso,std::size_t maxpoints)
  {
    {
      for (auto& f:fotos)
        cm->push_back(f.second->measure(mt,id(),dia_,x,minimalDistanceTissue,minimalDistanceVaso,maxpoints));
    }

  }


private:
  std::string id_;
  double dia_;

  // BaseClass interface
public:
  virtual std::string myClass() const override
  {
    return "TissueSections";
  }

  // BaseObject interface
public:
  virtual BaseObject *create() const override
  {
    return new TissueSection();
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"id",id_);
    writeField(s,"dia",dia_);
    s<<fotos.begin()->second->astr_.front().getHeader()<<"\n";
    for (auto e:fotos)
      e.second->write(s,false);
    return s;

  }
  virtual void clear() override{}
  virtual bool readBody(std::string &line, std::istream &s, std::ostream &logs) override
  {
    if (!readField(line,s,"id",id_,logs)) return false;
    if (!readField(line,s,"dia",dia_,logs)) return false;
    if (!readField(line,s,"ratones",fotos,logs)) return false;
    return true;
  }

protected:
  virtual void update() override{}
};




















inline
std::vector<double>& operator+=(
    std::vector<double> & x
    ,const std::vector<double> & y
    )
{

  for (unsigned i=0; i<x.size(); ++i)
    x[i]+=y[i];

  return x;
}





inline
bool addStep(
    std::vector<double> & x
    ,const std::vector<double> & y, double h
    )
{

  bool result=true;
  for (unsigned i=0; i<x.size(); ++i)
    {
      x[i]+=y[i]*h;
      if (x[i]<0)
        {
          result=false;
        }
    }
  return result;
  }



inline
std::vector<double> operator+(
    const std::vector<double> & x
    ,const std::vector<double> & y
    )
{
  std::vector<double> o(x);
  o+=y;

  return o;
}


inline
std::vector<double>& operator-=(
    std::vector<double> & x
    ,const std::vector<double> & y
    )
{

  for (unsigned i=0; i<x.size(); ++i)
    x[i]-=y[i];

  return x;
}


template<typename T>
std::vector<T>& initialize(std::vector<T>& o, T x)
{
  for (auto& a:o)
    a=x;
  return o;
}



inline
std::vector<double> operator-(
    const std::vector<double> & x
    ,const std::vector<double> & y
    )
{
  std::vector<double> o(x);
  o-=y;

  return o;
}










inline
std::vector<double>& operator*=(
    std::vector<double> & x
    ,double t)
{

  for (unsigned i=0; i<x.size(); ++i)
    x[i]*=t;

  return x;
}

inline
std::vector<double> operator*(
    const std::vector<double> & x
    ,double t)
{

  std::vector<double> o(x);
  o*=t;
  return o;
}


inline
std::vector<std::vector<double>>& operator+=(
    std::vector<std::vector<double>> & x
    ,const std::vector<std::vector<double>> & y
    )
{

  for (unsigned i=0; i<x.size(); ++i)
    for (unsigned j=0; j<x.front().size(); ++j)
      x[i][j]+=y[i][j];

  return x;
}


inline
bool addMStep(
    std::vector<std::vector<double>> & x
    ,const std::vector<std::vector<double>> & y
    ,double h)
{
  bool result=true;
  for (unsigned i=0; i<x.size(); ++i)
    {
    for (unsigned j=0; j<x.front().size(); ++j)
      {
        x[i][j]+=y[i][j]*h;
        if ((x[i][j])<0)
          result=false;
      }
    }
  return result;
}





inline
std::vector<std::vector<double>>& operator*=(
    std::vector<std::vector<double>> & x
    ,double t)
{

  for (unsigned i=0; i<x.size(); ++i)
    for (unsigned j=0; j<x.front().size(); ++j)
      x[i][j]*=t;

  return x;
}

inline
std::vector<std::vector<double>> operator*(
    const std::vector<std::vector<double>> & x
    ,double t)
{

  std::vector<std::vector<double>> o(x);
  o*=t;
  return o;
}


inline
std::vector<std::vector<double>> operator*(
    std::vector<std::vector<double>> x
    ,double t)
{

  std::vector<std::vector<double>> o(x);
  o*=t;
  return o;
}




inline
std::vector<std::vector<double>> operator+(
    const std::vector<std::vector<double>> & x
    ,const std::vector<std::vector<double>> & y
    )
{
  std::vector<std::vector<double>> o(x);
  o+=y;
  return o;
}



inline
std::vector<std::vector<std::vector<double>>>& operator+=(
    std::vector<std::vector<std::vector<double>>> & x
    ,const std::vector<std::vector<std::vector<double>>> & y
    )
{

  for (unsigned i=0; i<x.size(); ++i)
    for (unsigned j=0; j<x.front().size(); ++j)
      for (unsigned k=0; k<x.front().front().size(); ++k)
        x[i][j][k]+=y[i][j][k];

  return x;
}

inline
std::vector<std::vector<std::vector<double>>> operator+(
    const std::vector<std::vector<std::vector<double>>> & x
    ,const std::vector<std::vector<std::vector<double>>> & y
    )
{
  std::vector<std::vector<std::vector<double>>> o(x);
  o+=y;
  return o;
}

inline
std::vector<std::vector<std::vector<double>>>& operator*=(
    std::vector<std::vector<std::vector<double>>>& x,double t)
{
  for (unsigned i=0; i<x.size(); ++i)
    for (unsigned j=0; j<x.front().size(); ++j)
      for (unsigned k=0; k<x.front().front().size(); ++k)
        x[i][j][k]*=t;

  return x;

}

inline
std::vector<std::vector<std::vector<double>>> operator*(
    const std::vector<std::vector<std::vector<double>>> & x
    ,double t
    )
{
  std::vector<std::vector<std::vector<double>>> o(x);
  o*=t;
  return o;
}



inline
std::vector<std::vector<double>>
sum(const std::vector<std::vector<std::vector<double>>>& x)
{
  std::vector<std::vector<double>> o(x.size(),std::vector<double>(x.front().size(),0));
  for (unsigned i=0; i<x.size(); ++i)
    for (unsigned j=0; j<x.front().size(); ++j)
      for (unsigned k=0; k<x.front().front().size(); ++k)
        o[i][j]+=x[i][j][k];

  return o;


}



#endif // CORTEXSTATE

