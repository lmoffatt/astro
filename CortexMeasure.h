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

  Astrocyte(std::string name,double xpos,double ypos,unsigned type,double prob):
    id_(name)
  ,p_(xpos,ypos)
  ,type_(type)
  ,prob_(prob)
  ,distance_to_lession_(std::numeric_limits<double>::infinity())
  ,distance_to_tissue_(std::numeric_limits<double>::infinity()){}
private:


  std::string id_;
  position p_;
  unsigned type_;
  double prob_;  // of current type (1-prob) is of previous type,
  double distance_to_lession_;
  double distance_to_tissue_;
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


  tissueElement(std::string name,
                std::vector<position> limit):
    name_(name)
  ,limits_(limit){
  }



  tissueElement(){}
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

  const std::vector<double>& dx()const
  {
    return dx_;
  }
  double dx(unsigned idx)const
  {
    return dx_[idx];
  }

  const std::vector<double>& numAstro()const
  {
    return numAstro_;
  }

  double numAstro(unsigned idx)const
  {
    return numAstro_[idx];
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

    for (unsigned ix=0; ix<dx_.size(); ++ix)
      {
        s<<dx_[ix]<<"\t";
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

    for (unsigned ix=0; ix<dx_.size(); ++ix)
      {

        s<<dx_[ix]<<"\t"<<numAstro(ix)<<"\t";
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
                std::vector<double> dx
                ,std::vector<std::vector<double>> meanAstro
                ,std::vector<std::vector<std::vector<double>>> covAstro)
    :dia_(dia),h_(h),dx_(dx),
      numAstro_(std::vector<double>(meanAstro.size(),0))
    ,meanAstro_(meanAstro),covAstro_(covAstro){

    setId(id);
    for (unsigned i=0; i<meanAstro_.size(); ++i)
      {
        for (unsigned j=0; j<meanAstro_[i].size();++j)
          numAstro_[i]+=meanAstro_[i][j];
      }

  }






protected:
  void update()
  {}



private:
  double dia_;
  double h_; // in um
  std::vector<double> dx_;
  std::vector<double> numAstro_;
  std::vector<std::vector<double>> meanAstro_;
  std::vector<std::vector<std::vector<double>>> covAstro_;

  // BaseObject interface
public:
  virtual void clear()override
  {
    dx_.clear();
    numAstro_.clear();
    meanAstro_.clear();
    covAstro_.clear();
  }
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"dia",dia_);
    writeField(s,"tissue_width",h_);
    writeField(s,"x_pos",dx_);
    writeField(s,"total_number_of_Astrocytes",numAstro_);
    writeField(s,"number_of_Astrocytes_of_each_type",meanAstro_);
    writeField(s,"covariance_of_number_of_Astrocytes_of_each_type",covAstro_);

    return s;
  }
  bool readBody(std::string& line,std::istream &s) override
  {
    if (
        readField(line,s,"dia",dia_)&&
        readField(line,s,"tissue_width",h_)&&
        readField(line,s,"x_pos",dx_)&&
        readField(line,s,"total_number_of_Astrocytes",numAstro_)&&
        readField(line,s,"number_of_Astrocytes_of_each_type",meanAstro_)&&
        readField(line,s,"covariance_of_number_of_Astrocytes_of_each_type",covAstro_)
        )
      return true;
    else
      return false;

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


  const std::vector<double>& dx() const;
  std::vector<double> x_in_m() const;

  double h() const;

  std::size_t numMeasures()const;

  double tsim()const;

  double tMeas(unsigned i)const;





  Experiment(std::string ide, std::vector<CortexMeasure> mv);

  const CortexMeasure* getMeasure(unsigned i)const;

  Experiment(){}






private:
  std::vector<CortexMeasure> m_;
  std::vector<double> tMeasures_;
  double tsim_;


  // BaseObject interface
public:
  virtual std::ostream &writeBody(std::ostream &s) const override
  {
    writeField(s,"simulation_time",tsim_);
    writeField(s,"time_of_Measures",tMeasures_);
    writeField(s,"Measures",m_);
    return s;
  }

  void clear()override
  {
    m_.clear();
    tMeasures_.clear();
  }

  virtual bool readBody(std::string& line,std::istream &s) override
  {
       if (!readField(line,s,"simulation_time",tsim_))
         return false;
        else if (!readField(line,s,"time_of_Measures",tMeasures_))
         return false;
       else if(!readField(line,s,"Measures",m_))
         return false;
       else
      return true;
    }

protected:
  void update(){}

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


  void read(std::string& line,std::istream& s);


  void write(std::ostream& s);


  /// it makes a copy of other where the position
  /// of the elements is based on the coordinate s
  /// system of this
  bool align(TissuePhoto& other);

  void correctPosition(double dx, double dy);


  void include(TissuePhoto& other)
  {
    astr_.insert(astr_.end(),other.astr_.begin(),other.astr_.end());
    lt_.limits().insert(lt_.limits().end(),other.lt_.limits().begin(),other.lt_.limits().end());
    ll_.limits().insert(ll_.limits().end(),other.ll_.limits().begin(),other.ll_.limits().end());
    lf_.limits().insert(lf_.limits().end(),other.lf_.limits().begin(),other.lf_.limits().end());

    vasos_.insert(vasos_.end(),other.vasos_.begin(),other.vasos_.end());
    pines_.insert(other.pines_.end(),other.pines_.end());
  }


  void calculate_distances()
  {
    for (Astrocyte& a:astr_)
      {
        a.calculateDistances(*this);
      }

  }



  unsigned num;
  std::string id;
  double h_;
  std::vector<Astrocyte> astr_;

  LimiteTejido lt_;
  LimiteLesion ll_;
  LimiteFoto lf_;

  std::vector<LimiteVaso> vasos_;

  std::map<unsigned,Pin> pines_;



  CortexMeasure* measure(std::string id, double dia, std::vector<double> x);
};






class TissueSection
{
public:
  TissueSection (std::string name,double dia):
    id_(name),
    dia_(dia)

  {}


  std::string id()const{return id_;}

  std::map<unsigned,TissuePhoto> fotos;


  void align();

  void merge();

  void distances();


  CortexMeasure* measure(std::vector<double> x)
  {
    return fotos.begin()->second.measure(id(),dia_,x);
  }



private:
  std::string id_;
  double dia_;


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

