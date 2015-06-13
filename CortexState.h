#ifndef CORTEXSTATE
#define CORTEXSTATE
#include <vector>
#include <cmath>
#include <map>
#include <string>

struct position
{
  position(double xpos,double ypos):
    x(xpos),y(ypos){}

  position():x(),y(){}
  double x;
  double y;


};

class elementBase
{
public:
  virtual position pos()const=0;
  virtual std::string id()const=0;
  virtual ~elementBase(){}
};


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

    void correctPosition(double dx,double dy)
    {
      p_.x+=dx;
      p_.y+=dy;
    }

    unsigned type()const
    {
      return type_;
    }

    double prob()const
    {
      return prob_;
    }

    Astrocyte(std::string name,double xpos,double ypos,unsigned type,double prob):
      id_(name)
      ,p_(xpos,ypos)
    ,type_(type)
    ,prob_(prob){}
private:


    std::string id_;
    position p_;
    unsigned type_;
    double prob_;  // of current type (1-prob) is of previous type,



};




class tissueElement: public elementBase
{
public:
  virtual position pos() const
  {
    return p;
  }



  virtual std::vector<position> limits()
  {
    return limits_;
  }

  virtual std::string id() const
  {
    return name_;
  }

  void correctPosition(double dx,double dy);

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

class TissuePhoto
{
public:


  void read(std::string& line,std::istream& s);


  void write(std::ostream& s);


  /// it makes a copy of other where the position
  /// of the elements is based on the coordinate s
  /// system of this
  bool align( TissuePhoto& other);

  void correctPosition(double dx, double dy);


  unsigned num;
  std::string id;
  std::vector<Astrocyte> astr_;

  LimiteTejido lt_;
  LimiteLesion ll_;

  std::vector<LimiteVaso> vasos_;

  std::map<unsigned,Pin> pines_;



};

class TissueSection
{
public:
  TissueSection (std::string name):
    id_(name){}


  std::string id()const{return id_;}

  std::map<unsigned,TissuePhoto> fotos;


  void align();

private:
  std::string id_;


};





class CortexState
{


public:
  unsigned nstates;
  unsigned nboxes;
  double dt_;
  std::vector<double> dx_;
  std::vector<double> cAstro;

  std::vector<std::vector<double> > pAstro_;

  std::vector<double> damp_;
  std::vector<double> med_;

};


class Parameters
{

};


class CortexModelBase
{
public:
  virtual void loadParameters(Parameters p)=0;
  virtual CortexState next(CortexState& aCortex)=0;

  virtual ~CortexModelBase(){}

};

class SimplestModel: public CortexModelBase
{
public:
  SimplestModel() {}
  void loadParameters(Parameters p)=0;
  CortexState next(CortexState& c);
private:
  double Dd_;
  double Dm_;
  double Kd_;
  double Km_;

  std::vector<double> a,b,v;

  std::vector<double> fpos,fneg;

  std::vector<double> nhd,nhm,q,r,Ed,Em;

};





#endif // CORTEXSTATE

