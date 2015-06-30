#ifndef CORTEXSTATE
#define CORTEXSTATE
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <limits>
#include <fstream>






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


class CortexMeasure
{
public:
  std::string id()const
  {
    return id_;
  }


  std::vector<double> dx()const
  {
    return dx_;
  }
  double dx(unsigned idx)const
  {
    return dx_[idx];
  }
  double numAstro(unsigned idx)const
  {
    return numAstro_[idx];
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

  std::ostream& write(std::ostream& s)
  {
    s<<id()<<"\n";
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
                std::vector<double> dx
                ,std::vector<std::vector<double>> meanAstro
                ,std::vector<std::vector<std::vector<double>>> covAstro)
    :id_(id),dx_(dx),
      numAstro_(std::vector<double>(meanAstro.size(),0))
    ,meanAstro_(meanAstro),covAstro_(covAstro){

    for (unsigned i=0; i<meanAstro_.size(); ++i
         )
      {
        for (unsigned j=0; j<meanAstro_[i].size();++j)
          numAstro_[i]+=meanAstro_[i][j];
      }

  }

private:
  std::string id_;
  std::vector<double> dx_;
  std::vector<double> numAstro_;
  std::vector<std::vector<double>> meanAstro_;
  std::vector<std::vector<std::vector<double>>> covAstro_;


};




class CortexExperiment
{

  // grilla receptores damp
  std::vector<double> n_;
  std::vector<double>dn_;

  //grilla receptores mediadores
  std::vector<double> m_;
  std::vector<double>dm_;

  // grilla espacial
  std::vector<double> x_;

  std::vector<double>dx_;

  double dt_;
  double teq;
  double tsim_;




  /*

  grilla espacial
  [um]
  0:+50:4000:*1.5:20000

  dt
  [s]
  1

  tiempo equilibrio
  [s]
  7200

  tiempo total
  [s]
  700000

  grilla numero receptores DAMP
  [number of DAMP receptors]
  0:+1:2:*1.5:1E5

  grilla numero receptores mediadores
  [number of mediator receptors]
  0:+1:2:*1.5:1E6

  DAMP inyectado
  [nM]
  400  300  200  0 ...
  */
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
  std::vector<Astrocyte> astr_;

  LimiteTejido lt_;
  LimiteLesion ll_;
  LimiteFoto lf_;

  std::vector<LimiteVaso> vasos_;

  std::map<unsigned,Pin> pines_;



  CortexMeasure* measure(std::string id, std::vector<double> x);
};






class TissueSection
{
public:
  TissueSection (std::string name):
    id_(name){}


  std::string id()const{return id_;}

  std::map<unsigned,TissuePhoto> fotos;


  void align();

  void merge();

  void distances();


  CortexMeasure* measure(std::vector<double> x)
  {
    return fotos.begin()->second.measure(id(),x);
  }



private:
  std::string id_;


};






class CortexState
{
public:
  std::vector<double> x_;
  std::vector<double> dx_;

  std::vector<double> psi_;

  std::vector<double> omega_;
  std::vector<double> N_;

  std::vector<double> n_;

  std::vector<double> M_;
  std::vector<double> m_;


  std::vector<double> dn_;
  std::vector<double> dm_;


  std::vector<std::vector<double> > rho_;

  std::vector<std::vector<std::vector<double>> > rho_i_;

  std::vector<std::vector<std::vector<double>> > rho_j_;




  CortexState(const std::vector<double>& x
              ,const std::vector<double>& dx

              , const std::vector<double>&N
              , const std::vector<double>& n
              , const std::vector<double>& dn
              , const std::vector<double>& M
              ,const std::vector<double>& m
              ,const std::vector<double>& dm
              )
    :
      x_(x)
    ,dx_(dx)
    ,psi_(std::vector<double>(x.size(),0))
    ,omega_(std::vector<double> (x.size(),0))
    ,N_(N),n_(n),M_(M),m_(m),
      dn_(dn)
    , dm_(dm)
    ,rho_(std::vector<std::vector<double> > (x.size(),std::vector<double>(N.size(),0)))
    ,rho_i_(std::vector<std::vector<std::vector<
            double>>>(x.size(),std::vector<std::vector<
                      double>>(N.size(),std::vector<double>(n.size(),0))))

    ,rho_j_(std::vector<std::vector<std::vector<
            double>>>(x.size(),std::vector<std::vector<
                      double>>(N.size(),std::vector<double>(m.size(),0))))
  {}


  void addDamp(const std::vector<double> d)
  {
    for (unsigned i=0; i<d.size(); ++i)
    psi_[i]+=d[i];
  }

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



    class SimParam
    {
    public:
      std::vector<double> n_;
      std::vector<double> m_;
      std::vector<double> x_;
      std::vector<double>dn_;
      std::vector<double>dm_;
      std::vector<double>dx_;
      double dt_;
      double tsim_;
      double teq_;
      std::vector<double> damp_;

    };

    class Param
    {
    public:
      double Dpsi_;
      double Domega_;
      std::vector<double> kon_psi_;
      std::vector<double> kcat_psi_;

      std::vector<double> ksig_omega_;
      std::vector<double> kon_omega_;
      std::vector<double> kcat_omega_;
      std::vector<double> g_rev_;

      std::vector<double> g_max_omega_;
      std::vector<double> h_omega_;
      std::vector<double> n_50_omega_;

      std::vector<double> a_max_omega_;
      std::vector<double> ha_omega_;
      std::vector<double> na_50_omega_;

      std::vector<double> g_max_psi_;
      std::vector<double> h_psi_;
      std::vector<double> n_50_psi_;

      std::vector<double> a_max_psi_;
      std::vector<double> ha_psi_;
      std::vector<double> na_50_psi_;


      std::vector<double> N_;
      std::vector<double> M_;

      std::vector<double> nAstr_;


    };

    SimplestModel() {}
    void loadParameters(Parameters p)=0;
    CortexState nextEuler(const Param &p,const CortexState& c, double dt);

    CortexState init(const Param &p, const SimParam &s);

    CortexState injectDamp(const CortexState& c, double damp);







    std::vector<double> dPsi_dt(const Param &p, const CortexState& c);

    std::vector<double> dOmega_dt(const Param &p, const CortexState &c);

    std::vector<std::vector<std::vector<double> > >
    dRhoOmegaLigand_dt(const Param &p, const CortexState &c);
    std::vector<std::vector<std::vector<double> > >
    dRhoOmegaState_dt(const Param &p, const CortexState &c);

    std::vector<std::vector<std::vector<double> > >
    dRhoPsiLigand_dt(const Param &p, const CortexState &c);

    std::vector<std::vector<std::vector<double> > >
    dRhoPsiState_dt(const Param &p, const CortexState &c);


    double g_omega(const Param &p, unsigned k, double m);
    double a_omega(const Param &p, unsigned k, double m);
    double g_psi(const Param &p, unsigned k, double m);
    double a_psi(const Param &p, unsigned k, double m);

    void run(const SimParam& s,
             const Param& p);



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

