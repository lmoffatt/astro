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
public:
  // grilla espacial
  std::vector<double> x_;

  std::vector<double>dx_;

  double h_;


  double sample_time_;
  double teq_;
  double tsim_;

  void read(std::string& line,std::istream& s);


  void write(std::ostream& s);
 };

class CortexState
{
public:
  std::vector<double> x_;
  std::vector<double> dx_;

  double h_;


  std::vector<double> psi_;

  std::vector<double> omega_;

  std::vector<std::vector<double> > rho_;



  CortexState(const std::vector<double>& x
              ,const std::vector<double>& dx
              ,double h
              ,unsigned numK)
    :
      x_(x)
    ,dx_(dx)
    ,h_(h)
    ,psi_(std::vector<double>(x.size(),0))
    ,omega_(std::vector<double> (x.size(),0))
    ,rho_(std::vector<std::vector<double> > (x.size(),std::vector<double>(numK,0)))
   {}


  void addDamp(const std::vector<double> d)
  {
    for (unsigned i=0; i<d.size(); ++i)
    psi_[i]+=d[i]/dx_[i]/h_/h_/1000.0;
  }

};










class CortexSimulation
{
public:
CortexSimulation(){}

std::ostream& write(std::ostream& s);

std::ostream& write(std::ostream& s, const std::string& var, const std::string& par );



std::string id_;

 std::vector<double> x_;
 std::vector<double> dx_;

 std::vector<double> t_;
 std::vector<double> sdt_;

 std::vector<std::vector<double>> psi_;

 std::vector<std::vector<double>> omega_;

 std::vector<std::vector<std::vector<double>>> rho_;

 CortexSimulation(const CortexState& c,unsigned numSamples):
   x_(c.x_),dx_(c.dx_),
   t_(std::vector<double>(numSamples)),
   sdt_(std::vector<double>(numSamples)),
   psi_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.psi_.size()))),
   omega_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.omega_.size()))),
   rho_(std::vector<std::vector<std::vector<double>>>(
                                                      numSamples,std::vector<std::vector<double>>
                                                      (c.rho_.size(),
                                                       std::vector<double>(
                                                         c.rho_.front().size()))))

 {
     psi_[0]=c.psi_;
     omega_[0]=c.omega_;
     rho_[0]=c.rho_;

 }
 
 
 CortexSimulation(const std::string& id,unsigned numSamples,unsigned numNodes,unsigned numStates):
   id_(id)
   ,x_(std::vector<double>(numNodes)),
   dx_(std::vector<double>(numNodes)),
   t_(std::vector<double>(numSamples)),
   sdt_(std::vector<double>(numSamples)),
   psi_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
   omega_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
   rho_(std::vector<std::vector<std::vector<double>>>(
                                                      numSamples,std::vector<std::vector<double>>
                                                      (numNodes,
                                                       std::vector<double>(
                                                         numStates)))){}
 
 
 
 
 void read(std::string &line, std::istream &s);
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
















  class Parameters
  {
  public:
    double get(const std::string& name)const
    {
      auto it=m_.find(name);
      if (it!=m_.end())
        return it->second.first;
      else
        return std::numeric_limits<double>::quiet_NaN();
    }
    void push_back(const std::string& name, double val,std::string comment)
    {
      m_[name]=std::pair<double,std::string>(val,comment);
    }
    void push_back(const std::string& name, double val)
    {
      m_[name]=std::pair<double,std::string>(val,"");
    }

    unsigned size()const
    {
      return m_.size();
    }

    void read(std::string &line, std::istream &s);



   private:
    std::map<std::string, std::pair<double,std::string>> m_;

  };



  class SimplestModel
  {
  public:




    class Param
    {
    public:
      std::vector<double> damp_;
      double Dpsi_;
      double Domega_;
      double epsilon_;

      double kon_psi_;
      double kcat_psi_;
      double kon_omega_;
      double kcat_omega_;

      std::vector<double> ksig_omega_;
      std::vector<double> g_left_;
      std::vector<double> g_rigth_;

      std::vector<double> g_max_omega_;
      std::vector<double> g_max_psi_;



      std::vector<double> a_;

      std::vector<double> a_omega_;
      std::vector<double> a_psi_;


      std::vector<double> N_;
      std::vector<double> M_;

       double dens_Astr_;

       double dens_Neur_;


    };

     SimplestModel() {}
    CortexState nextEuler(const Param &p,const CortexState& c, double dt) const;

    CortexState init(const Param &p, const CortexExperiment &s) const;

    CortexState injectDamp(const CortexState& c, double damp)const;

    CortexSimulation simulate(const Param &p, const CortexExperiment &sp, double dt)const;






    std::vector<double> dPsi_dt(const Param &p, const CortexState& c) const;

    std::vector<double> dOmega_dt(const Param &p, const CortexState &c)const;





    std::vector<std::vector<double> > dRho_dt(const Param &p, const CortexState &c, bool hasOmega)const;
  };





  class BaseModel
  {public:

    virtual std::string id()const=0;


    virtual Parameters getParameters()const=0;

    virtual void loadParameters(const Parameters& p)=0;

    virtual CortexSimulation run(const CortexExperiment& e,double dt) const=0;


    static BaseModel* create(const Parameters& p);

    virtual ~BaseModel(){}

   private:
    static std::map<double,BaseModel*> models_;
    static std::map<double,BaseModel*> getModels();
  };


  class Model00:public BaseModel
  {
    SimplestModel m;




    class myParameters
    {
    public:
      double D_;
      double epsilon_;
      double kon_;
      double kcat_;
      double g_01_;
      double g_10_;
      double g_23_;
      double g_max_;
      double N_0_;
      double N_2_;
      double N_N_;
      double N_Astr_;
      double N_Neuron_;
      double a_2_;
      double a_factor_;
      double a_max_Neuron_;
      double DAMP_;
    };


    SimplestModel::Param toModelParameters(const myParameters& p)const
    {
      SimplestModel::Param s;
      s.damp_=std::vector<double>(1);
      s.damp_[0]=p.DAMP_;
      s.Dpsi_=p.D_;
      s.Domega_=0;
      s.epsilon_=p.epsilon_;

      s.kon_psi_=p.kon_;
      s.kcat_psi_=p.kcat_;
      s.kon_omega_=0;
      s.kcat_omega_=0;

       s.ksig_omega_=std::vector<double>(7,0);

      s.g_left_=std::vector<double> (7,0);
      s.g_left_[2]=p.g_10_;


      s.g_rigth_=std::vector<double> (7,0);
      s.g_rigth_[1]=p.g_01_;

      s.g_rigth_[3]=p.g_23_;
      s.g_rigth_[4]=p.g_23_;
      s.g_rigth_[5]=p.g_23_;

      s.g_max_omega_=std::vector<double> (7,0);


      s.g_max_psi_=std::vector<double> (7,0);

      s.g_max_psi_[2]=p.g_max_;


      s.a_=std::vector<double> (7,0);
      s.a_[3]=p.a_2_;
      s.a_[4]=p.a_2_*p.a_factor_;
      s.a_[5]=s.a_[4]*p.a_factor_;
      s.a_[6]=s.a_[5]*p.a_factor_;

      s.a_omega_=std::vector<double> (7,0);
      s.a_psi_=std::vector<double> (7,0);

      s.a_psi_[0]=p.a_max_Neuron_;


      s.N_=std::vector<double> (7,0);

      s.N_[0]=p.N_N_;
      s.N_[1]=p.N_0_;
      s.N_[2]=p.N_0_;
      s.N_[3]=p.N_2_;
      s.N_[4]=p.N_2_*1.5;
      s.N_[5]=p.N_2_*3;
      s.N_[6]=p.N_2_*6;



      s.M_=std::vector<double> (7,0);



      s.dens_Astr_=p.N_Astr_;



      s.dens_Neur_=p.N_Neuron_;


      return s;
    }


    myParameters p_;

    // BaseModel interface
  public:
    Model00(){}
    ~Model00(){}
    virtual std::string id() const
    {
      return "Model 0.0";
    }
    static double number()
    {
      return 0;
    }
    virtual Parameters getParameters() const
    {
      Parameters out;
      out.push_back("D",p_.D_);
      out.push_back("epsilon",p_.epsilon_);
      out.push_back("kon",p_.kon_);
      out.push_back("kcat", p_.kcat_);
      out.push_back("g_01",p_.g_01_);
      out.push_back("g_10",p_.g_10_ );
      out.push_back("g_23",p_.g_23_ );
      out.push_back("g_max",p_.g_max_ );
      out.push_back("N_0",p_.N_0_ );
      out.push_back("N_2",p_.N_2_ );
      out.push_back("N_N",p_.N_N_ );
      out.push_back("N_Astr",p_.N_Astr_);
      out.push_back("N_Neuron_",p_.N_Neuron_);
      out.push_back("a_2",p_.a_2_ );
      out.push_back("a_factor",p_.a_factor_ );

      out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
      out.push_back("DAMP",p_.DAMP_);





      return out;

    }
    virtual void loadParameters(const Parameters& p)
    {
      p_.D_=p.get("D");
      p_.epsilon_=p.get("epsilon");
      p_.kon_=p.get("kon");
      p_.kcat_=p.get("kcat");
      p_.g_01_=p.get("g_01");
      p_.g_10_=p.get("g_10");
      p_.g_23_=p.get("g_23");
      p_.g_max_=p.get("g_max");
      p_.N_0_=p.get("N_0");
      p_.N_2_=p.get("N_2");
      p_.N_N_=p.get("N_N");
      p_.a_2_=p.get("a_2");
      p_.DAMP_=p.get("DAMP");
      p_.N_Astr_=p.get("N_Astr");
      p_.N_Neuron_=p.get("N_Neuron");

      p_.a_factor_=p.get("a_factor");
      p_.a_max_Neuron_=p.get("a_max_Neuron");




    }

    virtual CortexSimulation run(const CortexExperiment& e,double dt) const
    {
      return m.simulate(toModelParameters(this->p_),e,dt);


    }

    Model00(const Parameters& p)
    {
      loadParameters(p);
     }


  };



  class Model10:public BaseModel
  {
    SimplestModel m;




    class myParameters
    {
    public:
      double D_;
      double epsilon_;
      double kon_;
      double kcat_;
      double g_01_;
      double g_10_;
      double g_23_;
      double g_max_;
      double N_0_;
      double N_2_;
      double N_N_;
      double N_Astr_;
      double N_Neuron_;
      double a_2_;
      double a_factor_;
      double a_max_Neuron_;
      double DAMP_;
      double k_sig_;
    };


    SimplestModel::Param toModelParameters(const myParameters& p)const
    {
      SimplestModel::Param s;
      s.damp_=std::vector<double>(1);
      s.damp_[0]=p.DAMP_;
      s.Dpsi_=p.D_;
      s.Domega_=p.D_;
      s.epsilon_=p.epsilon_;

      s.kon_psi_=p.kon_;
      s.kcat_psi_=p.kcat_;
      s.kon_omega_=p.kon_;
      s.kcat_omega_=p.kcat_;

       s.ksig_omega_=std::vector<double>(7,0);
       s.ksig_omega_[3]=p.k_sig_;
       s.ksig_omega_[4]=p.k_sig_*1.5;
       s.ksig_omega_[5]=p.k_sig_*3;
       s.ksig_omega_[6]=p.k_sig_*6;

      s.g_left_=std::vector<double> (7,0);
      s.g_left_[2]=p.g_10_;


      s.g_rigth_=std::vector<double> (7,0);
      s.g_rigth_[1]=p.g_01_;

      s.g_rigth_[3]=p.g_23_;
      s.g_rigth_[4]=p.g_23_;
      s.g_rigth_[5]=p.g_23_;

      s.g_max_omega_=std::vector<double> (7,0);


      s.g_max_psi_=std::vector<double> (7,0);

      s.g_max_psi_[2]=p.g_max_;

      s.g_max_omega_[2]=p.g_max_;


      s.a_=std::vector<double> (7,0);
      s.a_[3]=p.a_2_;
      s.a_[4]=p.a_2_*p.a_factor_;
      s.a_[5]=s.a_[4]*p.a_factor_;
      s.a_[6]=s.a_[5]*p.a_factor_;

      s.a_omega_=std::vector<double> (7,0);
      s.a_psi_=std::vector<double> (7,0);

      s.a_psi_[0]=p.a_max_Neuron_;


      s.N_=std::vector<double> (7,0);

      s.N_[0]=p.N_N_;
      s.N_[1]=p.N_0_;
      s.N_[2]=p.N_0_;
      s.N_[3]=p.N_2_;
      s.N_[4]=p.N_2_*1.5;
      s.N_[5]=p.N_2_*3;
      s.N_[6]=p.N_2_*6;



      s.M_=std::vector<double> (7,0);
      s.M_[0]=p.N_N_;
      s.M_[1]=p.N_0_;
      s.M_[2]=p.N_0_;
      s.M_[3]=p.N_2_;
      s.M_[4]=p.N_2_*1.5;
      s.M_[5]=p.N_2_*3;
      s.M_[6]=p.N_2_*6;





      s.dens_Astr_=p.N_Astr_;



      s.dens_Neur_=p.N_Neuron_;


      return s;
    }


    myParameters p_;

    // BaseModel interface
  public:
    Model10(){}
    ~Model10(){}
    virtual std::string id() const
    {
      return "Model 1.0";
    }
    static double number()
    {
      return 1;
    }
    virtual Parameters getParameters() const
    {
      Parameters out;
      out.push_back("D",p_.D_);
      out.push_back("epsilon",p_.epsilon_);
      out.push_back("kon",p_.kon_);
      out.push_back("kcat", p_.kcat_);
      out.push_back("g_01",p_.g_01_);
      out.push_back("g_10",p_.g_10_ );
      out.push_back("g_23",p_.g_23_ );
      out.push_back("g_max",p_.g_max_ );
      out.push_back("N_0",p_.N_0_ );
      out.push_back("N_2",p_.N_2_ );
      out.push_back("N_N",p_.N_N_ );
      out.push_back("N_Astr",p_.N_Astr_);
      out.push_back("N_Neuron_",p_.N_Neuron_);
      out.push_back("a_2",p_.a_2_ );
      out.push_back("a_factor",p_.a_factor_ );

      out.push_back("a_max_Neuron",p_.a_max_Neuron_ );
      out.push_back("DAMP",p_.DAMP_);
      out.push_back("k_sig",p_.k_sig_);






      return out;

    }
    virtual void loadParameters(const Parameters& p)
    {
      p_.D_=p.get("D");
      p_.epsilon_=p.get("epsilon");
      p_.kon_=p.get("kon");
      p_.kcat_=p.get("kcat");
      p_.g_01_=p.get("g_01");
      p_.g_10_=p.get("g_10");
      p_.g_23_=p.get("g_23");
      p_.g_max_=p.get("g_max");
      p_.N_0_=p.get("N_0");
      p_.N_2_=p.get("N_2");
      p_.N_N_=p.get("N_N");
      p_.a_2_=p.get("a_2");
      p_.DAMP_=p.get("DAMP");
      p_.N_Astr_=p.get("N_Astr");
      p_.N_Neuron_=p.get("N_Neuron");

      p_.a_factor_=p.get("a_factor");
      p_.a_max_Neuron_=p.get("a_max_Neuron");
      p_.k_sig_=p.get("k_sig");




    }

    virtual CortexSimulation run(const CortexExperiment& e,double dt) const
    {
      return m.simulate(toModelParameters(this->p_),e,dt);


    }

    Model10(const Parameters& p)
    {
      loadParameters(p);
     }


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

