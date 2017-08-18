#ifndef CORTEXSIMULATION
#define CORTEXSIMULATION


#include <iomanip>
#include <vector>
#include <map>
#include <limits>
#include <string>
#include <fstream>

#include "Parameters.h"
#include "MatrixBanded.h"

template<class C>
struct Di: public std::pair<C,std::vector<C>>
{
  typedef   std::pair<C,std::vector<C>> P;
  Di(std::size_t n):std::pair<C,std::vector<C>>
  {{},std::vector<C>(n)}{}


  Di(const C& e,const std::vector<C>& v): P(e,v){}


};



class CortexState
{
private:
  static double d(double x,double y, double unit)
  {
    if ((x>0)&&(y>0))
      return std::abs((x-y)/unit);
    else if((x>=0)&&(y>=0))
      return 0;
    else
      return std::numeric_limits<double>::infinity();
  }

  static double d(const std::vector<double>& x, const std::vector<double>& y)
  {
    double out=0;
    double max=x[0];
    for (std::size_t i=0; i<x.size(); ++i)
      {
        max=std::max(max,std::max(x[i],y[i]));
      }

    for (std::size_t i=0; i<x.size(); ++i)
      {
        double e=d(x[i],y[i],max);
        if (e>out)
          out=e;
      }
    return out;
  }
  static double d(const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& y)
  {
    double out=0;
    for (std::size_t i=0; i<x.size(); ++i)
      {
        out=std::max(d(x[i],y[i]),out);
        if (!std::isfinite(out)) return out;
      }
    return out;
  }


public:

  double maxlogdistance(const CortexState& other)const
  {
    if (!isValid_||!other.isValid_)
      return std::numeric_limits<double>::infinity();
    else
      {
        double out=d(rho_,other.rho_);
        out=std::max(d(psi_T_,other.psi_T_),out);
        out=std::max(d(psi_B_,other.psi_B_),out);
        if (hasOmega_)
          {
            out=std::max(d(omega_T_,other.omega_T_),out);
            out=std::max(d(omega_B_,other.omega_B_),out);
          }
        return out;
      }
  }


  struct Der
  {
    std::vector<std::vector<double>> dRho;
    std::vector<double> dPsi;
    std::vector<double> dOmega;
    bool hasOmega;

    std::vector<double> toVector()const
    {
      if (hasOmega)
        {
          std::size_t n=dPsi.size();
          std::size_t m=(2+dRho[0].size());
          std::vector<double> v(n*m);
          for (std::size_t i=0; i<dPsi.size(); ++i)
            {
              v[i*m]=dPsi[i];
              v[i*m+1]=dOmega[i];
              for (std::size_t j=0; j<m-2; ++j)
                v[i*m+2+j]=dRho[i][j];

            }
          return v;
        }
      else
        {
          std::size_t n=dPsi.size();
          std::size_t m=(1+dRho[0].size());
          std::vector<double> v(n*m);
          for (std::size_t i=0; i<dPsi.size(); ++i)
            {
              v[i*m]=dPsi[i];
              for (std::size_t j=0; j<m-1; ++j)
                v[i*m+1+j]=dRho[i][j];

            }
          return v;
        }

    }

    void fromVector(const std::vector<double>& v)
    {

      if (hasOmega)
        {
          std::size_t m=2+dRho[0].size();
          assert (v.size()==dPsi.size()*m);
          for (std::size_t i=0; i<dPsi.size(); ++i)
            {
              dPsi[i]=v[i*m];
              dOmega[i]=v[i*m+1];
              for (std::size_t j=0; j<m-2; ++j)
                dRho[i][j]=v[i*m+2+j];
            }
        }
      else
        {
          std::size_t m=1+dRho[0].size();
          assert (v.size()==dPsi.size()*m);
          for (std::size_t i=0; i<dPsi.size(); ++i)
            {
              dPsi[i]=v[i*m];
              for (std::size_t j=0; j<m-1; ++j)
                dRho[i][j]=v[i*m+1+j];

            }
        }

    }
  };

  struct dDer {

    dDer(std::size_t nx, std::size_t nks, bool hasOmega)
      :nx_(nx), nks_(nks),m_(hasOmega? 2+nks:1+nks),hasOmega_(hasOmega),
        d_(hasOmega? MatrixBanded((2+nks)*nx,2+nks,2*nks+3):MatrixBanded((1+nks)*nx,1+nks,1+2*nks))
    {}

    double& dpsi_dpsi(std::size_t ix, std::size_t jx)
    {
      return d_(psi_ix(ix),psi_ix(jx));
    }
    double const& dpsi_dpsi(std::size_t ix, std::size_t jx)const
    {
      return d_(psi_ix(ix),psi_ix(jx));
    }



    double& dpsi_domega(std::size_t ix, std::size_t jx)
    {
      return d_(psi_ix(ix),omega_ix(jx));
    }

    const double& dpsi_domega(std::size_t ix, std::size_t jx)const
    {
      return d_(psi_ix(ix),omega_ix(jx));
    }

    double& dpsi_drho(std::size_t ix, std::size_t jx, std::size_t ik)
    {
      return d_(psi_ix(ix),rho_ix_ik(jx,ik));
    }

    const double& dpsi_drho(std::size_t ix, std::size_t jx, std::size_t ik) const
    {
      return d_(psi_ix(ix),rho_ix_ik(jx,ik));
    }


    double& domega_dpsi(std::size_t ix, std::size_t jx)
    {
      return d_(omega_ix(ix),psi_ix(jx));
    }



    double const& domega_dpsi(std::size_t ix, std::size_t jx)const
    {
      return d_(omega_ix(ix),psi_ix(jx));
    }


    double& domega_domega(std::size_t ix, std::size_t jx)
    {
      return d_(omega_ix(ix),omega_ix(jx));
    }

    const double& domega_domega(std::size_t ix, std::size_t jx) const
    {
      return d_(omega_ix(ix),omega_ix(jx));
    }


    double& domega_drho(std::size_t ix, std::size_t jx, std::size_t ik)
    {
      return d_(omega_ix(ix),rho_ix_ik(jx,ik));
    }

    const double& domega_drho(std::size_t ix, std::size_t jx, std::size_t ik) const
    {
      return d_(omega_ix(ix),rho_ix_ik(jx,ik));
    }

    double& drho_dpsi(std::size_t ix, std::size_t ik, std::size_t jx)
    {
      return d_(rho_ix_ik(ix,ik),omega_ix(jx));
    }

    const double& drho_dpsi(std::size_t ix, std::size_t ik, std::size_t jx) const
    {
      return d_(rho_ix_ik(ix,ik),omega_ix(jx));
    }

    double& drho_domega(std::size_t ix, std::size_t ik, std::size_t jx)
    {
      return d_(rho_ix_ik(ix,ik),omega_ix(jx));
    }

    const double& drho_domega(std::size_t ix, std::size_t ik, std::size_t jx) const
    {
      return d_(rho_ix_ik(ix,ik),omega_ix(jx));
    }

    double& drho_drho(std::size_t ix, std::size_t ik, std::size_t jx, std::size_t jk)
    {
      return d_(rho_ix_ik(ix,ik),rho_ix_ik(jx,jk));
    }

    std::size_t psi_ix(std::size_t ix) const
    {
      return ix*m_;
    }
    std::size_t omega_ix(std::size_t ix) const
    {
      return ix*m_+1;
    }
    std::size_t rho_ix_ik(std::size_t ix, std::size_t ik) const
    {
      if (hasOmega_)
       return ix*m_+2+ik;
      else
        return ix*m_+1+ik;
    }

    MatrixBanded I_dt(double dt) const
    {
      MatrixBanded out(d_);
      out*=-dt/2;
      for (std::size_t i=0; i<out.n(); ++i)
        out(i,i)+=1.0;
      return out;
    }

   private:
      std::size_t nx_;
      std::size_t nks_;
      std::size_t m_;
      bool hasOmega_;
      MatrixBanded d_;



  };


  bool isValid_;
  bool hasOmega_;
  std::vector<double> x_;
  std::vector<double> dx_;

  double h_;

  double epsilon_;

  std::vector<double> psi_T_;

  std::vector<double> psi_B_;

  std::vector<double> omega_T_;

  std::vector<double> omega_B_;

  std::vector<std::vector<double> > rho_;



  CortexState()=default;


  ~CortexState(){

  }


  CortexState(const std::vector<double>& x
              ,const std::vector<double>& dx
              ,double h
              ,double epsilon
              ,unsigned numK
              ,bool hasOmega)
    :
      isValid_(true)
    ,hasOmega_(hasOmega)
    , x_(x)
    ,dx_(dx)
    ,h_(h)
    ,epsilon_(epsilon)
    ,psi_T_(std::vector<double>(x.size(),0))
    ,psi_B_(std::vector<double>(x.size(),0))
    ,omega_T_(hasOmega ? std::vector<double> (x.size(),0):std::vector<double> (0))
    ,omega_B_(hasOmega ? std::vector<double> (x.size(),0):std::vector<double>(0))
    ,rho_(std::vector<std::vector<double> > (x.size(),std::vector<double>(numK,0)))
  {}

  CortexState(const std::vector<double>& x
              ,double h
              ,double epsilon
              ,unsigned numK
              ,bool hasOmega)
    :
      isValid_(true)
    , hasOmega_(hasOmega)
    , x_(x)
    ,dx_(x.size())
    ,h_(h)
    ,epsilon_(epsilon)
    ,psi_T_(std::vector<double>(x.size(),0))
    ,psi_B_(std::vector<double>(x.size(),0))
    ,omega_T_(std::vector<double> (x.size(),0))
    ,omega_B_(std::vector<double> (x.size(),0))
    ,rho_(std::vector<std::vector<double> > (x.size(),std::vector<double>(numK,0)))
  {
    auto n=x.size();
    for (unsigned i=0; i<x_.size()-1; ++i)
      dx_[i]=x_[i+1]-x_[i];
    dx_[n-1]=dx_[n-2];
  }

  CortexState(const CortexState& other)=default;


  CortexState(CortexState&& other)=default;



  CortexState& operator =(const CortexState& other)=default;

  CortexState& operator =(CortexState&& other)=default;


};





class CortexSimulation
{
public:
  struct P
  {
    double dx;
    std::size_t nPoints_per_decade;
    double dtmax;
    double dtmin; double tequilibrio;
    double desiredError;
    double dtinf;
  };

  CortexSimulation()=default;

  CortexSimulation(const CortexSimulation& other)=default;

  CortexSimulation(CortexSimulation&& other)=default;

  CortexSimulation& operator =(CortexSimulation&& other)=default;

  CortexSimulation& operator =(const CortexSimulation& other)=default;


  std::ostream& write(std::ostream& s);

  std::ostream& write(std::ostream& s, const std::string& var, const std::string& par );
  std::ostream& writeHeaderDataFrame(std::ostream& os)const
  {
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<"psi_T...t.."<<t_[idt]<<"..x.."<<x_[idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<"psi_B...t.."<<t_[idt]<<"..x.."<<x_[idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<"omega_T...t.."<<t_[idt]<<"..x.."<<x_[idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<"omega_B...t.."<<t_[idt]<<"..x.."<<x_[idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            for (std::size_t type=0; type<rho_[idt][idx].size(); ++type)
              {
                os<<"rho...t.."<<t_[idt]<<"..x.."<<x_[idx]<<"..type.."<<type;
                if (!((type+1==rho_[idt][idx].size())&&(idx+1==x_.size())&&(idt+1==t_.size())))
                  os<<"\t";
              }
          }
      }
    return os;
  }

  std::ostream& writeRowDataFrame(std::ostream& os)const
  {
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<psi_T_[idt][idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<psi_B_[idt][idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<omega_T_[idt][idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            os<<omega_B_[idt][idx]<<"\t";
          }
      }
    for (std::size_t idt=0; idt<t_.size(); ++idt)
      {
        for (std::size_t idx=0; idx<x_.size(); ++idx)
          {
            double simVol_liter=h_*h_*(x_[idx+1]-x_[idx])*1e-6*1e3;
            for (std::size_t type=0; type<rho_[idt][idx].size(); ++type)
              {
                os<<rho_[idt][idx][type]/simVol_liter;
                if (!((type+1==rho_[idt][idx].size())&&(idx+1==x_.size())&&(idt+1==t_.size())))
                  os<<"\t";
              }
          }
      }

    return os;
  }

  void writeDataFrameHeader(std::ostream& os)const
  {
    os<<"x"<<"\t";
    os<<"t"<<"\t";
    os<<"type"<<"\t";
    os<<"psi_T"<<"\t";
    os<<"psi_B"<<"\t";
    os<<"omega_T"<<"\t";
    os<<"omega_B"<<"\t";
    os<<"rho"<<"\t";

  }
  void writeDataFrame(std::ostream& os,
                      std::size_t idt,
                      std::size_t idx,
                      std::size_t type) const
  {
    os<<x_[idx]<<"\t";
    os<<t_[idt]<<"\t";
    os<<type<<"\t";
    os<<psi_T_[idt][idx]<<"\t";
    os<<psi_B_[idt][idx]<<"\t";
    os<<omega_T_[idt][idx]<<"\t";
    os<<omega_B_[idt][idx]<<"\t";
    double simVol_liter=h_*h_*(x_[idx+1]-x_[idx])*1e-6*1e3;
    os<<rho_[idt][idx][type]/simVol_liter;
  }


  void writeDataFrame(std::ostream& os) const
  {
    writeDataFrameHeader(os);
    os<<"\n";

    for (std::size_t idt=0; idt<t_.size(); ++idt)
      for (std::size_t idx=0; idx<x_.size(); ++idx)
        for (std::size_t type=0; type<rho_[0][0].size(); ++type)
          {
            writeDataFrame(os,idt,idx,type);
            os<<"\n";
          }
    os<<"\n";

  }

  void writeRowDataFrame(std::ostream& os, std::string s) const
  {

    for (std::size_t idt=0; idt<t_.size(); ++idt)
      for (std::size_t idx=0; idx<x_.size(); ++idx)
        for (std::size_t type=0; type<rho_[0][0].size(); ++type)
          {
            os<<s;
            writeDataFrame(os,idt,idx,type);
            os<<"\n";
          }
  }

  bool isValid_=false;

  std::string id_;

  Parameters p_;

  double dt_;
  double h_;

  std::vector<double> x_;
  std::vector<double> dx_;

  std::vector<double> t_;
  std::vector<double> maxlogErrt_;

  std::vector<double> logError_;


  std::vector<std::vector<double>> psi_T_;

  std::vector<std::vector<double>> psi_B_;


  std::vector<std::vector<double>> omega_T_;

  std::vector<std::vector<double>> omega_B_;


  std::vector<std::vector<std::vector<double>>> rho_;




  CortexSimulation(const CortexState& c,unsigned numSamples):
    isValid_(false)
  ,id_()
  ,p_(),dt_(),h_(),
    x_(c.x_),dx_(c.dx_),
    t_(std::vector<double>(numSamples)),
    maxlogErrt_(std::vector<double>(numSamples)),
    psi_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.psi_T_.size()))),
    psi_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.psi_T_.size()))),
    omega_T_
    (  std::vector<std::vector<double>>(numSamples,std::vector<double>(c.omega_T_.size()))
                                       ),
    omega_B_
    (  std::vector<std::vector<double>>(numSamples,
                                        std::vector<double>(c.omega_B_.size()))
                                       ),
    rho_
    (std::vector<std::vector<std::vector<double>>>(numSamples,std::vector<std::vector<double>>
                                                   (c.rho_.size(),
                                                    std::vector<double>(
                                                      c.rho_.front().size()))))

  {
    psi_T_[0]=c.psi_T_;
    if (c.hasOmega_)
      omega_T_[0]=c.omega_T_;
    rho_[0]=c.rho_;

  }


  CortexSimulation(const std::string& id,unsigned numSamples,unsigned numNodes,unsigned numStates):
    isValid_(false)
  ,id_(id)
  ,p_(),dt_(),h_()
  ,x_(std::vector<double>(numNodes)),
    dx_(std::vector<double>(numNodes)),
    t_(std::vector<double>(numSamples)),
    maxlogErrt_(std::vector<double>(numSamples)),
    psi_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    psi_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    omega_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    omega_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    rho_(std::vector<std::vector<std::vector<double>>>(
                                                       numSamples,std::vector<std::vector<double>>
                                                       (numNodes,
                                                        std::vector<double>(
                                                          numStates)))){}






  void read(std::string &line, std::istream &s, std::ostream &logs);
};

class Experiment;


#endif // CORTEXSIMULATION

