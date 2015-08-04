#ifndef CORTEXSIMULATION
#define CORTEXSIMULATION


#include <iomanip>
#include <vector>
#include <map>
#include <limits>
#include <string>
#include <fstream>


class CortexState
{
public:
  std::vector<double> x_;
  std::vector<double> dx_;

  double h_;

  double epsilon_;

  std::vector<double> psi_T_;

  std::vector<double> psi_B_;

  std::vector<double> omega_T_;

  std::vector<double> omega_B_;

  std::vector<std::vector<double> > rho_;




  CortexState(const std::vector<double>& x
              ,const std::vector<double>& dx
              ,double h
              ,double epsilon
              ,unsigned numK)
    :
      x_(x)
    ,dx_(dx)
    ,h_(h)
    ,epsilon_(epsilon)
    ,psi_T_(std::vector<double>(x.size(),0))
    ,omega_T_(std::vector<double> (x.size(),0))
    ,rho_(std::vector<std::vector<double> > (x.size(),std::vector<double>(numK,0)))
  {}



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

  std::ostream& write(std::ostream& s)const
  {
    s<<"parameters \n";
    for (auto& elem:m_)
      {
        s<<elem.first<<"\t"<<elem.second.first<<"\t"<<elem.second.second<<"\n";
      }
    return s;
  }

private:

  std::map<std::string, std::pair<double,std::string>> m_;

};


class CortexSimulation
{
public:
  CortexSimulation(){}

  std::ostream& write(std::ostream& s);

  std::ostream& write(std::ostream& s, const std::string& var, const std::string& par );



  std::string id_;

  Parameters p_;

  double dt_;

  std::vector<double> x_;
  std::vector<double> dx_;

  std::vector<double> t_;
  std::vector<double> sdt_;

  std::vector<std::vector<double>> psi_T_;

  std::vector<std::vector<double>> psi_B_;


  std::vector<std::vector<double>> omega_T_;

  std::vector<std::vector<double>> omega_B_;


  std::vector<std::vector<std::vector<double>>> rho_;

  CortexSimulation(const CortexState& c,unsigned numSamples):
    id_()
  ,p_(),dt_(),
    x_(c.x_),dx_(c.dx_),
    t_(std::vector<double>(numSamples)),
    sdt_(std::vector<double>(numSamples)),
    psi_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.psi_T_.size()))),
    psi_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.psi_T_.size()))),
    omega_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.omega_T_.size()))),
    omega_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(c.omega_T_.size()))),
    rho_(std::vector<std::vector<std::vector<double>>>(
                                                       numSamples,std::vector<std::vector<double>>
                                                       (c.rho_.size(),
                                                        std::vector<double>(
                                                          c.rho_.front().size()))))

  {
    psi_T_[0]=c.psi_T_;
    omega_T_[0]=c.omega_T_;
    rho_[0]=c.rho_;

  }


  CortexSimulation(const std::string& id,unsigned numSamples,unsigned numNodes,unsigned numStates):
    id_(id)
  ,p_(),dt_()
  ,x_(std::vector<double>(numNodes)),
    dx_(std::vector<double>(numNodes)),
    t_(std::vector<double>(numSamples)),
    sdt_(std::vector<double>(numSamples)),
    psi_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    psi_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    omega_T_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    omega_B_(std::vector<std::vector<double>>(numSamples,std::vector<double>(numNodes))),
    rho_(std::vector<std::vector<std::vector<double>>>(
                                                       numSamples,std::vector<std::vector<double>>
                                                       (numNodes,
                                                        std::vector<double>(
                                                          numStates)))){}




  void read(std::string &line, std::istream &s);
};



#endif // CORTEXSIMULATION

