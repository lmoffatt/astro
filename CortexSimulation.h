#ifndef CORTEXSIMULATION
#define CORTEXSIMULATION


#include <iomanip>
#include <vector>
#include <map>
#include <limits>
#include <string>
#include <fstream>

#include "Parameters.h"


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

  CortexState(const std::vector<double>& x
              ,double h
              ,double epsilon
              ,unsigned numK)
    :
      x_(x)
    ,dx_(x.size())
    ,h_(h)
    ,epsilon_(epsilon)
    ,psi_T_(std::vector<double>(x.size(),0))
    ,omega_T_(std::vector<double> (x.size(),0))
    ,rho_(std::vector<std::vector<double> > (x.size(),std::vector<double>(numK,0)))
  {
    auto n=x.size();
    for (unsigned i=0; i<x_.size()-1; ++i)
      dx_[i]=x_[i+1]-x_[i];
    dx_[n-1]=dx_[n-2];
  }



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

class Experiment;

class Likelihood
{
public:
  Likelihood(const Experiment* e, const CortexSimulation* s):
    e_(e), s_(s){}


  double total();
  std::vector<double> partial();

  private:
   const Experiment* e_;
   const CortexSimulation* s_;

};

#endif // CORTEXSIMULATION

