#include "CortexLikelihood.h"
#include "CommandManager.h"






const Parameters &CortexModelLikelihood::getPrior() const
{
  return prior_;
}

std::vector<std::vector<double> > CortexModelLikelihood::p_exp(const Parameters &parameters) const
{
  std::vector<std::vector<double>> o(nstate_.size());

  m_->loadParameters(parameters);

  CortexSimulation s=m_->run(*e_,dt_,tequilibrio_);
  std::size_t ic=0;
  unsigned is=0;
  for (unsigned ie=0; ie<e_->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e_->getMeasure(ie);
      double t=cm->dia()*24*60*60;
      while (s.t_[is]<t
             &&is<s.t_.size())
        ++is;
      std::vector<double> rho_sim;
      for (unsigned ix=0; ix<cm->dx().size(); ++ix)
        {
          rho_sim=m_->getObservedProbFromModel(s.rho_[is][ix]);
          o[ic]=rho_sim;
          ++ic;
        }

    }
  return o;
}



CortexModelLikelihood::CortexModelLikelihood(std::string id, const Experiment *e, const Parameters &prior, double dt, double tequilibrio):
 prior_(prior)
,m_()
,e_(e)
,dt_(dt)
,tequilibrio_(tequilibrio)
,nstate_()
,ntot_()
{
  setId(id);
  update();

}

const std::vector<std::vector<double> > &CortexModelLikelihood::n_obs() const
{
 return nstate_;
}




const std::vector<double> &CortexModelLikelihood::n_obs(unsigned i) const
{
  return nstate_[i];
}

const std::vector<double> &CortexModelLikelihood::ntot_obs() const
{
  return ntot_;
}

std::ostream &CortexModelLikelihood::put(std::ostream &s) const
{

}

std::ostream &CortexModelLikelihood::write(std::ostream &s) const
{
  s<<"Likelihood \n";

  s<<"Experiment \t"<<e_->id()<<"\n";
  s<<"Prior  \t" <<getPrior().id();


}

std::vector<std::vector<double> > CortexModelLikelihood::getstate(const Experiment *e)
{
  std::vector<std::vector<double>> o;
  for (unsigned ie=0; ie<e->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e->getMeasure(ie);
      std::vector<double> rho_meas;
      for (unsigned ix=0; ix<cm->dx().size(); ++ix)
        {
          auto ob=cm->meanAstro(ix);
          rho_meas=m_->getObservedNumberFromData(ob);
          o.push_back(rho_meas);
        }
    }
  return o;
}

std::vector<double> CortexModelLikelihood::getNtot(const std::vector<std::vector<double> > nstate)
{
  std::vector<double> o(nstate.size(),0);
  for (std::size_t i=0; i<o.size(); ++i)
    for (std::size_t j=0; j<nstate[i].size(); ++j)
      o[i]+=nstate[i][j];
  return o;
}

bool CortexModelLikelihood::readBody(std::string &line, std::istream &s)
{

  if (!readField(line,s,"prior",prior_)) return false;
  std::string experimentName;
  if (!readField(line,s,"experimental_results",experimentName)) return false;
  else
    e_=cm_->getExperiment(experimentName);
  if (!readField(line,s,"sample_time",dt_)) return false;
  if (!readField(line,s,"tiempo_equilibrio",tequilibrio_)) return false;
  if (!readField(line,s,"num_Astrocitos_cada_estado",nstate_)) return false;
  if (!readField(line,s,"num_total_Astrocitos",ntot_)) return false;
  return true;
}


void CortexModelLikelihood::setPrior(const Parameters &parameters)
{
  prior_=parameters;
}

const ABC_Freq_obs& CortexModelLikelihood::getData() const
{
  return *this;
}

std::ostream &CortexLikelihoodEvaluation::writeBody(std::ostream &s) const
{
  writeField(s,"Likelihood_Model",CL_->id());
  writeField(s,"Evaluated_Parameters",p_);
  writeField(s,"Expected_State_Probabilities",p_exp_);

  return s;

}

std::ostream &CortexLikelihoodEvaluation::extract(std::ostream &s, const std::string &s1, const std::string &s2) const
{

  std::size_t numMeasures=CL_->getExperiment()->numMeasures();

  std::vector<double> x;
  for (std::size_t i=0; i<numMeasures; ++i)
    x.insert(x.begin(),CL_->getExperiment()->getMeasure(i)->dx().begin(),
             CL_->getExperiment()->getMeasure(i)->dx().end());

  std::vector<std::vector<double>> pexp=CL_->p_exp(p_);


  std::vector<double> ntot=CL_->getData().ntot_obs();

  std::vector<std::vector<double>> pobs=CL_->getData().n_obs();

  for (std::size_t i=0; i<pobs.size(); ++i)
    {
      if (ntot[i]>0)
      for (std::size_t j=0; j< pobs[i].size(); ++j)
        pobs[i][j]/=ntot[i];
    }

  auto logLcell=CL_->logLikCells(pexp);
  auto logLSam=CL_->logLikSamples(pexp);



  std::string xtitle="distance_to_lession[um]";

  std::vector<std::vector<std::string>> ytitles(5);

  ytitles[0].push_back("num_meas");
  for (std::size_t is=0; is<pexp[0].size(); ++is)
    {
      ytitles[1].push_back("p_meas_"+std::to_string(is));
      ytitles[2].push_back("p_pred_"+std::to_string(is));
      ytitles[3].push_back("logL_"+std::to_string(is));
    }
  ytitles[4].push_back("logL_total");

  std::vector<std::vector<std::vector<double>>> ytable;
  ytable.push_back(std::vector<std::vector<double>>(1,ntot));
  ytable.push_back(pobs);
  ytable.push_back(pexp);
  ytable.push_back(logLcell);
  ytable.push_back(std::vector<std::vector<double>>(1,logLSam));
  std::size_t i0=0;
  std::size_t ie=0;
  for (std::size_t ime=0; ime<numMeasures; ++ime)
    {
      i0=ie;
      ie=i0+CL_->getExperiment()->getMeasure(ime)->dx().size();
      std::string title=CL_->getExperiment()->getMeasure(ime)->id()+"_predicted_measured_Likelihdood";
      writeTable(s,title,xtitle,x,ytitles,ytable,i0,ie);
    }

  return s;
}
