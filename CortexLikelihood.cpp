#include "CortexLikelihood.h"
#include "CommandManager.h"





const Parameters &CortexLikelihood::getPrior() const
{
  return prior_;
}




CortexLikelihood::CortexLikelihood(std::string id, const Experiment *e, const Parameters &prior, double dt, double tequilibrio):
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

const std::vector<std::vector<double> > &CortexLikelihood::n_obs() const
{
  return nstate_;
}




const std::vector<double> &CortexLikelihood::n_obs(unsigned i) const
{
  return nstate_[i];
}

const std::vector<double> &CortexLikelihood::bin_dens() const
{
  return bin_dens_;
}

const std::vector<double> &CortexLikelihood::ntot_obs() const
{
  return ntot_;
}

std::ostream &CortexLikelihood::put(std::ostream &s) const
{

}

std::ostream &CortexLikelihood::write(std::ostream &s) const
{
  s<<"Likelihood \n";

  s<<"Experiment \t"<<e_->id()<<"\n";
  s<<"Prior  \t" <<getPrior().id();


}

std::vector<std::vector<double> > CortexLikelihood::getstate(const Experiment *e)
{
  std::vector<std::vector<double>> o;
  for (unsigned ie=0; ie<e->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e->getMeasure(ie);
      std::vector<double> rho_meas;
      for (unsigned ix=0; ix<cm->meanAstro().size(); ++ix)
        {
          auto ob=cm->meanAstro(ix);
          rho_meas=m_->getObservedNumberFromData(ob);
          o.push_back(rho_meas);
        }
    }
  return o;
}


std::vector<double>  CortexLikelihood::getNBins(const Experiment *e)
{
  std::vector<double> o;
  for (unsigned ie=0; ie<e->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e->getMeasure(ie);
      o.insert(o.end(),cm->areaAstro().begin(),cm->areaAstro().end());
    }
  return o;
}


std::vector<double> CortexLikelihood::getNtot(const std::vector<std::vector<double> > nstate)
{
  std::vector<double> o(nstate.size(),0);
  for (std::size_t i=0; i<o.size(); ++i)
    for (std::size_t j=0; j<nstate[i].size(); ++j)
      o[i]+=nstate[i][j];
  return o;
}

bool CortexLikelihood::readBody(std::string &line, std::istream &s)
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


void CortexLikelihood::setPrior(const Parameters &parameters)
{
  prior_=parameters;
}




const ABC_Freq_obs& CortexLikelihood::getData() const
{
  return *this;
}



//9----------------------------------------------------------



std::vector<std::vector<double> > CortexMultinomialLikelihood::f(const Parameters &parameters) const
{
  std::vector<std::vector<double>> o(nstate_.size());

  m_->loadParameters(parameters);

  CortexSimulation s=m_->run(*e_,dt_,tequilibrio_);
  if (s.dx_.empty())
    return {};
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
      for (unsigned ix=0; ix<cm->meanAstro().size(); ++ix)
        {
          rho_sim=m_->getObservedProbFromModel(s.rho_[is][ix]);
          o[ic]=rho_sim;
          ++ic;
        }

    }
  return o;
}










std::vector<std::vector<double> > CortexPoisonLikelihood::f(const Parameters &parameters) const
{
  std::vector<std::vector<double>> o(nstate_.size());

  double h=1e-5;

  m_->loadParameters(parameters);

  CortexSimulation s=m_->run(*e_,dt_,tequilibrio_);
  if (s.x_.empty())
    return{};
  std::size_t ic=0;
  unsigned is=0;
  double width_lesion=0;
  for (unsigned ie=0; ie<e_->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e_->getMeasure(ie);
      double t=cm->dia()*24*60*60;
      while (s.t_[is]<t
             &&is<s.t_.size())
        ++is;
      std::vector<double> rho_sim;
      std::vector<double> rho_sim_0;
      std::vector<double> rho_sim_00;
      // esto es un hack desesperado que solo sirve para un experimento de dos dias nomas.
      // lo que hace es descartar todos los bines anteriores a la lesion, la cual es un parametro del modelo

      width_lesion+=parameters.get("inj_width_"+std::to_string(ie));
      double landa;

      std::size_t numExtraBins=std::ceil(width_lesion/(e_->xpos()[1]-e_->xpos()[0]));
      if (numExtraBins>0)
        {
          rho_sim_0=m_->getObservedNumberFromModel(s.rho_[is][numExtraBins-1])
              *(1.0/cm->h()/cm->h()/s.dx_[numExtraBins-1]);

           landa=(s.x_[numExtraBins]-width_lesion*1e-6)/s.dx_[numExtraBins-1];
        }
      for (unsigned ix=numExtraBins; ix<cm->meanAstro().size()+numExtraBins; ++ix)
        {
          /// en celulas por cuadrado de h (100 micrones) por dx
          rho_sim=m_->getObservedNumberFromModel(s.rho_[is][ix])*(1.0/cm->h()/cm->h()/s.dx_[ix]); // c->h() es una variable de la simulacion, creo que esta alpedo.
          if (!rho_sim_0.empty())
            {
              rho_sim_00=rho_sim;
              rho_sim=rho_sim*(1-landa)+rho_sim_0*landa;
              rho_sim_0=rho_sim_00;
            }

            double measuredArea=cm->areaAstro()[ix-numExtraBins]*1e-12*h; // en um cuadrados  h indica el espesor del corte

          o[ic]=rho_sim*measuredArea;
          ++ic;
        }

    }
  return o;
}





















//--------------------------------------------------------------------------

std::ostream &CortexMultinomialLikelihoodEvaluation::writeBody(std::ostream &s) const
{
  writeField(s,"Likelihood_Model",CL_->id());
  writeField(s,"Evaluated_Parameters",p_);
  writeField(s,"Expected_State_Probabilities",p_exp_);

  return s;

}

std::ostream &CortexMultinomialLikelihoodEvaluation::extract(std::ostream &s, const std::string &s1, const std::string &s2) const
{

  std::size_t numMeasures=CL_->getExperiment()->numMeasures();

  std::vector<double> x;
  for (std::size_t i=0; i<numMeasures; ++i)
    x.insert(x.begin(),++CL_->getExperiment()->getMeasure(i)->xpos().begin(),
             CL_->getExperiment()->getMeasure(i)->xpos().end());

  std::vector<std::vector<double>> f=CL_->f(p_);



  std::vector<double> ntot=CL_->getData().ntot_obs();

  std::vector<double> binDens=CL_->getData().bin_dens();



  std::vector<std::vector<double>> nobs=CL_->getData().n_obs();





  auto logLcell=CL_->logLikCells(f);
  auto logLSam=CL_->logLikSamples(f);

  if (CL_->myClass()==CortexMultinomialLikelihood::ClassName())
    {
      for (std::size_t i=0; i<f.size(); ++i)
        {
          for (std::size_t j=0; j< f[i].size(); ++j)
            f[i][j]*=ntot[i];
        }

    }
  std::vector<std::vector<double>> dobs=nobs;
  std::vector<std::vector<double>> dpred=f;

  for (std::size_t i=0; i<f.size(); ++i)
    {
      for (std::size_t j=0; j< f[i].size(); ++j)
        {
          dobs[i][j]/=binDens[i]/1e4;
          dpred[i][j]/=binDens[i]/1e4;
        }

    }





  std::string xtitle="distance_to_lession[um]";

  std::vector<std::vector<std::string>> ytitles(8);

  ytitles[0].push_back("num_meas");
  ytitles[1].push_back("den_meas");

  for (std::size_t is=0; is<f[0].size(); ++is)
    {
      ytitles[2].push_back("d_meas_"+std::to_string(is));
      ytitles[3].push_back("d_pred_"+std::to_string(is));

      ytitles[4].push_back("n_meas_"+std::to_string(is));
      ytitles[5].push_back("n_pred_"+std::to_string(is));
      ytitles[6].push_back("logL_"+std::to_string(is));
    }
  ytitles[7].push_back("logL_total");

  std::vector<std::vector<std::vector<double>>> ytable;
  ytable.push_back(std::vector<std::vector<double>>(1,ntot));
  ytable.push_back(std::vector<std::vector<double>>(1,binDens));
  ytable.push_back(dobs);
  ytable.push_back(dpred);
  ytable.push_back(nobs);
  ytable.push_back(f);
  ytable.push_back(logLcell);
  ytable.push_back(std::vector<std::vector<double>>(1,logLSam));
  std::size_t i0=0;
  std::size_t ie=0;
  for (std::size_t ime=0; ime<numMeasures; ++ime)
    {
      i0=ie;
      ie=i0+CL_->getExperiment()->getMeasure(ime)->meanAstro().size();
      std::string title=CL_->getExperiment()->getMeasure(ime)->id()+"_predicted_measured_Likelihdood";
      writeTable(s,title,xtitle,x,ytitles,ytable,i0,ie);
    }

  return s;
}

std::ostream &CortexPoisonLikelihoodEvaluation::extract(std::ostream &s, const std::string &s1, const std::string &s2) const
{

  std::size_t numMeasures=CL_->getExperiment()->numMeasures();

  std::vector<double> x;
  for (std::size_t i=0; i<numMeasures; ++i)
    x.insert(x.begin(),++CL_->getExperiment()->getMeasure(i)->xpos().begin(),
             CL_->getExperiment()->getMeasure(i)->xpos().end());

  std::vector<std::vector<double>> f=CL_->f(p_);



  std::vector<double> ntot=CL_->getData().ntot_obs();

  std::vector<double> binDens=CL_->getData().bin_dens();



  std::vector<std::vector<double>> nobs=CL_->getData().n_obs();





  auto logLcell=CL_->logLikCells(f);
  auto logLSam=CL_->logLikSamples(f);

  if (CL_->myClass()==CortexMultinomialLikelihood::ClassName())
    {
      for (std::size_t i=0; i<f.size(); ++i)
        {
          for (std::size_t j=0; j< f[i].size(); ++j)
            f[i][j]*=ntot[i];
        }

    }
  std::vector<std::vector<double>> dobs=nobs;
  std::vector<std::vector<double>> dpred=f;

  for (std::size_t i=0; i<f.size(); ++i)
    {
      for (std::size_t j=0; j< f[i].size(); ++j)
        {
          dobs[i][j]/=binDens[i]/1e4;
          dpred[i][j]/=binDens[i]/1e4;
        }

    }





  std::string xtitle="distance_to_lession[um]";

  std::vector<std::vector<std::string>> ytitles(7);

  ytitles[0].push_back("num_meas");
  for (std::size_t is=0; is<f[0].size(); ++is)
    {
      ytitles[1].push_back("d_meas_"+std::to_string(is));
      ytitles[2].push_back("d_pred_"+std::to_string(is));

      ytitles[3].push_back("n_meas_"+std::to_string(is));
      ytitles[4].push_back("n_pred_"+std::to_string(is));
      ytitles[5].push_back("logL_"+std::to_string(is));
    }
  ytitles[6].push_back("logL_total");

  std::vector<std::vector<std::vector<double>>> ytable;
  ytable.push_back(std::vector<std::vector<double>>(1,ntot));
  ytable.push_back(dobs);
  ytable.push_back(dpred);
  ytable.push_back(nobs);
  ytable.push_back(f);
  ytable.push_back(logLcell);
  ytable.push_back(std::vector<std::vector<double>>(1,logLSam));
  std::size_t i0=0;
  std::size_t ie=0;
  for (std::size_t ime=0; ime<numMeasures; ++ime)
    {
      i0=ie;
      ie=i0+CL_->getExperiment()->getMeasure(ime)->meanAstro().size();
      std::string title=CL_->getExperiment()->getMeasure(ime)->id()+"_predicted_measured_Likelihdood";
      writeTable(s,title,xtitle,x,ytitles,ytable,i0,ie);
    }

  return s;
}

