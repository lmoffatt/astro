#include "CortexLikelihood.h"
#include "CommandManager.h"
#include "Splines.h"

/// x indica los intervalos de la simulacion
/// y indica el promedio dentro de estos intervalos
/// xint indica los intervalos interpolados
/// @return una matriz de una fila menos que xint con los promedios en los
/// intervalos xint a partir de x e y

std::vector<std::vector<double>> interpolateInjury(const std::vector<double>&x,
                                                   const std::vector<std::vector<double> >& y,
                                                   const std::vector<double>& xint,
                                                   double injLenth)
{
  std::vector<std::vector<double>> cumY(x.size(),std::vector<double>(y[0].size(),0));


  std::vector<std::vector<double>> o(xint.size());

  for (std::size_t i=0; i<x.size()-1; ++i)
    for (std::size_t j=0; j<y[0].size(); ++j)
      cumY[i+1][j]=cumY[i][j]+y[i][j];
  MIExpSpline sp(x,cumY);

  std::vector<double> cumrho0(y[0].size(),0);
  std::vector<double> cumrho;
  double currx;
  for (std::size_t i=0; i<xint.size(); ++i)
    {
      currx=xint[i]+injLenth;
      cumrho=sp.eval(currx);
      o[i]=(cumrho-cumrho0);
      cumrho0=cumrho;
    }
  return o;

}







const Parameters &CortexLikelihood::getPrior() const
{
  return prior_;
}




CortexLikelihood::CortexLikelihood(std::string id, const Experiment *e, const Parameters &prior, double dx, double dtmin, std::size_t nPoints_per_decade,double dtmax, double tequilibrio):
  prior_(prior)
,m_()
,e_(e)
,dx_(dx)
,dtmin_(dtmin)
,nPoints_per_decade_(nPoints_per_decade)
,dtmax_(dtmax)
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
  return s;
}

CortexLikelihood::~CortexLikelihood(){}

std::ostream &CortexLikelihood::write(std::ostream &s) const
{
  s<<"Likelihood \n";

  s<<"Experiment \t"<<e_->id()<<"\n";
  s<<"Prior  \t" <<getPrior().id();

  return s;
}

const Experiment *CortexLikelihood::getExperiment() const{return e_;}

const BaseModel *CortexLikelihood::getModel() const {return m_;}

std::vector<std::vector<double> > CortexLikelihood::getstate(const Experiment *e)
{
  std::vector<std::vector<double>> o;
  for (unsigned ie=0; ie<e->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e->getMeasure(ie);
      std::vector<double> rho_meas(m_->getNumberOfObservedStates(),0);
      if (cm->inj_Width()>0)
        {
          for (std::size_t i=0; i<m_->getNumberOfSimulatedStates(); ++i)
            o.push_back(rho_meas);
        }

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
      if (cm->inj_Width()>0)
        o.insert(o.end(),getModel()->getNumberOfSimulatedStates(),cm->inj_Area());
      o.insert(o.end(),cm->areaAstro().begin(),cm->areaAstro().end());
    }
  return o;
}

void CortexLikelihood::update(){
  m_=BaseModel::create(prior_);
    bin_dens_=getNBins(e_);
    nstate_=getstate(e_);
    ntot_=getNtot(nstate_);

}


std::vector<double> CortexLikelihood::getNtot(const std::vector<std::vector<double> > nstate)
{
  std::vector<double> o(nstate.size(),0);
  for (std::size_t i=0; i<o.size(); ++i)
    for (std::size_t j=0; j<nstate[i].size(); ++j)
      o[i]+=nstate[i][j];
  return o;
}

std::__cxx11::string CortexLikelihood::ClassName(){return "CortexLikelihood";}

std::__cxx11::string CortexLikelihood::myClass() const {return ClassName();}

CortexLikelihood::CortexLikelihood(){}

std::ostream &CortexLikelihood::writeBody(std::ostream &s) const
{
  writeField(s,"prior",prior_);
  writeField(s,"experimental_results",e_->id());
  writeField(s,"grid_legth",dx_);
  writeField(s,"min_sample_time",dtmin_);
  writeField(s,"prod_sample_time",nPoints_per_decade_);

  writeField(s,"max_sample_time",dtmax_);

  writeField(s,"tiempo_equilibrio",tequilibrio_);
  return s;
}

void CortexLikelihood::clear()
{
  prior_.clear();
  e_=nullptr;
  nstate_.clear();
  ntot_.clear();
}

bool CortexLikelihood::readBody(std::string &line, std::istream &s, std::ostream& logs)
{

  if (!readField(line,s,"prior",prior_,logs)) return false;
  std::string experimentName;
  if (!readField(line,s,"experimental_results",experimentName,logs)) return false;
  else
    e_=cm_->getExperiment(experimentName);
  if (!readField(line,s,"grid_length",dx_,logs)) return false;
  if (!readField(line,s,"min_sample_time",dtmin_,logs)) return false;
  if (!readField(line,s,"prod_sample_time",nPoints_per_decade_,logs)) return false;

  if (!readField(line,s,"max_sample_time",dtmax_,logs)) return false;
  if (!readField(line,s,"tiempo_equilibrio",tequilibrio_,logs)) return false;
  if (!readField(line,s,"num_Astrocitos_cada_estado",nstate_,logs)) return false;
  if (!readField(line,s,"num_total_Astrocitos",ntot_,logs)) return false;
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

  BaseModel* m=BaseModel::create(parameters);

  CortexSimulation s=m->run(*e_,dx_,dtmin_,nPoints_per_decade_,dtmax_,tequilibrio_);
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








CortexSimulation CortexPoisonLikelihood::simulate(const Parameters &parameters)const
{
  BaseModel * m=BaseModel::create(parameters);
  CortexSimulation s=m->run(*e_,dx_,dtmin_,nPoints_per_decade_,dtmax_,tequilibrio_);
  delete m;
  return s;
}



std::vector<std::vector<double> > CortexPoisonLikelihood::f(const Parameters &parameters) const
{
  std::vector<std::vector<double>> o(nstate_.size());

  double h=1e-5;

  BaseModel * m=BaseModel::create(parameters);
  CortexSimulation s=m->run(*e_,dx_,dtmin_,nPoints_per_decade_,dtmax_,tequilibrio_);
  if (!s.isValid_)
    return{};
  std::size_t ic=0;
  unsigned is=0;

  for (unsigned ie=0; ie<e_->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e_->getMeasure(ie);
      double t=cm->dia()*24*60*60;
      while (s.t_[is]<t
             &&is<s.t_.size())
        ++is;


      double currInjury=parameters.get("inj_width_"+cm->id());
      if (std::isnan(currInjury))
        currInjury=0;


      // voy a interpolar los resultados de la medicion a partir de la simulacion.
      // la simulacion pasa a ser independiente de la medicion

      // rho debe estar en celulas por litro
      auto rho=interpolateInjury(s.x_,s.rho_[is],cm->xpos()*1e-6,currInjury*1e-6);

      ///horrible hack for the lession:
      /// dedico una fila para la probabilidad de cada estado antes de la lesion

      if ( cm->inj_Width()>0)
        {
          double injVolume_liters=cm->inj_Area()*1e-12*h*1000;
          double simVol_liters=cm->inj_Width()*1e-6*cm->h()*cm->h()*1000;
          double f=injVolume_liters/simVol_liters;



          for (std::size_t istate=0; istate<rho[0].size(); ++istate)
            {
              o[ic]=std::vector<double>(5,rho[0][istate]/5.0*f);

              ++ic;
            }

        }
      for (std::size_t ix=0; ix<rho.size()-1;++ix)
        {
          /// en celulas por litro
          auto rho_sim=m_->getObservedNumberFromModel(rho[ix+1]);

          double measuredVolume_inLiters=cm->areaAstro()[ix]*1e-12*h*1000; // en um cuadrados  h indica el espesor del corte
          double simVol_liters=cm->h()*cm->h()*(cm->xpos()[ix+1]-cm->xpos()[ix])*1e-6*1000;
          double f=measuredVolume_inLiters/simVol_liters;

          o[ic]=rho_sim*f;
          ++ic;
        }

    }
  return o;
}



std::vector<std::vector<double> > CortexPoisonLikelihood::g(const Parameters& parameters,const CortexSimulation &s) const
{
  std::vector<std::vector<double>> o(nstate_.size());

  double h=1e-5;

  if (s.x_.empty())
    return{};
  std::size_t ic=0;
  unsigned is=0;

  for (unsigned ie=0; ie<e_->numMeasures(); ie++)
    {
      const CortexMeasure* cm=e_->getMeasure(ie);
      double t=cm->dia()*24*60*60;
      while (s.t_[is]<t
             &&is<s.t_.size())
        ++is;


      double currInjury=parameters.get("inj_width_"+cm->id());
      if (std::isnan(currInjury))
        currInjury=0;


      // voy a interpolar los resultados de la medicion a partir de la simulacion.
      // la simulacion pasa a ser independiente de la medicion

      // rho debe estar en celulas por litro
      auto rho=interpolateInjury(s.x_,s.rho_[is],cm->xpos()*1e-6,currInjury*1e-6);

      ///horrible hack for the lession:
      /// dedico una fila para la probabilidad de cada estado antes de la lesion

      if ( cm->inj_Width()>0)
        {
          double injVolume_liters=cm->inj_Area()*1e-12*h*1000;
          double simVol_liters=cm->inj_Width()*1e-6*cm->h()*cm->h()*1000;
          double f=injVolume_liters/simVol_liters;



          for (std::size_t istate=0; istate<rho[0].size(); ++istate)
            {
              o[ic]=std::vector<double>(5,rho[0][istate]/5.0*f);

              ++ic;
            }

        }
      for (std::size_t ix=0; ix<rho.size()-1;++ix)
        {
          /// en celulas por litro
          auto rho_sim=m_->getObservedNumberFromModel(rho[ix+1]);

          double measuredVolume_inLiters=cm->areaAstro()[ix]*1e-12*h*1000; // en um cuadrados  h indica el espesor del corte
          double simVol_liters=cm->h()*cm->h()*(cm->xpos()[ix+1]-cm->xpos()[ix])*1e-6*1000;
          double f=measuredVolume_inLiters/simVol_liters;

          o[ic]=rho_sim*f;
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

std::ostream &CortexMultinomialLikelihoodEvaluation::extract(std::ostream &s, const std::string &, const std::string &) const
{

  std::size_t numMeasures=CL_->getExperiment()->numMeasures();

  std::vector<double> x;
  for (std::size_t i=0; i<numMeasures; ++i)
    {
      if (CL_->getExperiment()->getMeasure(i)->inj_Width()>0)
        {
          double injW=this->p_.get("inj_width_"+CL_->getExperiment()->getMeasure(i)->id());
        x.insert(x.end(),this->CL_->getModel()->getNumberOfSimulatedStates(),-injW);
        }
      x.insert(x.end(),++CL_->getExperiment()->getMeasure(i)->xpos().begin(),
               CL_->getExperiment()->getMeasure(i)->xpos().end());
    }

  CortexSimulation simul=CL_->simulate(p_);
  std::vector<std::vector<double>> f=CL_->g(p_,simul);





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
      if (CL_->getExperiment()->getMeasure(ime)->inj_Width()>0)
        ie=i0+CL_->getExperiment()->getMeasure(ime)->meanAstro().size()
            +CL_->getModel()->getNumberOfSimulatedStates();
      else
        ie=i0+CL_->getExperiment()->getMeasure(ime)->meanAstro().size();

      std::string title=CL_->getExperiment()->getMeasure(ime)->id()+"_predicted_measured_Likelihdood";
      writeTable(s,title,xtitle,x,ytitles,ytable,i0,ie);
    }

  s<<"//------------------------Simulations------------------------//\n";

  x=simul.x_;

   ytitles.clear();
   ytitles.resize(3);



  ytitles[0].push_back("DAMP_tot");
  ytitles[0].push_back("DAMP_free");
  ytitles[0].push_back("DAMP_bound");

  ytitles[1].push_back("MED_tot");
  ytitles[1].push_back("MED_free");
  ytitles[1].push_back("MED_bound");

  ytitles[2].push_back("d_pred_Neur");

  for (std::size_t is=1; is<simul.rho_[0][0].size(); ++is)
    {
     ytitles[2].push_back("d_pred_A"+std::to_string(is-1));
    }



  std::size_t numSims=CL_->getExperiment()->numSimPoints();
  for (std::size_t isim=0; isim<numSims; ++isim)
    {
      ytable.clear();
      ytable.resize(3);
      std::vector<std::vector<double>> psi(simul.psi_B_[isim].size(),std::vector<double>(3));
      ytable[0]=psi;
      ytable[1]=psi;

      for (std::size_t ix=0; ix<simul.psi_B_[isim].size(); ++ix)
        {
          ytable[0][ix][0]=simul.psi_T_[isim][ix];
          ytable[0][ix][1]=(simul.psi_T_[isim][ix]-simul.psi_B_[isim][ix]);
          ytable[0][ix][2]=simul.psi_B_[isim][ix];

          ytable[1][ix][0]=simul.omega_T_[isim][ix];
          ytable[1][ix][1]=(simul.omega_T_[isim][ix]-simul.omega_B_[isim][ix]);
          ytable[1][ix][2]=simul.omega_B_[isim][ix];

        }


      ytable[2]=simul.rho_[isim];

      std::string title="simulations_at"+std::to_string(simul.t_[isim])+"_s";

      writeTable(s,title,xtitle,x,ytitles,ytable);

}
  simul.write(s);



  s<<"//---------------------Parameters------------------------//\n";
  s<<"model"<<this->p_.model()<<"\n";

  s<<"name \t fitmean \t dBfit \t prior mean \t dB prior\n";
  for (std::size_t i=0; i<p_.size(); ++i)
    {
      std::string name=p_.names()[i];
      s<<name<<"\t"<<p_.mean(name)<<"\t"<<p_.pStd(name)<<"\t"<<CL_->getPrior().mean(name)<<"\t"<<CL_->getPrior().pStd(name)<<"\n";
    }
return s;
}




std::ostream &CortexPoisonLikelihoodEvaluation::extract(std::ostream &s, const std::string &, const std::string &) const
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















