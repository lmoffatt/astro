#include "Models.h"
#include "CortexSimulation.h"
#include "CommandManager.h"
#include "LevenbergMarquardt.h"
#include "CortexLikelihood.h"
//#include "BayesIteration.h"
#include <sstream>
#include <iostream>


std::vector<double> label_to_sequence(std::string s)
{


  std::stringstream ss(s);

  std::vector<double> o;
  double val;
  double dx;
  char ch;
  while (ss>>val)
    {
      o.push_back(val);
      ss>>ch;
      if ((ch>='0')&&(ch<='9'))
        {
          ss.putback(ch);
        }
      else if (ch==':')
        {
          ss>>dx>>ch>>val;
          while (o.back()<val)
            o.push_back(o.back()+dx);
        }

    }

  return o;


}



int Script::run(char* filename)
{
  std::ifstream f(filename);
  while (f)
    {
      std::string line;
      safeGetline(f,line);
      removeComments(line);
      cm_->execute(line);
    }
  return 0;

}

int Script::runDefine(const std::string &filename, const std::vector<std::string> &label, const std::vector<std::string> &valueInplace)
{
  std::ifstream f(filename);
  while (f)
    {
      std::string line;
      safeGetline(f,line);
      removeComments(line);

      replaceLabel(line,label,valueInplace);
      cm_->execute(line);
    }
  return 0;

}



CommandManager::CommandManager()
{
  cmd_["read"]=new readCommand(this);
  cmd_["align"]=new AlignCommand(this);
  cmd_["write"]=new writeCommand(this);
  cmd_["merge"]=new MergeCommand(this);
  cmd_["distances"]=new DistancesCommand(this);
  cmd_["histogram"]=new HistogramCommand(this);
  cmd_["simulate"]=new SimulateCommand(this);
  cmd_["experiment"]=new ExperimentCommand(this);
  cmd_["likelihood"]=new LikelihoodCommand(this);
  cmd_["optimize"]=new OptimizeCommand(this);
  cmd_["evidence"]=new EvidenceCommand(this);
  cmd_["tempering"]=new TemperingCommand(this);




}

void CommandManager::execute(std::string line)
{

  std::stringstream ss(line);
  std::string command;
  ss>>command;
  CommandBase* c=Command(command);
  if (c!=0)
    c->run(line);

}

CommandBase *CommandManager::Command(std::string commandName)
{
  std::map<std::string,CommandBase*>::iterator it=cmd_.find(commandName);
  if (it!=cmd_.end())
    return it->second;
  else
    return 0;


}

TissueSection *CommandManager::getSection(std::string idSection)
{
  auto it=sections.find(idSection);
  if(it!=sections.end())
    return it->second;
  else
    return nullptr;
}

void CommandManager::push_back(BaseModel *model)
{
  models[model->id()]=model;
}

void CommandManager::push_back(CortexSimulation *simulation)
{
  simulations[simulation->id_]=simulation;
}

void CommandManager::push_back(CortexMeasure *measure)
{
  measures[measure->id()]=measure;
}

void CommandManager::push_back(Experiment *experiment)
{
  experiments[experiment->id()]=experiment;
}

//void CommandManager::push_back(mcmc *mcmc)
//{
//  mcmcs[mcmc->id()]=mcmc;
//}

void CommandManager::push_back(CortexLikelihood *likelihood)
{
  likelihood->setCommandManager(this);
  likelihoods[likelihood->id()]=likelihood;
}



CommandManager::~CommandManager()
{
  for (auto elem:cmd_)
    delete elem.second;
  for (auto elem: sections)
    delete elem.second;
  for (auto elem: measures)
    delete elem.second;
  for (auto elem: models)
    delete elem.second;
  for (auto elem: simulations)
    delete elem.second;
  for (auto elem: experiments)
    delete elem.second;
  for (auto elem: mcmcs)
    delete elem.second;
  for (auto elem: likelihoods)
    delete elem.second;


}

BaseModel *CommandManager::getModel(std::string idModel)
{
  auto it=models.find(idModel);
  if(it!=models.end())
    return it->second;
  else
    return nullptr;
}
CortexSimulation *CommandManager::getSimulation(std::string idSimulation)
{
  auto it=simulations.find(idSimulation);
  if(it!=simulations.end())
    return it->second;
  else
    return nullptr;
}


CortexMeasure *CommandManager::getMeasure(std::string id)
{
  auto it=measures.find(id);
  if(it!=measures.end())
    return it->second;
  else
    return nullptr;
}

Experiment *CommandManager::getExperiment(const std::string& id)
{
  auto it=experiments.find(id);
  if(it!=experiments.end())
    return it->second;
  else
    return nullptr;
}


mcmc *CommandManager::getMcmc(const std::string& id)
{
  auto it=mcmcs.find(id);
  if(it!=mcmcs.end())
    return it->second;
  else
    return nullptr;
}




CortexLikelihood *CommandManager::getLikelihood(const std::string &id)
{
  auto it=likelihoods.find(id);
  if(it!=likelihoods.end())
    return it->second;
  else
    return nullptr;

}




#include <sstream>



void AlignCommand::run(const std::string& line)
{

  std::string alignName, dataName;
  std::stringstream ss(line);
  ss>>alignName>>dataName;

  TissueSection* s=cm_->getSection(dataName);
  if (s!=nullptr)
    {
      s->align();
    }




}


void MergeCommand::run(const std::string& line)
{
  std::string cName, dataName;
  std::stringstream ss(line);
  ss>>cName>>dataName;

  TissueSection* s=cm_->getSection(dataName);
  if (s!=nullptr)
    {
      s->merge();
    }



}








void DistancesCommand::run(const std::string& line)
{
  std::string cName, dataName;
  std::stringstream ss(line);
  ss>>cName>>dataName;

  TissueSection* s=cm_->getSection(dataName);
  if (s!=nullptr)
    {
      s->distances();
    }



}




void HistogramCommand::run(const std::string& line)
{

  //histogram 3dpl d_les  0:+100:3000
  // histogram 3dpl d_les 1E7 1 0  0 [0 25 50 100 200 500:500:4000]   3dpl 7dpl

  std::size_t maxnumpoints,initseed;
  double minDistance_Tissue, minDistance_Vaso;

  std::string cName, dataName,kind;
  std::stringstream ss(line);
  ss>>cName>>dataName>>kind>>maxnumpoints>>initseed
      >>minDistance_Tissue>>minDistance_Vaso;

  std::string ch;
  while ((ss>>ch)&&(ch!="[")){}
  std::string interval;
  while ((ss>>ch)&&(ch!="]"))
    interval+=ch+" ";

  TissueSection* s=cm_->getSection(dataName);

  std::vector<double> limits=label_to_sequence(interval);

  std::mt19937_64 mt;
  if (initseed==0)
    {
      std::random_device rd;
      auto seed=rd();
      mt.seed(seed);
      std::cout<<"experiment uses seed="<<seed<<std::endl;
    }
  else if (initseed==1)
    {
      mt.seed();
      std::cout<<"experiment uses default seed"<<std::endl;
    }
  else
    {
      mt.seed(initseed);
      std::cout<<"experiment uses provided seed="<<initseed<<std::endl;

    }


  if (s!=nullptr)
    {
      s->measure(cm_,mt,limits,minDistance_Tissue, minDistance_Vaso,maxnumpoints);

    }


}


void ExperimentCommand::run(const std::string& line)

{
  // experiment experiment1 1E7 seed  0  0 [0 25 50 100 250 500:500:4000]   3dpl 7dpl

  std::size_t maxnumpoints, initseed;
  double minDistance_Tissue, minDistance_Vaso;
  std::string cName,expName, currsection;
  std::vector<std::string> sections;
  std::stringstream ss(line);
  //  experiment experiment_0 10000000 1 0  0  [0 25 50 100 250 500:500:4000]  3dplCL 3dpl 7dpl2


  ss>>cName>>expName>>maxnumpoints>>initseed>>minDistance_Tissue>>minDistance_Vaso;

  std::string ch;
  while ((ss>>ch)&&(ch!="[")){}
  std::string interval;
  while ((ss>>ch)&&(ch!="]"))
    interval+=ch+" ";

  while(ss>>currsection)
    sections.push_back(currsection);

  if (initseed==0)
    {
      std::random_device dv;
      initseed=dv();

    }


  std::vector<CortexMeasure> vCM;
  for (auto s:sections)
    {
      cm_->execute("read "+s);
      //cm_->execute("align "+s);
      //cm_->execute("merge "+s);
      cm_->execute("distances "+s);
      cm_->execute("write "+s);

      cm_->execute("histogram "+s+ " d_les "
                   +std::to_string(maxnumpoints)+"  "
                   +std::to_string(initseed)+"  "
                   +std::to_string(minDistance_Tissue)+"  "
                   +std::to_string(minDistance_Vaso)+"  [ "
                   +interval +" ]");
      vCM.push_back(*cm_->getMeasure(s));

    }

  Experiment* e=new Experiment(expName,vCM,{});
  cm_->push_back(e);

  e->store(expName+"_"+std::to_string(initseed));

}


void SimulateCommand::run(const std::string& line)
{

  //simulate opt experiment1 parameters_10 parameters_10 10 100

  std::string simulateS, simName, experimentName, priorName, paramName;
  double dtmin,dtmax, dx, tequilibrio=100000, maxduration;
  double factor=0;
  std::mt19937_64::result_type initseed=0;
  std::size_t nPoints_per_decade,niter,nseeds=0;
  std::stringstream ss(line);

  ss>>simulateS>>simName>>experimentName>>priorName>>paramName>>dx>>dtmin>>nPoints_per_decade>>dtmax>>niter>>maxduration>>factor>>nseeds>>initseed;

  Experiment *e=new Experiment;
  e->load(experimentName);
  if (e->numMeasures()==0)
    {
      std::cerr<<"Experiment "<<experimentName<<" not found\n";
      return;
    }

  std::string filename=priorName;
  std::fstream f;

  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          std::cerr<<"Parameters file "<<filename<<" or "<<filenaExt<<" not found"<<std::endl;
          f.close();
          return;
        }
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters prior;

  if (!prior.read(line2,f))
    {
      std::cerr<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();

  filename=paramName;
  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  if (!f)
    {
      std::cerr<<"Parameters file "<<filename<<" not found"<<std::endl;
      f.close();
      return;
    }


  safeGetline(f,line2);
  Parameters p;

  if (!p.read(line2,f))
    {
      std::cerr<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();

  BaseModel*m=BaseModel::create(prior);
  if (m!=nullptr)
    {
      cm_->push_back(m);

      CortexPoisonLikelihood  CL(simName+"_L",e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio);


      CortexMultinomialLikelihoodEvaluation CE(CL,p);
      std::ofstream fo;
      std::string fnameout=simName+"_lik.txt";
      fo.open(fnameout.c_str());
      CE.extract(fo);
      fo.close();



    }




}





void readCommand::run(const std::string& rline)
{
  std::stringstream ss(rline);
  std::string commandName;
  std::string filename;
  ss>>commandName>>filename;

  std::ifstream f(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          f.close();
          std::string below="../"+filenaExt;
          f.open(below);
        }
    }
  std::string line;
  safeGetline(f,line);
  double dia=0;
  std::size_t rata=0;
  if (line.find("experiment")!=line.npos)
    {
      safeGetline(f,line);
      while (f && line.empty())
        safeGetline(f,line);

      if (line.find("dia")!=line.npos)
        {
          std::stringstream ss(line);
          std::string diastr;

          ss>>diastr>>dia;
          safeGetline(f,line);

        }


      auto  s=new TissueSection(filename,dia);

      while (f)
        {

          if (line.find("rata")!=line.npos)
            {
              std::stringstream ss(line);
              std::string ratastr;

              ss>>ratastr>>rata;
              TissuePhoto foto;

              foto.read(line,f);
              s->fotos[foto.rata]=foto;

            }
          else if (line.find("foto")!=line.npos)
            {
              TissuePhoto foto;

              foto.read(line,f);
              s->fotos[foto.num]=foto;
            }
          else
            safeGetline(f,line);

        }
      cmd_->push_back(s);
    }
  else  if (line.find("simulation")!=line.npos)
    while (f)
      {
        auto sim=new CortexSimulation;
        sim->read(line,f);
        if (sim->t_.size()>0)
          {
            cmd_->push_back(sim);
          }

      }

}








void writeCommand::run(const std::string& line)
{
  std::string writeName, dataName;
  std::stringstream ss(line);
  ss>>writeName>>dataName;

  TissueSection* s=cmd_->getSection(dataName);
  if (s!=nullptr)
    {
      std::string filename=dataName+"_astro.txt";
      std::ofstream f;
      f.open(filename.c_str(),std::ofstream::out);

      for (auto fotop:s->fotos)
        fotop.second.write(f);

      f.close();

    }


  CortexMeasure* m=cmd_->getMeasure(dataName);
  if (m!=nullptr)
    {
      std::string filename=dataName+"_histo.txt";
      std::ofstream f;
      f.open(filename.c_str(),std::ofstream::out);

      m->print(f);
      f.close();

    }

  CortexSimulation* sim=cmd_->getSimulation(dataName);
  if (sim!=nullptr)
    {
      std::string var;
      ss>>var;
      std::string par;
      std::getline(ss,par);
      if (var.empty())
        {
          std::string newName;
          std::string filename=dataName+"_sim.txt";
          std::ifstream fi;
          fi.open(filename.c_str(),std::ios_base::in);
          if (fi)
            {
              unsigned num=0;
              newName=dataName+"_"+std::to_string(num);
              filename=newName+"_sim.txt";
              fi.close();
              fi.open(filename);
              while (fi)
                {
                  ++num;
                  newName=dataName+"_"+std::to_string(num);
                  filename=newName+"_sim.txt";
                  fi.close();
                  fi.open(filename);

                }
            }
          std::ofstream f;
          f.open(filename.c_str(),std::ofstream::out);
          sim->id_=newName;
          sim->write(f);
          f.close();

        }
      else
        {
          dataName=sim->id_;
          std::string filename=dataName+"_"+var+".txt";
          std::ifstream fi;
          fi.open(filename.c_str(),std::ios_base::in);
          if (fi)
            {
              unsigned num=0;
              filename=dataName+"_"+std::to_string(num)+"_"+var+".txt";
              fi.close();
              fi.open(filename);
              while (fi)
                {
                  ++num;
                  filename=dataName+"_"+std::to_string(num)+"_"+var+".txt";
                  fi.close();
                  fi.open(filename);

                }
            }
          std::ofstream f;
          f.open(filename.c_str(),std::ofstream::out);

          sim->write(f,var,par);
          f.close();

        }

    }



  //  CortexModelLikelihood* lik=cmd_->getLikelihood(dataName);
  //  if (lik!=nullptr)
  //    {
  //      std::string filename=lik+"_lik.txt";
  //      std::ofstream f;
  //      f.open(filename.c_str(),std::ofstream::out);


  //    }


}







void LikelihoodCommand::run(const std::string& line)
{

  //likelihood lik experiment1 parameters_10 parameters_10 0.5 50 1000

  std::string cName, lName, eName, priorName, paramName;
  std::size_t nPoints_per_decade;
  double dtmin,dtmax,dx, tequilibrio;
  std::stringstream ss(line);
  ss>>cName>>lName>>eName>>priorName>>paramName>>dtmin>>nPoints_per_decade>>dtmax>>dx>>tequilibrio;

  Experiment* e=cm_->getExperiment(eName);
  if (e==nullptr)
    {
      std::cerr<<"Experiment "<<eName<<" not found\n";
      return;
    }

  std::string filenamePrior=priorName;
  std::ifstream fp;

  fp.open(filenamePrior.c_str(),std::ios_base::in);
  if (!fp)
    {
      std::string filenaExt=filenamePrior+".txt";
      fp.open(filenaExt.c_str(),std::ios_base::in);
    }
  std::string line2p;
  safeGetline(fp,line2p);
  Parameters prior;

  prior.read(line2p,fp);

  fp.close();

  std::string filename=paramName;
  std::fstream f;

  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters p;

  p.read(line2,f);

  f.close();


  CortexLikelihood* CL= new CortexPoisonLikelihood(lName,e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio);
  cm_->push_back(CL);

  CortexMultinomialLikelihoodEvaluation CE(*CL,p);
  std::ofstream fo;
  std::string fnameout=lName+"_lik.txt";
  fo.open(fnameout.c_str(),std::ofstream::out);
  CE.extract(fo);
  fo.close();

}






std::string LikelihoodCommand::id() const
{
  return "likelihood";
}


void OptimizeCommand::run(const std::string& line)
{

  //optimize opt experiment1 parameters_10 parameters_10 10 100

  std::string optimizeS, optName, experimentName, priorName, paramName;
  double dtmin,dtmax, dx, tequilibrio=100000, maxduration;
  double factor=0;
  std::mt19937_64::result_type initseed=0;
  std::size_t nPoints_per_decade,niter,nseeds=0;
  std::stringstream ss(line);

  ss>>optimizeS>>optName>>experimentName>>priorName>>paramName>>dx>>dtmin>>nPoints_per_decade>>dtmax>>niter>>maxduration>>factor>>nseeds>>initseed;

  Experiment *e=new Experiment;
  e->load(experimentName);
  if (e->numMeasures()==0)
    {
      std::cerr<<"Experiment "<<experimentName<<" not found\n";
      return;
    }

  std::string filename=priorName;
  std::fstream f;

  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          std::cerr<<"Parameters file "<<filename<<" or "<<filenaExt<<" not found"<<std::endl;
          f.close();
          return;
        }
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters prior;

  if (!prior.read(line2,f))
    {
      std::cerr<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();

  filename=paramName;
  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  if (!f)
    {
      std::cerr<<"Parameters file "<<filename<<" not found"<<std::endl;
      f.close();
      return;
    }


  safeGetline(f,line2);
  Parameters p;

  if (!p.read(line2,f))
    {
      std::cerr<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();

  BaseModel*m=BaseModel::create(prior);
  if (m!=nullptr)
    {
      cm_->push_back(m);

      CortexPoisonLikelihood  CL(optName+"_lik",e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio);



      LevenbergMarquardtDistribution LM(&CL,p,niter,maxduration,optName);


      LM.optimize(optName,factor,nseeds,initseed);

    }




}



void EvidenceCommand::run(const std::__cxx11::string& line)
{

  // evidence evi experiment1 paramters_10
  //optimize opt experiment1 parameters_10 parameters_10 10 100
  typedef Landa AP;

  std::string evidenceS, eviName, experimentName, priorName;
  double dtmin,dtmax, dx, tequilibrio=100000, maxduration;
  double landa0,v;
  std::size_t nmaxloop;

  std::mt19937_64::result_type initseed=0;
  std::size_t nPoints_per_decade,niter;
  std::stringstream ss(line);
  M_Matrix<double> betas;
  M_Matrix<AP> aps;
  std::vector<std::vector<double>> apsPar;

  M_Matrix<std::size_t> samples, nskip;

  ss>>evidenceS>>eviName>>experimentName>>priorName>>dx>>dtmin>>nPoints_per_decade>>dtmax>>niter>>maxduration>>landa0>>v>>nmaxloop>>initseed>>aps>>apsPar>>betas>>samples>>nskip;


  std::cout<<"evidenceS: "<<evidenceS;
  std::cout<<" eviName: "<<eviName;
  std::cout<<" experimentName:"<<experimentName;
  std::cout<<" priorName: "<<priorName;
  std::cout<<" dx: "<<dx;
  std::cout<<" dtmin: "<<dtmin;
  std::cout<<" nPoints_per_decade: "<<nPoints_per_decade;
  std::cout<<" dtmax: "<<dtmax;
  std::cout<<" niter: "<<niter;
  std::cout<<" maxduration "<<maxduration;
  std::cout<<" landa0 "<<landa0;
  std::cout<<" v "<<v;
  std::cout<<" nmaxloop "<<nmaxloop;

  std::cout<<" initseed "<<initseed;
  std::cout<<" betas "<<betas;
  std::cout<<AP::ClassName()<<" "<<aps;
  for (std::size_t i=0; i<apsPar.size(); ++i)
    {
      std::cout<<AP::ParName(i)<<" "<<apsPar;
      for (std::size_t j=0; j<apsPar[i].size(); ++j)
        std::cout<<apsPar[i][j]<<" ";
    }

  std::cout<<" samples "<<samples;
  std::cout<<" nskip "<<nskip;

  Adaptive_discrete<AP> landa(aps,apsPar);

  std::cout<<landa;

  std::cerr<<line;

  Experiment *e=new Experiment;
  e->load(experimentName);
  if (e->numMeasures()==0)
    {
      std::cerr<<"Experiment "<<experimentName<<" not found\n";
      return;
    }

  std::string filename=priorName;
  std::fstream f;

  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          std::cerr<<"Parameters file "<<filename<<" or "<<filenaExt<<" not found"<<std::endl;
          f.close();
          return;
        }
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters prior;

  if (!prior.read(line2,f))
    {
      std::cerr<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();


  BaseModel*m=BaseModel::create(prior);
  if (m!=nullptr)
    {
      cm_->push_back(m);

      auto CL=new CortexPoisonLikelihood(eviName+"_lik",e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio);

      MyModel<MyData> m(CL);
      MyData d(CL);
      Metropolis_Hastings_mcmc<
          MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa> mcmc;
      LevenbergMarquardt_step<MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa> LMLik;
      Poisson_DLikelihood<MyData,MyModel> DLik;
      TI ti;
      std::mt19937_64 mt;
      std::random_device rd;

      if (initseed==0)
        {
          std::mt19937_64::result_type seed=rd();

          mt.seed(seed);
          eviName+=time_now()+"_"+std::to_string(seed);
          std::cerr<<"\n random seed =\n"<<seed<<"\n";
        }
      else
        {
          std::mt19937_64::result_type seed=initseed;

          mt.seed(seed);
          eviName+=time_now()+"_"+std::to_string(seed);

          std::cerr<<"\n provided seed =\n"<<seed<<"\n";

        }

      std::vector<std::tuple<double, std::size_t,std::size_t>>
          beta=getBeta(betas,samples,nskip);

      std::string eviNameLog=eviName+"_log.txt";
      std::ofstream flog;

      flog.open(eviNameLog.c_str(), std::ofstream::out | std::ofstream::app);
      flog<<line<<"\n";
      flog<<"evidenceS: "<<evidenceS<<"\n";
      flog<<" eviName: "<<eviName<<"\n";
      flog<<" experimentName:"<<experimentName<<"\n";
      flog<<" priorName: "<<priorName<<"\n";
      flog<<" dx: "<<dx<<"\n";
      flog<<" dtmin: "<<dtmin<<"\n";
      flog<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
      flog<<" dtmax: "<<dtmax<<"\n";
      flog<<" niter: "<<niter<<"\n";
      flog<<" maxduration "<<maxduration<<"\n";
      flog<<" landa0 "<<landa0<<"\n";
      flog<<" v "<<v<<"\n";
      flog<<" nmaxloop "<<nmaxloop<<"\n";

      flog<<" initseed "<<initseed<<"\n";
      flog<<" betas "<<betas<<"\n";
      flog<<AP::ClassName()<<" "<<aps<<"\n";
      flog<<" samples "<<samples<<"\n";
      flog<<" nskip "<<nskip<<"\n";
      flog.flush();


      std::chrono::steady_clock::time_point startTime=std::chrono::steady_clock::now();
      double timeOpt=0;


      typename TI::myEvidence * ev= ti.run(mcmc,LMLik,DLik,m,d,landa,beta,mt,flog,startTime,timeOpt);
      std::cout<<*ev;
      flog<<*ev;
      flog.close();
      std::ofstream fout;
      fout.open(eviName.c_str(), std::ofstream::out | std::ofstream::app);
      fout<<line<<"\n";
      fout<<"evidenceS: "<<evidenceS<<"\n";
      fout<<" eviName: "<<eviName<<"\n";
      fout<<" experimentName:"<<experimentName<<"\n";
      fout<<" priorName: "<<priorName<<"\n";
      fout<<" dx: "<<dx<<"\n";
      fout<<" dtmin: "<<dtmin<<"\n";
      fout<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
      fout<<" dtmax: "<<dtmax<<"\n";
      fout<<" niter: "<<niter<<"\n";
      fout<<" maxduration "<<maxduration<<"\n";
      fout<<" landa0 "<<landa0<<"\n";
      fout<<" v "<<v<<"\n";
      fout<<" nmaxloop "<<nmaxloop<<"\n";

      fout<<" initseed "<<initseed<<"\n";
      fout<<" betas "<<betas<<"\n";
      fout<<" samples "<<samples<<"\n";
      fout<<" nskip "<<nskip<<"\n";

      fout<<*ev<<"\n";
      fout.close();
    }




}




void TemperingCommand::run(const std::__cxx11::string& line)
{

  // evidence evi experiment1 paramters_10
  //optimize opt experiment1 parameters_10 parameters_10 10 100
  typedef Landa AP;

  std::string temperateS, eviName, experimentName, priorName;
  double dtmin,dtmax, dx, tequilibrio=100000, maxduration;
  double landa0,v;
  double pTjump;
  std::size_t nmaxloop;

  std::mt19937_64::result_type initseed=0;
  std::size_t nPoints_per_decade,niter;
  std::stringstream ss(line);
  std::size_t N_betas;
  double beta_min;
  double nu_beta;
  std::size_t nsamples_50_beta;
  M_Matrix<AP> aps;
  std::vector<std::vector<double>> apsPar;
  double maxTime;
  std::size_t samples, nskip;

  ss>>temperateS>>eviName>>experimentName>>priorName>>dx>>dtmin>>nPoints_per_decade>>dtmax>>niter>>maxduration>>landa0>>v>>nmaxloop>>initseed>>aps>>apsPar>>N_betas>>beta_min>>maxTime>>samples>>nskip>>nu_beta
      >>nsamples_50_beta>>pTjump;

  Adaptive_Beta
      aBeta(N_betas,beta_min,nu_beta,nsamples_50_beta);



  std::cout<<"evidenceS: "<<temperateS;
  std::cout<<" eviName: "<<eviName;
  std::cout<<" experimentName:"<<experimentName;
  std::cout<<" priorName: "<<priorName;
  std::cout<<" dx: "<<dx;
  std::cout<<" dtmin: "<<dtmin;
  std::cout<<" nPoints_per_decade: "<<nPoints_per_decade;
  std::cout<<" dtmax: "<<dtmax;
  std::cout<<" niter: "<<niter;
  std::cout<<" maxduration "<<maxduration;
  std::cout<<" landa0 "<<landa0;
  std::cout<<" v "<<v;
  std::cout<<" nmaxloop "<<nmaxloop;

  std::cout<<" initseed "<<initseed;
  std::cout<<" adaptive beta "<<aBeta;
  std::cout<<AP::ClassName()<<" "<<aps;
  for (std::size_t i=0; i<apsPar.size(); ++i)
    {
      std::cout<<AP::ParName(i)<<" "<<apsPar;
      for (std::size_t j=0; j<apsPar[i].size(); ++j)
        std::cout<<apsPar[i][j]<<" ";
    }
  std::cout<<" maxTime (h) "<<maxTime;

  std::cout<<" samples "<<samples;
  std::cout<<" nskip "<<nskip;
  std::cout<<" pTjump "<<pTjump;


  Adaptive_discrete<AP> landa(aps,apsPar);

  std::cout<<landa;

  std::cerr<<line;

  Experiment *e=new Experiment;
  e->load(experimentName);
  if (e->numMeasures()==0)
    {
      std::cerr<<"Experiment "<<experimentName<<" not found\n";
      return;
    }

  std::string filename=priorName;
  std::fstream f;

  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          std::cerr<<"Parameters file "<<filename<<" or "<<filenaExt<<" not found"<<std::endl;
          f.close();
          return;
        }
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters prior;

  if (!prior.read(line2,f))
    {
      std::cerr<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();


  BaseModel*m=BaseModel::create(prior);
  if (m!=nullptr)
    {
      cm_->push_back(m);

      auto CL=new CortexPoisonLikelihood(eviName+"_lik",e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio);

      MyModel<MyData> m(CL);
      MyData d(CL);
      Metropolis_Hastings_mcmc<
          MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa> mcmc;
      LevenbergMarquardt_step<MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa> LMLik;
      Poisson_DLikelihood<MyData,MyModel> DLik;
      TT tt;
      std::mt19937_64 mt;
      std::random_device rd;

      if (initseed==0)
        {
          std::mt19937_64::result_type seed=rd();

          mt.seed(seed);
          eviName+=time_now()+"_"+std::to_string(seed);
          std::cerr<<"\n random seed =\n"<<seed<<"\n";
        }
      else
        {
          std::mt19937_64::result_type seed=initseed;

          mt.seed(seed);
          eviName+=time_now()+"_"+std::to_string(seed);

          std::cerr<<"\n provided seed =\n"<<seed<<"\n";

        }


      std::string eviNameLog=eviName+"_log.txt";
      std::ofstream flog;

      flog.open(eviNameLog.c_str(), std::ofstream::out | std::ofstream::app);
      flog<<line<<"\n";
      flog<<"temperateS: "<<temperateS<<"\n";
      flog<<" eviName: "<<eviName<<"\n";
      flog<<" experimentName:"<<experimentName<<"\n";
      flog<<" priorName: "<<priorName<<"\n";
      flog<<" dx: "<<dx<<"\n";
      flog<<" dtmin: "<<dtmin<<"\n";
      flog<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
      flog<<" dtmax: "<<dtmax<<"\n";
      flog<<" niter: "<<niter<<"\n";
      flog<<" maxduration "<<maxduration<<"\n";
      flog<<" landa0 "<<landa0<<"\n";
      flog<<" v "<<v<<"\n";
      flog<<" nmaxloop "<<nmaxloop<<"\n";

      flog<<" initseed "<<initseed<<"\n";
      flog<<" adaptive beta "<<aBeta<<"\n";
      flog<<AP::ClassName()<<" "<<aps<<"\n";
      flog<<" maxTime (h) "<<maxTime<<"\n";

      flog<<" samples "<<samples<<"\n";
      flog<<" nskip "<<nskip<<"\n";
      flog<<" pTjump "<<pTjump<<"\n";

      flog.flush();


      std::chrono::steady_clock::time_point startTime=std::chrono::steady_clock::now();
      double timeOpt=0;


      typename TT::myEvidence * ev= tt.run(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,pTjump,mt,flog,startTime,timeOpt);
      std::cout<<*ev;
      flog<<*ev;
      flog.close();
      std::ofstream fout;
      fout.open(eviName.c_str(), std::ofstream::out | std::ofstream::app);
      fout<<line<<"\n";
      fout<<"evidenceS: "<<temperateS<<"\n";
      fout<<" eviName: "<<eviName<<"\n";
      fout<<" experimentName:"<<experimentName<<"\n";
      fout<<" priorName: "<<priorName<<"\n";
      fout<<" dx: "<<dx<<"\n";
      fout<<" dtmin: "<<dtmin<<"\n";
      fout<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
      fout<<" dtmax: "<<dtmax<<"\n";
      fout<<" niter: "<<niter<<"\n";
      fout<<" maxduration "<<maxduration<<"\n";
      fout<<" landa0 "<<landa0<<"\n";
      fout<<" v "<<v<<"\n";
      fout<<" nmaxloop "<<nmaxloop<<"\n";

      fout<<" initseed "<<initseed<<"\n";
      fout<<" adaptive beta "<<aBeta<<"\n";
      fout<<" maxTime (h) "<<maxTime<<"\n";

      fout<<" samples "<<samples<<"\n";
      fout<<" nskip "<<nskip<<"\n";

      fout<<*ev<<"\n";
      fout.close();
    }




}
