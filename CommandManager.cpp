#include "Models.h"
#include "CortexSimulation.h"
#include "CommandManager.h"
#include "LevenbergMarquardt.h"
#include "CortexLikelihood.h"
#include <sstream>
#include <iostream>

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

int Script::runDefine(const std::string &filename, const std::string &label, const std::string &valueInplace)
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

Experiment *CommandManager::getExperiment(std::string id)
{
  auto it=experiments.find(id);
  if(it!=experiments.end())
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



void AlignCommand::run(const std::string line)
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


void MergeCommand::run(const std::string line)
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








void DistancesCommand::run(const std::string line)
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




void HistogramCommand::run(const std::string line)
{

  //histogram 3dpl d_les  0:+100:3000

  std::string cName, dataName,kind;
  double start,dx,end;
  char colon;
  std::stringstream ss(line);
  ss>>cName>>dataName>>kind>>start>>colon>>dx>>colon>>end;

  TissueSection* s=cm_->getSection(dataName);

  std::vector<double> limits;

  double current=start;
  limits.push_back(current);
  while(current<end)
    {
      current+=dx;
      limits.push_back(current);

    }

  if (s!=nullptr)
    {
      CortexMeasure* m= s->measure(limits);
      cm_->push_back(m);
    }


}


void ExperimentCommand::run(const std::string line)

{
  // experiment experiment1 0:+200:4000   3dpl 7dpl

  std::string cName,expName, currsection;
  std::vector<std::string> sections;
  double start,dx,end;
  char colon;
  std::stringstream ss(line);
  ss>>cName>>expName>>start>>colon>>dx>>colon>>end;
  while(ss>>currsection)
    sections.push_back(currsection);
  std::vector<CortexMeasure> vCM;
  for (auto s:sections)
    {
      cm_->execute("read "+s);
      cm_->execute("align "+s);
      cm_->execute("merge "+s);
      cm_->execute("distances "+s);
      cm_->execute("histogram "+s+ " d_les "+std::to_string(start)+colon+std::to_string(dx)+colon+std::to_string(end) );
      vCM.push_back(*cm_->getMeasure(s));

    }
  Experiment* e=new Experiment(expName,vCM);
  cm_->push_back(e);

  std::ofstream fo;
  fo.open("out.txt");
  e->write(fo);
  fo.close();
  Experiment e2;
  std::ifstream fi;
  fi.open("out.txt");
  std::string line2;
  e2.read(line2,fi);

}


void SimulateCommand::run(const std::string line)
{
  std::string cName, eName, paramName,simName;
  double dt;
  std::stringstream ss(line);
  ss>>cName>>eName>>paramName>>dt>>simName;

  std::string filename=eName;
  std::ifstream f(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  std::string line2;
  safeGetline(f,line2);
  CortexExperiment e;

  e.read(line2,f);

  f.close();

  filename=paramName;
  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  safeGetline(f,line2);
  Parameters p;

  p.read(line2,f);

  f.close();

  BaseModel*m=BaseModel::create(p);
  if (m!=nullptr)
    {
      cm_->push_back(m);

      CortexSimulation* s=new CortexSimulation;
      *s=m->run(e,dt);
      if (simName.empty())
        {
          s->id_="sim_";
          s->id_+=paramName;
        }
      else
        s->id_=simName;

      cm_->push_back(s);
      std::string c="write  ";
      c+=s->id_;

      cm_->execute(c);

    }




}





void readCommand::run(const std::string rline)
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
    }
  std::string line;
  safeGetline(f,line);
  double dia=0;
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

          if (line.find("foto")!=line.npos)
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








void writeCommand::run(const std::string line)
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







void LikelihoodCommand::run(const std::string line)
{

  //likelihood lik experiment1 parameters_10 parameters_10 0.5 50 1000

  std::string cName, lName, eName, priorName, paramName;
  double dt,dx, tequilibrio;
  std::stringstream ss(line);
  ss>>cName>>lName>>eName>>priorName>>paramName>>dt>>dx>>tequilibrio;

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


  CortexLikelihood* CL= new CortexPoisonLikelihood(lName,e,prior,dt,tequilibrio);
  cm_->push_back(CL);

  CortexMultinomialLikelihoodEvaluation CE(*CL,p);
  std::ofstream fo;
  std::string fnameout=lName+"_lik.txt";
  fo.open(fnameout.c_str());
  CE.extract(fo);
  fo.close();

  }






std::string LikelihoodCommand::id() const
{
  return "likelihood";
}


void OptimizeCommand::run(const std::string line)
{

  //optimize opt experiment1 parameters_10 parameters_10 10 100

  std::string optimizeS, optName, experimentName, priorName, paramName;
  double dt,dx=50, tequilibrio=100000;
  double factor=0,probParChange=1;
  unsigned int initseed=0;
  std::size_t niter,nseeds=0;
  std::stringstream ss(line);

  ss>>optimizeS>>optName>>experimentName>>priorName>>paramName>>dt>>niter>>factor>>nseeds>>probParChange>>initseed;

  Experiment* e=cm_->getExperiment(experimentName);
  if (e==nullptr)
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
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters prior;

  prior.read(line2,f);

  f.close();

  filename=paramName;
  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  safeGetline(f,line2);
  Parameters p;

  p.read(line2,f);

  f.close();

  BaseModel*m=BaseModel::create(prior);
  if (m!=nullptr)
    {
      cm_->push_back(m);

      CortexPoisonLikelihood  CL(optName+"_lik",e,prior,dt,tequilibrio);



      LevenbergMarquardtDistribution LM(&CL,p,niter,optName);


      LM.optimize(optName,factor,nseeds,probParChange,initseed);

      Parameters* opt=new Parameters(LM.OptimParameters());
      std::string optfname=opt->save(optName);
      CortexMultinomialLikelihoodEvaluation CE(CL,LM.OptimParameters());
      std::ofstream fo;
      std::string fnameout=optfname;
      fnameout.insert(optfname.size()-5,"_lik.txt");
      fo.open(fnameout.c_str());
      CE.extract(fo);
      fo.close();





    }




}
