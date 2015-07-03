#ifndef COMMANDMANAGER
#define COMMANDMANAGER
#include <string>
#include <map>
#include "CortexState.h"


class CommandBase
{
public:
  virtual void run(const std::string line)=0;

  virtual std::string id()const=0;
  virtual ~CommandBase(){}



};





class CommandManager
{
public:
  CommandManager();


  void execute(std::string line);

  CommandBase* Command(std::string commandName);
  TissueSection* getSection(std::string idSection);

  void push_back(TissueSection* section)
  {
    sections[section->id()]=section;
  }

  void push_back(BaseModel* model)
  {
    models[model->id()]=model;
  }

  void push_back(CortexSimulation* simulation)
  {
    simulations[simulation->id_]=simulation;
  }


  void push_back(CortexMeasure* measure)
  {
    measures[measure->id()]=measure;
  }


  ~CommandManager()
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

  }

  CortexMeasure *getMeasure(std::string id);
  BaseModel *getModel(std::string idModel);
  CortexSimulation *getSimulation(std::string idSimulation);
private:
  std::map <std::string, CommandBase*> cmd_;
  std::map <std::string,TissueSection*> sections;
  std::map <std::string,CortexMeasure*> measures;


  std::map <std::string,BaseModel*> models;

  std::map <std::string,CortexSimulation*> simulations;



};


#endif // COMMANDMANAGER

