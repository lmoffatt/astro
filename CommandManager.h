#ifndef COMMANDMANAGER
#define COMMANDMANAGER
#include <string>
#include <map>
#include "BaseClass.h"
#include "CortexMeasure.h"

/*! @file CommandManager.h   Management of the commands
 *
 *
 * */







class BaseModel;
class CortexSimulation;
class CortexModelLikelihood;


inline std::string& removeComments(std::string& line)
{
  auto pos=line.find("//");
  if (pos!=line.npos)
    line.erase(pos);
  return line;
}



class CommandBase: public BaseClass
{
public:
  static std::string ClassName(){
    return "Command";
  }

   virtual std::string myClass()const override
  {
    return ClassName();
  }

  CommandBase(const std::string& id){setId(id);}

  virtual ~CommandBase(){}

  virtual void run(const std::string line)=0;

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

  void push_back(BaseModel* model);

  void push_back(CortexSimulation* simulation);


  void push_back(CortexMeasure* measure);

  void push_back(CortexModelLikelihood* likelihood);


  ~CommandManager();

  CortexMeasure *getMeasure(std::string id);
  BaseModel *getModel(std::string idModel);
  CortexSimulation *getSimulation(std::string idSimulation);
  Experiment *getExperiment(std::string id);

  CortexModelLikelihood* getLikelihood(const std::string& idLik);
  void push_back(Experiment *experiment);
private:
  std::map <std::string, CommandBase*> cmd_;
  std::map <std::string,TissueSection*> sections;

  std::map <std::string,CortexMeasure*> measures;

  std::map <std::string,Experiment*> experiments;

  std::map <std::string,BaseModel*> models;

  std::map <std::string,CortexSimulation*> simulations;

  std::map <std::string,CortexModelLikelihood*> likelihoods;


};


class Script
{
public:
  Script(CommandManager* cm):cm_(cm){}

  int run(char* filenme);

 private:
  CommandManager* cm_;


};


class readCommand: public CommandBase
{


  // CommandBase interface
public:
  virtual void run(const std::string rline);

  readCommand(CommandManager* cm):
    CommandBase("read"),
    cmd_(cm){}


  ~readCommand(){}

private: CommandManager* cmd_;
};



class writeCommand: public CommandBase
{


  // CommandBase interface
public:
  virtual void run(const std::string line);

  writeCommand(CommandManager* cm):
    CommandBase("write"),
    cmd_(cm){}

  ~writeCommand(){}


private:
  CommandManager* cmd_;
};







class AlignCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "align";
  }
  AlignCommand(CommandManager* cm):
    CommandBase("align"),
    cm_(cm){}

  ~AlignCommand(){}

private:
  CommandManager* cm_;

};

class MergeCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line) override;
  MergeCommand(CommandManager* cm):
    CommandBase("merge"),
    cm_(cm){}
  ~MergeCommand(){}

private:
  CommandManager* cm_;

};


class DistancesCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "distances";
  }
  DistancesCommand(CommandManager* cm):
    CommandBase("distances"),
    cm_(cm){}
  ~DistancesCommand(){}
private:
  CommandManager* cm_;

};



class HistogramCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  HistogramCommand(CommandManager* cm):
    CommandBase("histogram"),
    cm_(cm){}
  ~HistogramCommand(){}

private:
  CommandManager* cm_;

};


class ExperimentCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "experiment";
  }
  ExperimentCommand(CommandManager* cm):
    CommandBase("experiment"),
    cm_(cm){}
~ExperimentCommand(){}
private:
  CommandManager* cm_;

};



class SimulateCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);

  SimulateCommand(CommandManager* cm):
    CommandBase("simulate"),
    cm_(cm){}

  ~SimulateCommand(){}
private:
  CommandManager* cm_;

};


class LikelihoodCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const;
  LikelihoodCommand(CommandManager* cm):
    CommandBase("likelihood"),
    cm_(cm){}

private:
  CommandManager* cm_;

};


class OptimizeCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "optimize";
  }
  OptimizeCommand(CommandManager* cm):
    CommandBase("optimize"),
    cm_(cm){}

  ~OptimizeCommand(){}
private:
  CommandManager* cm_;

};




inline
std::string operator+(const std::string& one,const  std::string& two)
{
  return one+two;

}








#endif // COMMANDMANAGER

