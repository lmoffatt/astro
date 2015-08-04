#ifndef COMMANDMANAGER
#define COMMANDMANAGER
#include <string>
#include <map>
#include "CortexMeasure.h"


/*! @file CommandManager.h   Management of the commands
 *
 *
 * */



class BaseModel;
class CortexSimulation;



inline std::string& removeComments(std::string& line)
{
  auto pos=line.find("//");
  if (pos!=line.npos)
    line.erase(pos);
  return line;
}



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

  void push_back(BaseModel* model);

  void push_back(CortexSimulation* simulation);


  void push_back(CortexMeasure* measure)
  {
    measures[measure->id()]=measure;
  }


  ~CommandManager();

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


class Script
{
public:
  Script(CommandManager* cm):cm_(cm){}

  int run(char* filenme);

  void execute (std::string line);

private:
  CommandManager* cm_;


};


class readCommand: public CommandBase
{


  // CommandBase interface
public:
  virtual void run(const std::string rline);

  readCommand(CommandManager* cm):
    cmd_(cm){}

  virtual std::string id() const
  {
    return "read";
  }
private: CommandManager* cmd_;
};



class writeCommand: public CommandBase
{


  // CommandBase interface
public:
  virtual void run(const std::string line);

  writeCommand(CommandManager* cm):
    cmd_(cm){}

  virtual std::string id() const
  {
    return "write";
  }
private:
  CommandManager* cmd_;
};


std::istream& safeGetline(std::istream& is, std::string& t);





class AlignCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "align";
  }
  AlignCommand(CommandManager* cm):cm_(cm){}

private:
  CommandManager* cm_;

};

class MergeCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "merge";
  }
  MergeCommand(CommandManager* cm):cm_(cm){}

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
  DistancesCommand(CommandManager* cm):cm_(cm){}

private:
  CommandManager* cm_;

};



class HistogramCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "histogram";
  }
  HistogramCommand(CommandManager* cm):cm_(cm){}

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
  ExperimentCommand(CommandManager* cm):cm_(cm){}

private:
  CommandManager* cm_;

};



class SimulateCommand:public CommandBase
{
 // CommandBase interface
public:
  virtual void run(const std::string line);
  virtual std::string id() const
  {
    return "simulate";
  }
  SimulateCommand(CommandManager* cm):cm_(cm){}

private:
  CommandManager* cm_;

};











#endif // COMMANDMANAGER

