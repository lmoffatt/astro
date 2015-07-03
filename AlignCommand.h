#ifndef ALIGNCOMMAND
#define ALIGNCOMMAND
#include "CommandManager.h"

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






#endif // ALIGNCOMMAND

