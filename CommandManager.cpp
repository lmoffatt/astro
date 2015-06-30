#include "CommandManager.h"
#include <sstream>
#include "read.h"
#include "AlignCommand.h"

CommandManager::CommandManager()
{
 cmd_["read"]=new readCommand(this);
 cmd_["align"]=new AlignCommand(this);
 cmd_["write"]=new writeCommand(this);
 cmd_["merge"]=new MergeCommand(this);
 cmd_["distances"]=new DistancesCommand(this);
 cmd_["histogram"]=new HistogramCommand(this);


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

CortexMeasure *CommandManager::getMeasure(std::string id)
{
 auto it=measures.find(id);
 if(it!=measures.end())
   return it->second;
 else
   return nullptr;
}
