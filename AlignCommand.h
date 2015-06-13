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


#endif // ALIGNCOMMAND

