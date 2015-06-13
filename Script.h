#ifndef SCRIPT
#define SCRIPT
#include <string>
#include "CommandManager.h"

class Script
{
public:
  Script(CommandManager* cm):cm_(cm){}

  int run(char* filenme);

  void execute (std::string line);

private:
  CommandManager* cm_;


};


#endif // SCRIPT

