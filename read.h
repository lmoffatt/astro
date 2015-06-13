#ifndef READ
#define READ
#include "CommandManager.h"


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



#endif // READ

