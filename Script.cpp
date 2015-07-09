

#include "Script.h"
#include "read.h"
#include <fstream>
#include <sstream>



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





