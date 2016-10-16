#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "CommandManager.h"

int main(int argc, char **argv)
{

 CommandManager cm;

  switch (argc) {
    case 1:

      break;
    case 2:
      {
        Script s(&cm);
        return s.run(argv[1]);
      }
      break;
    case 4:
      {
        Script s(&cm);
        return s.runDefine(argv[1],{argv[2]},{argv[3]});
      }
      break;

    default:
      {
        std::size_t nlabels=argc/2-1;
        std::vector<std::string> labels(nlabels);
        std::vector<std::string> valuesInPlace(nlabels);
        for (std::size_t i=0; i<nlabels; ++i)
          {
            labels[i]=argv[2*(i+1)];
            valuesInPlace[i]=argv[2*(i+1)+1];
          }
        Script s(&cm);
        return s.runDefine(argv[1],labels,valuesInPlace);


      }
      break;
    }



  return 0;
}

