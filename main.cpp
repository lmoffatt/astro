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
        return s.runDefine(argv[1],argv[2],argv[3]);
      }
      break;

    default:
      break;
    }



  return 0;
}

