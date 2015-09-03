#ifndef POLYGON
#define POLYGON
#include "CortexMeasure.h"



class Polygon
{

  Polygon(const std::vector<position>& limits):limits_(limits){}

  bool isInside(const position& x)
  {

  }

private:
  const std::vector<position>* limits_;

};


#endif // POLYGON

