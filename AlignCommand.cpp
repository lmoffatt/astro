#include "AlignCommand.h"
#include <sstream>



void AlignCommand::run(const std::string line)
{

  std::string alignName, dataName;
  std::stringstream ss(line);
  ss>>alignName>>dataName;

  TissueSection* s=cm_->getSection(dataName);
  if (s!=nullptr)
    {
      s->align();
    }




}


void MergeCommand::run(const std::string line)
{
  std::string cName, dataName;
  std::stringstream ss(line);
  ss>>cName>>dataName;

  TissueSection* s=cm_->getSection(dataName);
  if (s!=nullptr)
    {
      s->merge();
    }



}








void DistancesCommand::run(const std::string line)
{
  std::string cName, dataName;
  std::stringstream ss(line);
  ss>>cName>>dataName;

  TissueSection* s=cm_->getSection(dataName);
  if (s!=nullptr)
    {
      s->distances();
    }



}




void HistogramCommand::run(const std::string line)
{

  //histogram 3dpl d_les  0:+100:3000

  std::string cName, dataName,kind;
  double start,dx,end;
  char colon;
  std::stringstream ss(line);
  ss>>cName>>dataName>>kind>>start>>colon>>dx>>colon>>end;

  TissueSection* s=cm_->getSection(dataName);

  std::vector<double> limits;

  double current=start;
  limits.push_back(current);
  while(current<end)
    {
      current+=dx;
      limits.push_back(current);

    }

  if (s!=nullptr)
    {
      CortexMeasure* m= s->measure(limits);
      cm_->push_back(m);
    }


}


void ExperimentCommand::run(const std::string line)

{

  //experiment  sim_zero


  std::string cName, dataName;
  std::stringstream ss(line);
  ss>>cName>>dataName;

  TissueSection* s=cm_->getSection(dataName);

  std::vector<double> limits;

  double current=start;
  limits.push_back(current);
  while(current<end)
    {
      current+=dx;
      limits.push_back(current);

    }

  if (s!=nullptr)
    {
      CortexMeasure* m= s->measure(limits);
      cm_->push_back(m);
    }


}
