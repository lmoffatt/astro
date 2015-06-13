#include "read.h"
#include "CortexState.h"
#include <fstream>
#include <sstream>




void readCommand::run(const std::string rline)
{
  std::stringstream ss(rline);
  std::string commandName;
  std::string filename;
  ss>>commandName>>filename;

  std::ifstream f(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
    }
  auto  s=new TissueSection(filename);
  std::string line;
  safeGetline(f,line);

  while (f)
    {

      if (line.find("foto")!=line.npos)
        {
          TissuePhoto foto;

          foto.read(line,f);
          s->fotos[foto.num]=foto;
        }
      else
        safeGetline(f,line);

    }
   cmd_->push_back(s);





}








void writeCommand::run(const std::string line)
{
  std::string writeName, dataName;
  std::stringstream ss(line);
  ss>>writeName>>dataName;

  TissueSection* s=cmd_->getSection(dataName);
  if (s!=nullptr)
    {
      std::string filename=dataName+"_astro.txt";
      std::ofstream f;
      f.open(filename.c_str(),std::ofstream::out);

          for (auto fotop:s->fotos)
            fotop.second.write(f);

      f.close();

    }
  }




std::istream &safeGetline(std::istream &is, std::string &t)
{
  std::getline(is,t);
  auto it=t.find('\r');
  if (it!=t.npos)
    t.erase(it);
  return is;
}
