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
  std::string line;
  safeGetline(f,line);
  if (line.find("foto")!=line.npos)
    {auto  s=new TissueSection(filename);

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
  else  if (line.find("simulation")!=line.npos)
    while (f)
      {
        auto sim=new CortexSimulation;
        sim->read(line,f);
        if (sim->t_.size()>0)
          {
            cmd_->push_back(sim);
            std::string c="write  ";
            c+=sim->id_;

            cmd_->execute(c);
          }

      }

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


  CortexMeasure* m=cmd_->getMeasure(dataName);
  if (m!=nullptr)
    {
      std::string filename=dataName+"_histo.txt";
      std::ofstream f;
      f.open(filename.c_str(),std::ofstream::out);

      m->write(f);
      f.close();

    }

  CortexSimulation* sim=cmd_->getSimulation(dataName);
  if (sim!=nullptr)
    {
      std::string var;
      ss>>var;
      std::string par;
      std::getline(ss,par);
      if (var.empty())
        {

          std::string filename=dataName+"_sim.txt";
          std::ofstream f;
          f.open(filename.c_str(),std::ofstream::out);

          sim->write(f);
          f.close();
        }
      else
        {
          std::string filename=dataName+"_"+var+".txt";
          std::ifstream fi;
fi.open(filename.c_str(),std::ios_base::in);
          if (fi)
            {
              unsigned num=0;
              filename=dataName+"_"+std::to_string(num)+"_"+var+".txt";
              fi.close();
              fi.open(filename);
              while (fi)
                {
                  ++num;
                  filename=dataName+"_"+std::to_string(num)+"_"+var+".txt";
                  fi.close();
                  fi.open(filename);

                }
            }
          std::ofstream f;
          f.open(filename.c_str(),std::ofstream::out);

          sim->write(f,var,par);
          f.close();

        }

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
