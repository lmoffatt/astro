#include "CortexSimulation.h"
#include "CommandManager.h"

#include <sstream>





std::ostream &CortexSimulation::write(std::ostream &s)
{
  s<<"simulation\n";
  s<<id_<<"\n";

  p_.write(s);

  s<<"dt simulation"<<"\t"<<dt_<<"\n";

  s<<"numSamples"<<"\t"<<t_.size()<<"\t";
  s<<"numNodes"<<"\t"<<x_.size()<<"\t";
  s<<"numStates"<<"\t"<<rho_.front().front().size()<<"\n"<<"\n";

  s<<"x"<<"\n";
  for (unsigned i=0; i<x_.size(); ++i)
    {
      s<<x_[i]<<"\t";
    }
  s<<"\n";

  s<<"dx"<<"\n";
  for (unsigned i=0; i<dx_.size(); ++i)
    {
      s<<dx_[i]<<"\t";
    }
  s<<"\n";



  s<<"psi total"<<"\n";
  for (unsigned i=0; i<t_.size(); ++i)
    {
      s<<t_[i]<<"\t";
      for (unsigned j=0; j<psi_T_[i].size(); ++j)
        {
          s<<psi_T_[i][j]<<"\t";
        }
      s<<"\n";
    }
  s<<"\n";


  s<<"psi bound"<<"\n";
  for (unsigned i=0; i<t_.size(); ++i)
    {
      s<<t_[i]<<"\t";
      for (unsigned j=0; j<psi_B_[i].size(); ++j)
        {
          s<<psi_B_[i][j]<<"\t";
        }
      s<<"\n";
    }
  s<<"\n";

  s<<"omega total"<<"\n";
  for (unsigned i=0; i<t_.size(); ++i)
    {
      s<<t_[i]<<"\t";
      for (unsigned j=0; j<omega_T_[i].size(); ++j)
        {
          s<<omega_T_[i][j]<<"\t";
        }
      s<<"\n";
    }
  s<<"\n";

  s<<"omega bound"<<"\n";
  for (unsigned i=0; i<t_.size(); ++i)
    {
      s<<t_[i]<<"\t";
      for (unsigned j=0; j<omega_B_[i].size(); ++j)
        {
          s<<omega_B_[i][j]<<"\t";
        }
      s<<"\n";
    }
  s<<"\n";



  for (unsigned k=0;k<rho_.front().front().size(); ++k)
    {
      s<<"rho"<<"\t"<<k<<"\n";
      for (unsigned i=0; i<t_.size(); ++i)
        {
          s<<t_[i]<<"\t";
          for (unsigned j=0; j<rho_[i].size(); ++j)
            {
              s<<rho_[i][j][k]<<"\t";
            }
          s<<"\n";
        }
      s<<"\n";
    }




  return s;
}

std::ostream &CortexSimulation::write(std::ostream &s, const std::string &var, const std::string &par)
{
  std::vector<double> xs;
  std::vector<double> ts;
  std::vector<double> ks;

  unsigned numK=rho_.front().front().size();


  if (par.empty())
    {
      xs=x_;
      ts=t_;

    }
  else
    {

      std::stringstream ss(par);
      double number;
      std::string var2;
      char equal;
      while (ss>>var2>>equal>>number)
        {
          if (var2=="t")
            {
              ts.push_back(number);
            }
          else if (var2=="x")
            {
              xs.push_back(number*1e-6);
            }
          else if (var2=="k")
            {
              ks.push_back(number);
            }

        }


    }

  if (var=="psi")
    {
      if (!xs.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<xs.size(); ++j)
            {
              double x=xs[j];
              while ((x_[jr]<x)&&jr<x_.size())
                ++jr;
              if (jr<x_.size())
                js.push_back(jr);
            }

          s<<"psi total and bound vs time at different positions"<<"\n";

          s<<"time"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<"T"<<x_[js[jjs]]<<"\tB"<<x_[js[jjs]]<<"\tF"<<x_[js[jjs]]<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<t_.size(); ++i)
            {
              s<<t_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<psi_T_[i][js[jjs]]<<"\t"<<psi_B_[i][js[jjs]]<<"\t";
                  s<<psi_T_[i][js[jjs]]-psi_B_[i][js[jjs]]<<"\t";

                }
              s<<"\n";
            }
          s<<"\n";

        }
      if (!ts.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<ts.size(); ++j)
            {
              double t=ts[j];
              while ((t_[jr]<t)&&jr<t_.size())
                ++jr;
              if (jr<t_.size())
                js.push_back(jr);
            }
          s<<"psi total and bound vs position at different times"<<"\n";

          s<<"position"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<"T"<<t_[js[jjs]]<<"\tB"<<t_[js[jjs]]<<"\tF"<<t_[js[jjs]]<<"\t";
            }
          s<<"\n";

          for (unsigned i=0; i<x_.size(); ++i)
            {
              s<<x_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<psi_T_[js[jjs]][i]<<"\t"<<psi_B_[js[jjs]][i]<<"\t";
                  s<<psi_T_[js[jjs]][i]-psi_B_[js[jjs]][i]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }

    }
  else if (var=="omega")
    {
      if (!xs.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<xs.size(); ++j)
            {
              double x=xs[j];
              while ((x_[jr]<x)&&jr<x_.size())
                ++jr;
              if (jr<x_.size())
                js.push_back(jr);
            }

          s<<"omega total and bound vs time at different positions"<<"\n";

          s<<"time"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<"T"<<x_[js[jjs]]<<"\tB"<<x_[js[jjs]]<<"\tF"<<x_[js[jjs]]<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<t_.size(); ++i)
            {
              s<<t_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<omega_T_[i][js[jjs]]<<"\t"<<omega_B_[i][js[jjs]]<<"\t";
                  s<<omega_T_[i][js[jjs]]-omega_B_[i][js[jjs]]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
      if (!ts.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<ts.size(); ++j)
            {
              double t=ts[j];
              while ((t_[jr]<t)&&jr<t_.size())
                ++jr;
              if (jr<t_.size())
                js.push_back(jr);
            }

          s<<"omega total and bound vs position at different times"<<"\n";

          s<<"position"<<"\t";
          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {

              s<<"T"<<t_[js[jjs]]<<"\tB"<<t_[js[jjs]]<<"\tF"<<t_[js[jjs]]<<"\t";
            }
          s<<"\n";

          for (unsigned i=0; i<x_.size(); ++i)
            {
              s<<x_[i]<<"\t";
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {

                  s<<omega_T_[js[jjs]][i]<<"\t"<<omega_B_[js[jjs]][i]<<"\t";
                  s<<omega_T_[js[jjs]][i]-omega_B_[js[jjs]][i]<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }

    }
  else if (var=="rho")
    {
      if (!xs.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<xs.size(); ++j)
            {
              double x=xs[j];
              while ((x_[jr]<x)&&jr<x_.size())
                ++jr;
              if (jr<x_.size())
                js.push_back(jr);
            }

          if (ks.empty())
            for (unsigned k=0;k<numK; ++k)
              ks.push_back(k);


          s<<"rho vs time at different states and positions"<<"\n";

          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {
              s<<"time"<<"\t";
              for (unsigned ik=0; ik<ks.size(); ++ik)
                {
                  unsigned k=ks[ik];

                  s<<"k="<<k<<",x="<<x_[js[jjs]]<<"\t";
                }
              s<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<t_.size(); ++i)
            {
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {
                  s<<t_[i]<<"\t";
                  for (unsigned ik=0; ik<ks.size(); ++ik)
                    {
                      unsigned k=ks[ik];
                      s<<rho_[i][js[jjs]][k]<<"\t";
                    }
                  s<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
      if (!ts.empty())
        {
          // find the set of js
          std::vector<unsigned> js;
          unsigned jr=0;
          for (unsigned j=0; j<ts.size(); ++j)
            {
              double t=ts[j];
              while ((t_[jr]<t)&&jr<t_.size())
                ++jr;
              if (jr<t_.size())
                js.push_back(jr);
            }

          if (ks.empty())
            for (unsigned k=0;k<numK; ++k)
              ks.push_back(k);


          s<<"rho vs position at different states and times"<<"\n";

          for (unsigned jjs=0; jjs<js.size(); ++jjs)
            {
              s<<"position"<<"\t";
              for (unsigned ik=0; ik<ks.size(); ++ik)
                {
                  unsigned k=ks[ik];

                  s<<"k="<<k<<",t="<<t_[js[jjs]]<<"\t";
                }
              s<<"\t";
            }
          s<<"\n";
          for (unsigned i=0; i<x_.size(); ++i)
            {
              for (unsigned jjs=0; jjs<js.size(); ++jjs)
                {
                  s<<x_[i]<<"\t";
                  for (unsigned ik=0; ik<ks.size(); ++ik)
                    {
                      unsigned k=ks[ik];
                      s<<rho_[js[jjs]][i][k]<<"\t";
                    }
                  s<<"\t";
                }
              s<<"\n";
            }
          s<<"\n";

        }
    }

  return s;
}



void CortexSimulation::read(std::string& line, std::istream &s, std::ostream& logs)
{
  std::string name;
  std::stringstream ss(line);
  ss>>name;
  if (name=="simulation")
    {
      name.clear();
      while (name.empty()&&safeGetline(s,line))
        {
          ss.str(line);
          ss.clear();
          ss>>name;
        }
      id_=name;
      name.clear();
      unsigned numSamples=0;
      unsigned numNodes;
      unsigned numStates;

      while (true)
        {
          while (name.empty()&&safeGetline(s,line))
            {
              ss.str(line);
              ss.clear();
              ss>>name;
            }
          if (name.find("numSamples")!=name.npos)
            {
              if  (ss>>numSamples>>name>>numNodes>>name>>numStates)
                {
                  *this=CortexSimulation(id_,numSamples,numNodes,numStates);
                }
              name.clear();

            }
          else if (name=="dt")
            {
              std::string sim;
              double dt;
              if(ss>>sim>>dt)
                dt_=dt;
            }
          else if (name=="parameters")
            {
              p_.read(line,s,logs);
            }
          else if (name=="x")
            {
              name.clear();
              safeGetline(s,line);
              ss.clear();
              ss.str(line);
              unsigned i=0;
              double myx;
              while (ss>>myx)
                {
                  x_[i]=myx;
                  ++i;
                }

            }
          else if (name.find("dx")!=name.npos)
            {
              name.clear();
              safeGetline(s,line);
              ss.clear();
              ss.str(line);
              unsigned i=0;
              while (ss>>dx_[i])
                {
                  ++i;
                }
            }
          else if (name.find("psi total")!=name.npos)
            {
              name.clear();
              unsigned i=0;
              while (true)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while (ss>>psi_T_[i][j])
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else if (name.find("psi bound")!=name.npos)
            {
              name.clear();
              unsigned i=0;
              while (true)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while (ss>>psi_B_[i][j])
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else if (name.find("omega total")!=name.npos)
            {
              name.clear();
              unsigned i=0;
              while (true)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while (ss>>omega_T_[i][j])
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else if (name.find("omega bound")!=name.npos)
            {
              name.clear();
              unsigned i=0;
              while (true)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while (ss>>omega_B_[i][j])
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else if (name.find("rho")!=name.npos)
            {
              name.clear();
              unsigned k;
              ss>>k;
              unsigned i=0;

              while (i<numSamples)
                {
                  safeGetline(s,line);
                  ss.clear();
                  ss.str(line);
                  ss>>t_[i];
                  unsigned j=0;
                  while ((ss>>rho_[i][j][k])&&j<numNodes)
                    {
                      ++j;
                    }
                  if (j>0)
                    ++i;
                  else break;
                }

            }
          else break;
        }


     isValid_=true;
    }
  else isValid_=false;
}



