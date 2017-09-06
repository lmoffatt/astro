#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "CommandManager.h"

int main(int argc, char **argv)
{
  std::cerr<<argv[0]<<"\n";
  std::cerr<<argv[1]<<"\n";


  typedef myCommandManager<Cls,Cs<bool,CortexMeasure>,Cs<TissueSection,CortexMeasure,Experiment,
      BaseModel,CortexSimulation,CortexLikelihood>> CM;

  CM cm;

  /* cmd_["read"]=new readCommand(this);
  cmd_["align"]=new AlignCommand(this);
  cmd_["write"]=new writeCommand(this);
  cmd_["merge"]=new MergeCommand(this);
  cmd_["distances"]=new DistancesCommand(this);
  cmd_["histogram"]=new HistogramCommand(this);
  cmd_["simulate"]=new SimulateCommand(this);
  cmd_["experiment"]=new ExperimentCommand(this);
  cmd_["likelihood"]=new LikelihoodCommand(this);
  cmd_["optimize"]=new OptimizeCommand(this);
  cmd_["evidence"]=new EvidenceCommand(this);
  cmd_["tempering"]=new TemperingCommand(this);
*/

  cm.push_command("read",
                  make_Command
                  (C<CM>(),C<void>(),Co<Cls>(),read<CM>(),
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::string,std::string>{"","filename"},
                   std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"}));

  cm.push_command("write",
                  make_Command
                  (C<CM>(),C<bool>(),Co<Cls>(),write<CM>(),
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::string,std::string>{"","id"},
                   std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"},
                   std::pair<std::string,std::string>{"","filename"},
                   std::pair<bool,std::string>{false,"append"}
                   ));
  cm.push_command("dataFrame",
                  make_Command
                  (C<CM>(),C<bool>(),Co<Cls>(),dataFrame<CM>(),
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::string,std::string>{"","id"},
                   std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"},
                   std::pair<std::string,std::string>{"","filename"}
                   ));


  cm.push_command("measure",
                  make_Command
                  (C<CM>(),C<void>(),Co<Cls>(),
                   Measure<CM>(),
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::vector<TissueSection*>,std::string>{{},"input"},
                   std::pair<bool,std::string>{true,"average"},
                   std::pair<std::size_t,std::string>{0,"maxnumpoints"},
                   std::pair<std::size_t,std::string>{0,"initseed"},
                   std::pair<double,std::string>{0,"minDistance_Tissue"},
                   std::pair<double,std::string>{0,"minDistance_Vaso"},
                   std::pair<std::vector<double>,std::string>{0,"limits"},
                   std::pair<std::string,std::string>{"","experimentName"}
                   ));


  cm.push_command("tempering",
                  make_Command
                  (C<CM>(),C<void>(),Co<Cls>(),
                   Tempering<CM>(),
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::string,std::string>{"","eviName"},
                   std::pair<std::string,std::string>{"","experimentName"},
                   std::pair<std::string,std::string>{"","priorName"},
                   std::pair<std::string,std::string>{"","stateFile"},
                   std::pair<double,std::string>{0,"dtmin0"},
                   std::pair<double,std::string>{0,"dtmin"},
                   std::pair<double,std::string>{0,"dtmax"},
                   std::pair<double,std::string>{0,"dx"},
                   std::pair<bool,std::string>{false,"adapt_dt"},
                   std::pair<double,std::string>{0,"maxlogError"},
                   std::pair<double,std::string>{0,"dtinf"},
                   std::pair<bool,std::string>{false,"CrankNicholson"},
                   std::pair<double,std::string>{0,"f_maxlogError"},
                   std::pair<std::size_t,std::string>{0,"maxloop"},
                   std::pair<bool,std::string>{false,"UseDerivative"},
                   std::pair<double,std::string>{0,"tequilibrio"},
                   std::pair<double,std::string>{0,"maxduration"},
                   std::pair<double,std::string>{0,"landa0"},
                   std::pair<double,std::string>{0,"v"},
                   std::pair<double,std::string>{0,"pTjump"},
                   std::pair<double,std::string>{0,"slogL_max"},
                   std::pair<std::size_t,std::string>{0,"ndts_max"},
                   std::pair<std::size_t,std::string>{0,"nmaxloop"},
                   std::pair<std::mt19937_64::result_type ,std::string>{0,"initseed"},
                   std::pair<std::size_t,std::string>{0,"nPoints_per_decade"},
                   std::pair<std::size_t,std::string>{0,"niter"},
                   std::pair<std::size_t,std::string>{0,"N_betas"},
                   std::pair<double,std::string>{0,"beta_min"},
                   std::pair<std::size_t,std::string>{0,"N_betas_2"},
                   std::pair<double,std::string>{0,"beta_infimo"},
                   std::pair<std::string,std::string>{"","Landa_algorithm"},
                   std::pair<double,std::string>{0,"targetProb"},
                   std::pair<M_Matrix<Landa>,std::string>{{},"aps"},
                   std::pair<std::vector<std::vector<double>>,std::string>{{},"apsPar"},
                   std::pair<double,std::string>{0,"maxTime"},
                   std::pair<std::size_t,std::string>{0,"samples"},
                   std::pair<std::size_t,std::string>{0,"nskip"},
                   std::pair<std::size_t,std::string>{10000000,"maxSimFileSize"},
                   std::pair<bool,std::string>{false,"does_stdout"},

                   std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"}));




  switch (argc) {
    case 1:

      break;
    case 2:
      {
        myScript<CM> s(&cm);
        return s.run(argv[1], std::cerr);
      }
      break;
    case 4:
      {
        myScript<CM> s(&cm);
        return s.runDefine(argv[1],{argv[2]},{argv[3]}, std::cerr);
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
      myScript<CM> s(&cm);
      return s.runDefine(argv[1],labels,valuesInPlace, std::cerr);


    }
    break;
}



return 0;
}

