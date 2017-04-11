#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "CommandManager.h"

int main(int argc, char **argv)
{


  typedef myCommandManager<Cls,TissueSection,CortexMeasure,Experiment,
      BaseModel,CortexSimulation,CortexLikelihood> CM;

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
                  (C<CM>(),C<void>(),Co<Cls>(),&read<CM>,
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::string,std::string>{"","filename"},
                   std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"}));



  cm.push_command("tempering",
                  make_Command
                  (C<CM>(),C<void>(),Co<Cls>(),
                   &Tempering<CM>,
                   std::pair<CM*,std::string>{&cm,"CommandManager"},
                   std::pair<std::string,std::string>{"","eviName"},
                   std::pair<std::string,std::string>{"","experimentName"},
                   std::pair<std::string,std::string>{"","priorName"},
                   std::pair<double,std::string>{0,"dtmin"},
                   std::pair<double,std::string>{0,"dtmax"},
                   std::pair<double,std::string>{0,"dx"},
                   std::pair<double,std::string>{0,"tequilibrio"},
                   std::pair<double,std::string>{0,"maxduration"},
                   std::pair<double,std::string>{0,"landa0"},
                   std::pair<double,std::string>{0,"v"},
                   std::pair<double,std::string>{0,"pTjump"},
                   std::pair<std::size_t,std::string>{0,"nmaxloop"},
                   std::pair<std::mt19937_64::result_type ,std::string>{0,"initseed"},
                   std::pair<std::size_t,std::string>{0,"nPoints_per_decade"},
                   std::pair<std::size_t,std::string>{0,"niter"},
                   std::pair<std::size_t,std::string>{0,"N_betasInit"},
                   std::pair<std::size_t,std::string>{0,"N_betasInit"},
                   std::pair<double,std::string>{0,"beta_min"},
                   std::pair<double,std::string>{0,"mL0"},
                   std::pair<double,std::string>{0,"sL0"},
                   std::pair<double,std::string>{0,"mlogdL"},
                   std::pair<double,std::string>{0,"slogdL"},
                   std::pair<double,std::string>{0,"mmlogdelta"},
                   std::pair<double,std::string>{0,"smlogdelta"},
                   std::pair<double,std::string>{0,"mlogslogdelta"},
                   std::pair<double,std::string>{0,"slogslogdelta"},
                   std::pair<M_Matrix<Landa>,std::string>{{},"aps"},
                   std::pair<std::vector<std::vector<double>>,std::string>{{},"apsPar"},
                   std::pair<double,std::string>{0,"maxTime"},
                   std::pair<std::size_t,std::string>{0,"samples"},
                   std::pair<std::size_t,std::string>{0,"nskip"},
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

