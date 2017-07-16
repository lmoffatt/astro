#ifndef COMMANDMANAGER
#define COMMANDMANAGER
#include <string>
#include <map>
#include "BaseClass.h"
#include "CortexMeasure.h"
#include "CortexSimulation.h"
#include "myCommandManagement.h"
#include "CortexLikelihood.h"
#include "Evidence.h"


#include <experimental/type_traits>
/*! @file CommandManager.h   Management of the commands
 *
 *
 * */






class BaseModel;
class CortexSimulation;
class CortexLikelihood;
//class mcmc;

inline std::string& removeComments(std::string& line)
{
  auto pos=line.find("//");
  if (pos!=line.npos)
    line.erase(pos);
  return line;
}


inline std::string& replaceLabel(std::string& line,
                                 const std::string& label,
                                 const std::string& replacement)
{
  std::size_t i=0;
  while (i!=line.npos)
    {
      i=line.find(label,i);
      if (i!=line.npos)
        line.replace(i,label.size(),replacement);
    }
  return line;
}


inline std::string& replaceLabel(std::string& line,
                                 const std::vector<std::string>& label,
                                 const std::vector<std::string>& replacement)
{
  for (std::size_t j=0; j<label.size(); ++j)
    {
      std::size_t i=0;
      while (i!=line.npos)
        {
          i=line.find(label[j],i);
          if (i!=line.npos)
            line.replace(i,label[j].size(),replacement[j]);
        }
    }
  return line;
}


class CommandBase: public BaseClass
{
public:
  static std::string ClassName(){
    return "Command";
  }

  virtual std::string myClass()const override
  {
    return ClassName();
  }

  CommandBase(const std::string& id){setId(id);}

  virtual ~CommandBase(){}

  virtual void run(const std::string& line, std::ostream& logs)=0;

};


class CommandManager
{
public:
  CommandManager();


  void execute(std::string line, std::ostream& logs);

  CommandBase* Command(std::string commandName);
  TissueSection* getSection(std::string idSection);

  void push_back(TissueSection* section)
  {
    sections[section->id()]=section;
  }

  void push_back(BaseModel* model);

  void push_back(CortexSimulation* simulation);


  void push_back(CortexMeasure* measure);

  void push_back(CortexLikelihood* likelihood);


  ~CommandManager();

  CortexMeasure *getMeasure(std::string id);
  BaseModel *getModel(std::string idModel);
  CortexSimulation *getSimulation(std::string idSimulation);
  Experiment *getExperiment(const std::string &id);

  CortexLikelihood* getLikelihood(const std::string& idLik);
  void push_back(Experiment *experiment);
  // void push_back(mcmc *mcmc);
  //mcmc *getMcmc(const std::string& id);
private:
  std::map <std::string, CommandBase*> cmd_;
  std::map <std::string,TissueSection*> sections;

  std::map <std::string,CortexMeasure*> measures;

  std::map <std::string,Experiment*> experiments;

  //std::map <std::string,mcmc*> mcmcs;

  std::map <std::string,BaseModel*> models;

  std::map <std::string,CortexSimulation*> simulations;

  std::map <std::string,CortexLikelihood*> likelihoods;

};


struct has_ClassName_tag{};

struct has_not_ClassName_tag{};


template<class, typename T=void>
struct has_ClassName_traits
{
  typedef has_not_ClassName_tag tag;
};

template<class T>
struct has_ClassName_traits<T,decltype (T::ClassName())>
{
  typedef has_ClassName_tag tag;

};



template<typename>
std::string ClassName_imp(has_not_ClassName_tag)
{
  return "unknown";
}



template<typename T>
auto ClassName_imp(has_ClassName_tag)->decltype (T::ClassName())
{
  return T::ClassName();
}





inline
bool read_from_stream(std::istream&,std::ostream&,... )
{
  return true;
}




template<typename T>
auto read_from_stream(std::istream& is,std::ostream& /*logstream*/,T& x)
->decltype (bool(is>>x))
{

  if (is>>x)
    return true;
  else
    return false;
}

template<typename T>
auto read_from_stream(std::istream& is,std::ostream& logstream,T*& x)
->decltype (bool(is>>*x))
{

  return bool(is>>*x);
}


template<typename T>
auto read_from_stream(std::istream& is,std::ostream& logstream,T*& x)
->decltype (x->read(std::string(),is,logstream))
{
  std::string s;
  return x->read(s,is,logstream);
}

template<typename T>
auto read_from_stream(std::istream& is,std::ostream& logstream,T& x)
->decltype (x.read(std::string(),is,logstream))
{
  std::string s;
  return x.read(s,is,logstream);
}


inline
bool write_to_stream(std::ostream&,std::ostream&,... )
{
  return true;
}




template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& /*logstream*/,T& x)
->decltype (bool(os<<x))
{

  if (os<<x) return true;
  else return false;
}

template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& /*logstream*/,T*& x)
->decltype (bool(os<<*x))
{

  return bool(os<<*x);
}


template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& logstream,T*& x)
->decltype (x->write(std::string(),os,logstream))
{
  std::string s;
  return x->write(s,os,logstream);
}

template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& logstream,T& x)
->decltype (x.write(std::string(),os,logstream))
{
  std::string s;
  return x.write(s,os,logstream);
}







template<class T>
class Cls
{
public:
  static std::string name()
  {
    return ClassName_imp<T>(typename has_ClassName_traits<T>::tag());
  }
  static
  bool read(std::istream& is,T& x,std::ostream& logstream)
  {
    return read_from_stream(is,logstream,x);
  }

  static
  bool write (std::ostream& os, const T& x, std::ostream& logstream)
  {
    return write_to_stream(os,logstream,x);

  }


};



template<class T>
class Cls<T*>
{
public:
  static std::string name()
  {
    return ClassName_imp<T>(typename has_ClassName_traits<T>::tag());
  }
  static
  bool read(std::istream& is,T*& x,std::ostream& logstream)
  {
    return read_from_stream(is,logstream,x);
  }
  static
  bool write (std::ostream& os, const T* x, std::ostream& logstream)
  {
    return write_to_stream(os,logstream,x);

  }

};



class Script
{
public:
  Script(CommandManager* cm):cm_(cm){}

  int run(char* filename, std::ostream &logs);


  int runDefine(const std::string& filename,
                const std::vector<std::__cxx11::string> &label,
                const std::vector<std::__cxx11::string> &valueInplace,
                std::ostream &logs);


private:
  CommandManager* cm_;


};


template<class Command_Manager>
class myScript
{
public:
  myScript(Command_Manager* cm):cm_(cm){}

  int run(char* filename, std::ostream &logs)
  {
    std::ifstream f(filename);
    while (f)
      {
        std::string line;
        safeGetline(f,line);
        removeComments(line);
        std::cerr<<line;
        cm_->execute(line,logs);
      }
    return 0;

  }


  int runDefine(const std::string& filename,
                const std::vector<std::__cxx11::string> &label,
                const std::vector<std::__cxx11::string> &valueInplace,
                std::ostream &logs)
  {
    std::ifstream f(filename);
    while (f)
      {
        std::string line;
        safeGetline(f,line);
        removeComments(line);

        replaceLabel(line,label,valueInplace);
        cm_->execute(line,logs);
      }
    return 0;

  }


private:
  Command_Manager* cm_;


};








class readCommand: public CommandBase
{


  // CommandBase interface
public:
  virtual void run(const std::string& rline, std::ostream& logs) override;

  readCommand(CommandManager* cm):
    CommandBase("read"),
    cmd_(cm){}


  ~readCommand(){}

private: CommandManager* cmd_;
};



class writeCommand: public CommandBase
{


  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs) override;

  writeCommand(CommandManager* cm):
    CommandBase("write"),
    cmd_(cm){}

  ~writeCommand(){}


private:
  CommandManager* cmd_;
};







class AlignCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs);
  virtual std::string id() const
  {
    return "align";
  }
  AlignCommand(CommandManager* cm):
    CommandBase("align"),
    cm_(cm){}

  ~AlignCommand(){}

private:
  CommandManager* cm_;

};

class MergeCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs) override;
  MergeCommand(CommandManager* cm):
    CommandBase("merge"),
    cm_(cm){}
  ~MergeCommand(){}

private:
  CommandManager* cm_;

};


class DistancesCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs);
  virtual std::string id() const
  {
    return "distances";
  }
  DistancesCommand(CommandManager* cm):
    CommandBase("distances"),
    cm_(cm){}
  ~DistancesCommand(){}
private:
  CommandManager* cm_;

};




class ExperimentCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs);
  virtual std::string id() const
  {
    return "experiment";
  }
  ExperimentCommand(CommandManager* cm):
    CommandBase("experiment"),
    cm_(cm){}
  ~ExperimentCommand(){}
private:
  CommandManager* cm_;

};



class SimulateCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs);

  SimulateCommand(CommandManager* cm):
    CommandBase("simulate"),
    cm_(cm){}

  ~SimulateCommand(){}
private:
  CommandManager* cm_;

};


class LikelihoodCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs);
  virtual std::string id() const;
  LikelihoodCommand(CommandManager* cm):
    CommandBase("likelihood"),
    cm_(cm){}

private:
  CommandManager* cm_;

};


class OptimizeCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::string& line, std::ostream& logs);

  virtual std::string id() const
  {
    return "optimize";
  }
  OptimizeCommand(CommandManager* cm):
    CommandBase("optimize"),
    cm_(cm){}

  ~OptimizeCommand(){}
private:
  CommandManager* cm_;

};


class EvidenceCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::__cxx11::string &line, std::ostream& logs);

  virtual std::string id() const
  {
    return "evidence";
  }
  EvidenceCommand(CommandManager* cm):
    CommandBase("evidence"),
    cm_(cm){}

  ~EvidenceCommand(){}
private:
  CommandManager* cm_;

};
template<typename T>
class WriteIt
{
public:
  static
  void apply(const T& x,
             const std::string& idname,
             std::ostream* f,
             std::ostream* logs)
  {
    *f<<idname<<"\n"<<Cls<T>::name()<<"\n";
    Cls<T>::write(*f,x,*logs);
    *f<<"\n";
  }

};


template<typename,class=void>
struct doesWriteDataFrame: std::false_type{};

template< class ... > using void_t = void;


template<typename T>
struct doesWriteDataFrame<T, void_t<decltype(std::declval<T>().writeDataFrame(std::declval<std::ostream&>())) >> : std::true_type{};



template <bool, typename T>
struct writeDataFrame
{
  static void write(const T&,std::ostream*, std::ostream* log)
  {
    *log<<Cls<T>::name()<<" writeDataFrame is not implemented\n";

  }

};

template <typename T>
struct writeDataFrame<true,T>
{
  static void write(const T& x,std::ostream* f, std::ostream* log)
  {
    x.writeDataFrame(*f);
    *f<<"\n";
  }

};

template <typename T>
struct writeDataFrame<true,T*>
{
  static void write(const T* x,std::ostream* f, std::ostream* /*log*/)
  {
    x->writeDataFrame(*f);
    *f<<"\n";

  }

};





template <typename T>
class DataFrameIt
{
public:
  static void apply(const T& x,std::ostream* f,std::ostream* log)
  {
    writeDataFrame<doesWriteDataFrame<std::remove_pointer_t<T>>::value,T>::write(x,f,log);
  }

};








template<class Cm>
struct read
{
void operator()(Cm* cm,const std::string& filename, std::ostream* logs  )const
{
  std::cerr<<"here\n";
  std::ifstream f(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          *logs<<"could not open file"<<filenaExt<<"\n";
          f.close();
          std::string below="../"+filenaExt;
          f.open(below.c_str(),std::ios_base::in);
          if (!f)
            {
              *logs<<"could not open file"<<below<<"\n";
            }
          else
            {
              *logs<<"file opened successfully "<<below<<"\n";

            }
        }
      *logs<<"file opened successfully "<<filenaExt<<"\n";

    }
  std::string line;
  safeGetline(f,line);
  double dia=0;
  std::size_t rata=0;
  if (line.find("experiment")!=line.npos)
    {
      safeGetline(f,line);
      while (f && line.empty())
        safeGetline(f,line);

      if (line.find("dia")!=line.npos)
        {
          std::stringstream ss(line);
          std::string diastr;

          ss>>diastr>>dia;
          safeGetline(f,line);

        }


      auto  s=new TissueSection(filename,dia);

      while (f)
        {

          if (line.find("rata")!=line.npos)
            {
              std::stringstream ss(line);
              std::string ratastr;

              ss>>ratastr>>rata;
              TissuePhoto* foto=new TissuePhoto;
              foto->dia_=dia;
              foto->read(line,f,*logs);
              s->fotos[foto->rata]=foto;

            }
          else if (line.find("foto")!=line.npos)
            {
              TissuePhoto* foto=new TissuePhoto;;
              foto->dia_=dia;
              foto->read(line,f,*logs);
              s->fotos[foto->num]=foto;
            }
          else
            safeGetline(f,line);

        }
      cm->push_back(s->id(),s);

    }
  else  if (line.find("simulation")!=line.npos)
    while (f)
      {
        auto sim=new CortexSimulation;
        sim->read(line,f,*logs);
        if (sim->t_.size()>0)
          {
            cm->push_back(sim->id_,sim);
          }

      }

}

};

template<class Cm>
struct write{
bool operator()(Cm* cm,const std::string& idname, std::ostream* logs,  std::string fname, bool append) const


{
  if (fname.empty())
    fname=idname+"_out.txt";
  std::ofstream f;
  if (append)
    f.open(fname, std::ofstream::app | std::ofstream::out);
  else
    f.open(fname,  std::ofstream::out);

  if (!f)
    {
      *logs<<"could not open file "+fname+" for writing";
      return false;
    }
  else
    {
      f<<idname<<"\n";
      cm->apply(Co<WriteIt>(),idname,idname,&f,logs);
      if (f)
        {
          *logs<<idname<<" written in file "<<fname;
          f.close();
          return true;
        }
      else
        {
          *logs<<idname<<" something wrong when written in file "<<fname;
          f.close();
          return false;
        }

    }
}

};

template<class Cm>
struct dataFrame
{
bool operator()(Cm* cm,const std::string& idname, std::ostream* logs,  std::string fname) const
{
  if (fname.empty())
    fname=idname+"_data_frame.txt";
  std::ofstream f;
  f.open(fname,  std::ofstream::out);

  if (!f)
    {
      *logs<<"could not open file "+fname+" for writing";
      return false;
    }
  else
    {
      f<<idname<<"\n";
      cm->apply(Co<DataFrameIt>(),idname,&f,logs);
      if (f)
        {
          *logs<<idname<<" written in file "<<fname;
          f.close();
          return true;
        }
      else
        {
          *logs<<idname<<" something wrong when written in file "<<fname;
          f.close();
          return false;
        }

    }
}
};



template<class Cm>
struct Measure
{
void operator()(Cm* cm_,
             std::vector<TissueSection *> sections,
             bool average,
             std::size_t maxnumpoints,
             std::size_t initseed,
             double minDistance_Tissue,
             double minDistance_Vaso,
             std::vector<double> limits, std::string expName) const
{
  std::mt19937_64 mt;
  if (initseed==0)
    {
      std::random_device rd;
      auto seed=rd();
      mt.seed(seed);
      std::cout<<"experiment uses seed="<<seed<<std::endl;
    }
  else if (initseed==1)
    {
      mt.seed();
      std::cout<<"experiment uses default seed"<<std::endl;
    }
  else
    {
      mt.seed(initseed);
      std::cout<<"experiment uses provided seed="<<initseed<<std::endl;

    }


  std::vector<CortexMeasure> vCM;
  for (auto s:sections)
    {
      auto p=s->measure(mt,limits,minDistance_Tissue, minDistance_Vaso,maxnumpoints);
      if (average)
        {
          for (auto it= p.begin(); it!= p.end(); ++it)
            {
              auto it2= vCM.begin();
              while(it2!=vCM.end())
                {
                  if (it2->add(*it))
                    break;
                  else
                    ++it2;
                }
              if (it2==vCM.end())
                vCM.push_back(*it);
            }
        }
      else
        {
          vCM.insert(vCM.begin(),p.begin(),p.end());
        }
    }


  auto e=new Experiment(expName,vCM,{});
  std::ofstream fo;
  fo.open(expName.c_str());
  e->write(fo);
  cm_->push_back(e->id(),e);
}

};



template<class Cm>
struct Tempering
{
void operator ()(Cm* cm_,
               std::string eviName,
               std::string experimentName,
               std::string priorName,
               double dtmin,
               double dtmax,
               double dx,
               double tequilibrio,
               double  maxduration,
               double landa0,
               double v,
               double pTjump,
               std::size_t nmaxloop,
               std::mt19937_64::result_type initseed,
               std::size_t nPoints_per_decade,
               std::size_t niter,
               std::size_t N_betasInit,
               std::size_t N_betas_max,
               double beta_min,
               double mL0,
               double sL0,
               double mmlogdL,
               double smlogdL,

               double mlogslogdL,
               double slogslogdL,

               double mloglandalogdL,
               double sloglandalogdL,

               double mlogepsilonlogdL,
               double slogepsilonlogdL,

               double mmlogdelta,
               double smlogdelta,
               double mlogslogdelta,
               double slogslogdelta,
               M_Matrix<Landa> aps,
               std::vector<std::vector<double>> apsPar,
               double maxTime,
               std::size_t samples,
               std::size_t nskip
               , std::ostream *logs)const
{

  // evidence evi experiment1 paramters_10
  //optimize opt experiment1 parameters_10 parameters_10 10 100
  typedef Landa AP;


  Master_Tempering_Likelihood::Prior MTLP
      ({},mL0,sL0,mmlogdL,smlogdL,mlogslogdL,slogslogdL,mloglandalogdL,
       sloglandalogdL,mlogepsilonlogdL,slogepsilonlogdL,
       mmlogdelta,smlogdelta,mlogslogdelta,slogslogdelta);
  Master_Adaptive_Beta
      aBeta(MTLP,N_betasInit,N_betas_max,beta_min,1);



  std::cout<<"\n eviName: "<<eviName;
  std::cout<<"\n experimentName:"<<experimentName;
  std::cout<<"\n priorName: "<<priorName;
  std::cout<<"\n dx: "<<dx;
  std::cout<<"\n dtmin: "<<dtmin;
  std::cout<<"\n nPoints_per_decade: "<<nPoints_per_decade;
  std::cout<<"\n dtmax: "<<dtmax;
  std::cout<<"\n niter: "<<niter;
  std::cout<<"\n maxduration "<<maxduration;
  std::cout<<"\n landa0 "<<landa0;
  std::cout<<"\n v "<<v;
  std::cout<<"\n nmaxloop "<<nmaxloop;

  std::cout<<"\n initseed "<<initseed;
  std::cout<<"\n adaptive beta "<<aBeta;
  std::cout<<AP::ClassName()<<" "<<aps;
  for (std::size_t i=0; i<apsPar.size(); ++i)
    {
      std::cout<<AP::ParName(i)<<" "<<apsPar;
      for (std::size_t j=0; j<apsPar[i].size(); ++j)
        std::cout<<apsPar[i][j]<<" ";
    }
  std::cout<<"\n maxTime (h) "<<maxTime;

  std::cout<<"\n samples "<<samples;
  std::cout<<"\n nskip "<<nskip;
  std::cout<<"\n pTjump "<<pTjump;


  Adaptive_discrete<AP> landa(aps,apsPar);

  std::cout<<landa;


  Experiment *e=new Experiment;
  e->load(experimentName,*logs);
  if (e->numMeasures()==0)
    {
      *logs<<"Experiment "<<experimentName<<" not found\n";
      return;
    }

  std::string filename=priorName;
  std::fstream f;

  f.open(filename.c_str());
  if (!f)
    {
      std::string filenaExt=filename+".txt";
      f.open(filenaExt.c_str());
      if (!f)
        {
          *logs<<"Parameters file "<<filename<<" or "<<filenaExt<<" not found"<<std::endl;
          f.close();
          return;
        }
    }
  std::string line2;
  safeGetline(f,line2);
  Parameters prior;

  if (!prior.read(line2,f,*logs))
    {
      *logs<<"File "<<filename<<" is not a Parameters file"<<std::endl;
      f.close();
      return;
    }

  f.close();


  BaseModel*m=BaseModel::create(prior);
  if (m!=nullptr)
    {
      cm_->push_back(m->id(),m);

      auto CL=new CortexPoisonLikelihood(eviName+"_lik",e,prior,dx,dtmin,nPoints_per_decade,dtmax,tequilibrio);

      MyModel<MyData> m(CL);
      MyData d(CL);
      Metropolis_Hastings_mcmc<
          MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa> mcmc;
      LevenbergMarquardt_step<MyData,MyModel,Poisson_DLikelihood,LM_MultivariateGaussian,Landa> LMLik;
      Poisson_DLikelihood<MyData,MyModel> DLik;
      TT tt;

      std::random_device rd;

      std::mt19937_64::result_type seed;
      if (initseed==0)
        {
           seed=rd();

          eviName+=time_now()+"_"+std::to_string(seed);
          *logs<<"\n random seed =\n"<<seed<<"\n";
        }
      else
        {
          seed=initseed;

          eviName+=time_now()+"_"+std::to_string(seed);

          *logs<<"\n provided seed =\n"<<seed<<"\n";

        }

      std::string eviNameLog=eviName+"_log.txt";
      std::ofstream flog;

      flog.open(eviNameLog.c_str(), std::ofstream::out | std::ofstream::app);
      flog<<" eviName: "<<eviName<<"\n";
      flog<<" experimentName:"<<experimentName<<"\n";
      flog<<" priorName: "<<priorName<<"\n";
      flog<<" dx: "<<dx<<"\n";
      flog<<" dtmin: "<<dtmin<<"\n";
      flog<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
      flog<<" dtmax: "<<dtmax<<"\n";
      flog<<" niter: "<<niter<<"\n";
      flog<<" maxduration "<<maxduration<<"\n";
      flog<<" landa0 "<<landa0<<"\n";
      flog<<" v "<<v<<"\n";
      flog<<" nmaxloop "<<nmaxloop<<"\n";

      flog<<" initseed "<<initseed<<"\n";
      flog<<" adaptive beta "<<aBeta<<"\n";
      flog<<AP::ClassName()<<" "<<aps<<"\n";
      flog<<" maxTime (h) "<<maxTime<<"\n";

      flog<<" samples "<<samples<<"\n";
      flog<<" nskip "<<nskip<<"\n";
      flog<<" pTjump "<<pTjump<<"\n";

      flog.flush();


      std::chrono::steady_clock::time_point startTime=std::chrono::steady_clock::now();
      double timeOpt=0;


      typename TT::myEvidence * ev= tt.run(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,pTjump,seed,eviName,flog,startTime,timeOpt);
      std::cout<<*ev;
      flog<<*ev;
      flog.close();
      std::ofstream fout;
      fout.open(eviName.c_str(), std::ofstream::out | std::ofstream::app);
      fout<<" eviName: "<<eviName<<"\n";
      fout<<" experimentName:"<<experimentName<<"\n";
      fout<<" priorName: "<<priorName<<"\n";
      fout<<" dx: "<<dx<<"\n";
      fout<<" dtmin: "<<dtmin<<"\n";
      fout<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
      fout<<" dtmax: "<<dtmax<<"\n";
      fout<<" niter: "<<niter<<"\n";
      fout<<" maxduration "<<maxduration<<"\n";
      fout<<" landa0 "<<landa0<<"\n";
      fout<<" v "<<v<<"\n";
      fout<<" nmaxloop "<<nmaxloop<<"\n";

      fout<<" initseed "<<initseed<<"\n";
      fout<<" adaptive beta "<<aBeta<<"\n";
      fout<<" maxTime (h) "<<maxTime<<"\n";

      fout<<" samples "<<samples<<"\n";
      fout<<" nskip "<<nskip<<"\n";

      fout<<*ev<<"\n";
      fout.close();
    }




}

};

/*
class TemperingCommand:public CommandBase
{
  // CommandBase interface
public:
  virtual void run(const std::__cxx11::string &line, std::ostream& logs);

  virtual std::string id() const
  {
    return "tempering";
  }
  TemperingCommand(CommandManager* cm):
    CommandBase("tempering"),
    cm_(cm){}

  ~TemperingCommand(){}
private:
  CommandManager* cm_;

};


*/



inline
std::string operator+(const std::string& one,const  std::string& two)
{
  return one+two;

}








#endif // COMMANDMANAGER

