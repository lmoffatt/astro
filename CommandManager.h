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

class MyModel
{
public:
  M_Matrix<double> sample(std::mt19937_64& mt)const
  {
    Parameters sample=CL_->getPrior().randomSample(mt,1);
    return M_Matrix<double> (1,sample.size(),sample.trMeans());
  }

  const Parameters& getPrior()const
  {
    return CL_->getPrior();
  }


  mcmc_prior<double> prior(const MyData&  , M_Matrix<double> param)const
  {
    mcmc_prior<double> out;
    out.param=param;
    auto& p=CL_->getPrior();
    //    std::size_t npar=param.size();

    Parameters Pa=p.toParameters(param.toVector());

    out.logPriorLikelihood=p.logProb(Pa);
    out.D_prior.H=M_Matrix<double>(M_Matrix<double>::SYMMETRIC,p.getInvCovariance());
    M_Matrix<double> d(param.nrows(),param.ncols(),p.trMeans());
    d=param - d;


    out.D_prior.G=d*out.D_prior.H;

    return out;

  }
  M_Matrix<double> f(const MyData& e, M_Matrix<double> param, std::pair<std::vector<double>,std::vector<std::size_t>>& dts) const
  {
    auto ff=CL_->f(e.myExperiment(),param.toVector(),dts);
    std::size_t nrows= ff.size();
    if (nrows>0)
      {
        std::size_t ncols=ff[0].size();
        M_Matrix<double> out(nrows,ncols);
        for (std::size_t i=0; i<nrows; ++i)
          for (std::size_t j=0; j<ncols; ++j)
            out(i,j)=ff[i][j];
        return out;
      }
    else return {};

  }

  M_Matrix<double> f(const MyData& e, M_Matrix<double> param,
                     double& dtmin,std::size_t& nper10,double& dtmax,
                     std::size_t dts_max,
                     std::pair<std::vector<double>,std::vector<std::size_t>>& dts)
  const{
    std::tuple<double,std::size_t,double> dtmin_n_dtmax{dtmin,nper10,dtmax};
    auto ff=CL_->f(e.myExperiment(),param.toVector(),dtmin_n_dtmax,dts_max,dts);
    dtmin=std::get<0>(dtmin_n_dtmax);
    nper10=std::get<1>(dtmin_n_dtmax);
    dtmax=std::get<2>(dtmin_n_dtmax);

    std::size_t nrows= ff.size();
    if (nrows>0)
      {
        std::size_t ncols=ff[0].size();
        M_Matrix<double> out(nrows,ncols);
        for (std::size_t i=0; i<nrows; ++i)
          for (std::size_t j=0; j<ncols; ++j)
            out(i,j)=ff[i][j];
        return out;
      }
    else return {};

  }



  M_Matrix<double> logLanda(const MyData& data, M_Matrix<double> param, std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const{
    M_Matrix<double> out=f(data,param,dts);
    out=out.apply([](double x){return log10_guard(x);});
    return out;

  }

  std::string id()const {return CL_->getPrior().id();}




  Parameters getParameter(const M_Matrix<double>& par)const
  {
    return CL_->getPrior().toParameters(par.toVector());
  }

  CortexSimulation getSimulation(const Experiment* e,const Parameters& par,const std::pair<std::vector<double>,std::vector<std::size_t>>& dts)const
  {
    return CL_->simulate(e,par,dts);
  }

  CortexSimulation getSimulation(const Experiment* e,const M_Matrix<double>& par,std::pair<std::vector<double>, std::vector<std::size_t>> dts)const
  {
    auto ec= const_cast<Experiment*> (e);
    ec->setSimulation();
    auto out=CL_->simulate(e,getParameter(par),dts);
    ec->setMeasure();
    return out;
  }

  const CortexLikelihood& getLikelihood()const
  {
    return *CL_;
  }


  MyModel(CortexLikelihood* CL):CL_(CL){}
private:
  const CortexLikelihood* CL_;
};



template <class Ad>
using TI=
Thermodynamic_Integration_mcmc<Ad,
MyData,MyModel,Poisson_DLikelihood<MyData,MyModel>,LM_MultivariateGaussian<double>,Landa,LevenbergMarquardt_step>;

template <class Ad, class MyData, class Lik>
using TT=
Template_Tempering_mcmc<Ad,
MyData,MyModel,Lik,LM_MultivariateGaussian<typename Lik::E>,Landa,LevenbergMarquardt_step> ;





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
auto read_from_stream(std::istream& is,std::ostream& /*logstream*/,T*& x)
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
  static void write(const T& x,std::ostream* f, std::ostream*/* log*/)
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

inline
std::pair<std::size_t, std::string>
extract_Seed(const std::string& s)
{
   auto last=s.find("_state");
   auto first=s.find_last_of('_',last-1);
   auto v= s.substr(first+1,last-first-1);
   auto val=std::stoull(v);
   auto eviName=s.substr(0,last);
   return {val,eviName};
}

template<class Cm>
struct Tempering
{
  void operator ()(Cm* cm_,
                   std::string eviName,
                   std::string experimentName,
                   std::string priorName,
                   std::string state_file,
                   double dtmin0,
                   double dtmin,
                   double dtmax,
                   double dx,
                   bool adapt_dt,
                   double maxlogError,
                   double dtinf,
                   bool CrankNicholson,
                   double f_maxlogError,
                   std::size_t maxloop,
                   bool UseDerivative,
                   double tequilibrio,
                   double  maxduration,
                   double landa0,
                   double v,
                   double pTjump,
                   double slogL_max,
                   std::size_t ndts_max,
                   std::size_t nmaxloop,
                   std::mt19937_64::result_type initseed,
                   std::size_t nPoints_per_decade,
                   std::size_t niter,
                   std::size_t N_betasInit,
                   double beta_min,
                   std::size_t N_beta_2,
                   double beta_infimo,
                   const std::string Landa_algorithm,
                   double targetProb,
                   M_Matrix<Landa> aps,
                   std::vector<std::vector<double>> apsPar,
                   std::size_t gainMoment,
                   bool unInformativePriorPar,
                   double maxTime,
                   std::size_t samples,
                   std::size_t nskip,
                   std::size_t nAdapt,

                   std::size_t maxSimFileSize,
                   bool does_stdout
                   , std::ostream *logs)const
  {

    // evidence evi experiment1 paramters_10
    //optimize opt experiment1 parameters_10 parameters_10 10 100
    typedef Landa AP;


    Master_Adaptive_Beta_New
        aBeta(N_betasInit,beta_min,N_beta_2,beta_infimo);
    std::random_device rd;

    std::mt19937_64::result_type seed;
    if (does_stdout)
      {
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
      }

    bool isContinuation=!state_file.empty();

    if (!isContinuation)
      {
        if (initseed==0)
          {
            seed=rd();

            eviName+="_"+time_now()+"_"+std::to_string(seed);
            *logs<<"\n random seed =\n"<<seed<<"\n";
          }
        else
          {
            seed=initseed;

            eviName+="_"+time_now()+"_"+std::to_string(seed);

            *logs<<"\n provided seed =\n"<<seed<<"\n";

          }
      }
    else
      {
        auto o=extract_Seed(state_file);
        seed=o.first;
        eviName="R_"+o.second;
        *logs<<"\n seed of previous run=\n"<<seed<<"\n";
      }



    Experiment esim;



    auto e=Experiment::loadList(experimentName,logs);

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

        CortexPoisonLikelihood* CL;
        if (!CrankNicholson)
          {
            if (!adapt_dt)
              CL=new CortexPoisonLikelihood
                  (eviName+"_lik",prior,dx,dtmin0,dtmin,
                   nPoints_per_decade,dtmax,tequilibrio);
            else
              CL=new CortexPoisonLikelihood(eviName+"_lik",prior,dx,dtmin0,dtmin,nPoints_per_decade,dtmax,tequilibrio,maxlogError,dtinf);
          }
        else
          {
            if (!adapt_dt)
              CL=new CortexPoisonLikelihood
                  (eviName+"_lik",prior,dx,dtmin0,dtmin,nPoints_per_decade,
                   dtmax,tequilibrio,f_maxlogError,maxloop);
            else
              CL=new CortexPoisonLikelihood
                  (eviName+"_lik",prior,dx,dtmin0,dtmin,nPoints_per_decade,
                   dtmax,tequilibrio,maxlogError,f_maxlogError,dtinf,maxloop,UseDerivative);
          }
        MyModel m(CL);
        if (e.size()==1)
          {
            MyData d(CL,&e[0]);
            LevenbergMarquardt_step
                <
                MyData,MyModel,Poisson_DLikelihood<MyData,MyModel>
                ,LM_MultivariateGaussian<double>,Landa
                > LMLik;
            Poisson_DLikelihood<MyData,MyModel> DLik;



            std::string eviNameLog0=eviName+"_log.txt";
            std::string eviNameLog=eviNameLog0+".0";
            std::ofstream flog;

            flog.open(eviNameLog.c_str(), std::ofstream::out | std::ofstream::app);
            flog<<" eviName: "<<eviName<<"\n";
            flog<<" experimentName:"<<experimentName<<"\n";
            flog<<" priorName: "<<priorName<<"\n";
            flog<<" dx: "<<dx<<"\n";
            flog<<" dtmin: "<<dtmin<<"\n";
            flog<<" dtmax: "<<dtmax<<"\n";
            flog<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
            flog<<" adapt_dt: "<<adapt_dt<<"\n";
            flog<<" maxlogError: "<<maxlogError<<"\n";
            flog<<" dtinf: "<<dtinf<<"\n";
            flog<<" CrankNicholson: "<<CrankNicholson<<"\n";
            flog<<" f_maxlogError: "<<f_maxlogError<<"\n";
            flog<<" maxloop: "<<maxloop<<"\n";
            flog<<" UseDerivative: "<<UseDerivative<<"\n";
            flog<<" slogL_max: "<<slogL_max<<"\n";
            flog<<" ndts_max: "<<ndts_max<<"\n";
            flog<<" nmaxloop: "<<nmaxloop<<"\n";
            flog<<" dtinf: "<<dtinf<<"\n";
            flog<<" N_betasInit: "<<N_betasInit<<"\n";
            flog<<" beta_min: "<<beta_min<<"\n";
            flog<<" N_beta_2: "<<N_beta_2<<"\n";
            flog<<" beta_infimo: "<<beta_infimo<<"\n";
            flog<<" Landa_algorithm: "<<Landa_algorithm<<"\n";
            flog<<" gainMoment: "<<gainMoment<<"\n";
            flog<<" unInformativePriorPar: "<<unInformativePriorPar<<"\n";
            flog<<" targetProb: "<<targetProb<<"\n";
            flog<<" maxSimFileSize: "<<maxSimFileSize<<"\n";
            flog<<" beta_min: "<<beta_min<<"\n";
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
            flog<<" nAdapt "<<nAdapt<<"\n";
            flog<<" pTjump "<<pTjump<<"\n";




            flog.close();


            std::chrono::steady_clock::time_point startTime=std::chrono::steady_clock::now();
            double timeOpt=0;

            if (Landa_algorithm.find("Probability")!=std::string::npos)
              {
                typedef
                Adaptive_probability<One<AP>,TargetProb,AP> Ad;
                TargetProb t(targetProb);
                Ad landa(t,aps);
                if (does_stdout)
                  std::cout<<landa;

                Metropolis_Hastings_mcmc<Ad,
                    MyData,MyModel,Poisson_DLikelihood<MyData,MyModel>,LM_MultivariateGaussian<double>,Landa> mcmc;
                TT<Ad,MyData,Poisson_DLikelihood<MyData,MyModel> > tt;

                tt.run<double>(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,nAdapt,pTjump,slogL_max,ndts_max,seed,eviName,eviNameLog0,startTime,timeOpt,maxSimFileSize,does_stdout,state_file );
                flog.close();

              }
            else if (Landa_algorithm.find("Velocity")!=std::string::npos)
              {
                typedef
                Adaptive_probability<AP::ExpectedVelocity,PascalProb,AP> Ad;
                Ad landa(aps);
                if (does_stdout)
                  std::cout<<landa;

                Metropolis_Hastings_mcmc<Ad,
                    MyData,MyModel,Poisson_DLikelihood<MyData,MyModel>,LM_MultivariateGaussian<double>,Landa> mcmc;
                TT<Ad, MyData,Poisson_DLikelihood<MyData,MyModel>> tt;

                tt.run<double>(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,nAdapt,pTjump,slogL_max,ndts_max,seed,eviName,eviNameLog0,startTime,timeOpt,maxSimFileSize,does_stdout,state_file );
                flog.close();

              }
            else
              {
                Adaptive_parameterized<AP> landa(aps,apsPar,gainMoment);

                if (does_stdout)
                  std::cout<<landa;

                Metropolis_Hastings_mcmc<Adaptive_parameterized<AP>,
                    MyData,MyModel,Poisson_DLikelihood<MyData,MyModel>,LM_MultivariateGaussian<double>,Landa> mcmc;
                TT<Adaptive_parameterized<AP>,MyData,Poisson_DLikelihood<MyData,MyModel>> tt;

                tt.run<double>(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,nAdapt,pTjump,slogL_max,ndts_max,seed,eviName,eviNameLog0,startTime,timeOpt,maxSimFileSize,does_stdout,state_file );
                flog.close();
              }


          }

        else  //Random effects

          {
            std::vector<MyData> d(e.size());
            for (std::size_t i=0; i<e.size(); ++i) d[i]= MyData(CL,&e[i]);
            typedef Random_Effects_Likelihood<Poisson_DLikelihood<MyData,MyModel>> RE;
            LevenbergMarquardt_step<std::vector<MyData>,MyModel,RE,LM_MultivariateGaussian<typename RE::E>,Landa> LMLik;
            RE DLik;


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

            std::string eviNameLog0=eviName+"_log.txt";
            std::string eviNameLog=eviNameLog0+".0";
            std::ofstream flog;

            flog.open(eviNameLog.c_str(), std::ofstream::out | std::ofstream::app);
            flog<<" eviName: "<<eviName<<"\n";
            flog<<" experimentName:"<<experimentName<<"\n";
            flog<<" priorName: "<<priorName<<"\n";
            flog<<" dx: "<<dx<<"\n";
            flog<<" dtmin: "<<dtmin<<"\n";
            flog<<" dtmax: "<<dtmax<<"\n";
            flog<<" nPoints_per_decade: "<<nPoints_per_decade<<"\n";
            flog<<" adapt_dt: "<<adapt_dt<<"\n";
            flog<<" maxlogError: "<<maxlogError<<"\n";
            flog<<" dtinf: "<<dtinf<<"\n";
            flog<<" CrankNicholson: "<<CrankNicholson<<"\n";
            flog<<" f_maxlogError: "<<f_maxlogError<<"\n";
            flog<<" maxloop: "<<maxloop<<"\n";
            flog<<" UseDerivative: "<<UseDerivative<<"\n";
            flog<<" slogL_max: "<<slogL_max<<"\n";
            flog<<" ndts_max: "<<ndts_max<<"\n";
            flog<<" nmaxloop: "<<nmaxloop<<"\n";
            flog<<" dtinf: "<<dtinf<<"\n";
            flog<<" N_betasInit: "<<N_betasInit<<"\n";
            flog<<" beta_min: "<<beta_min<<"\n";
            flog<<" N_beta_2: "<<N_beta_2<<"\n";
            flog<<" beta_infimo: "<<beta_infimo<<"\n";
            flog<<" Landa_algorithm: "<<Landa_algorithm<<"\n";
            flog<<" gainMoment: "<<gainMoment<<"\n";
            flog<<" unInformativePriorPar: "<<unInformativePriorPar<<"\n";
            flog<<" targetProb: "<<targetProb<<"\n";
            flog<<" maxSimFileSize: "<<maxSimFileSize<<"\n";
            flog<<" beta_min: "<<beta_min<<"\n";
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
            flog<<" nAdapt "<<nAdapt<<"\n";
            flog<<" pTjump "<<pTjump<<"\n";




            flog.close();


            std::chrono::steady_clock::time_point startTime=std::chrono::steady_clock::now();
            double timeOpt=0;
            typedef M_Matrix<M_Matrix<double>> E;

            if (Landa_algorithm.find("Probability")!=std::string::npos)
              {
                typedef
                Adaptive_probability<One<AP>,TargetProb,AP> Ad;


                TargetProb t(targetProb);
                Ad landa(t,aps);
                if (does_stdout)
                  std::cout<<landa;

                Metropolis_Hastings_mcmc<Ad,
                    std::vector<MyData>,MyModel,RE,
                    LM_MultivariateGaussian<typename RE::E>,Landa> mcmc;
                TT<Ad,std::vector<MyData>,RE> tt;

                tt.run<E>(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,nAdapt,pTjump,slogL_max,ndts_max,seed,eviName,eviNameLog0,startTime,timeOpt,maxSimFileSize,does_stdout,state_file );
                flog.close();

              }
            else if (Landa_algorithm.find("Velocity")!=std::string::npos)
              {
                typedef
                Adaptive_probability<AP::ExpectedVelocity,PascalProb,AP> Ad;
                Ad landa(aps);
                if (does_stdout)
                  std::cout<<landa;

                Metropolis_Hastings_mcmc<Ad,
                    std::vector<MyData>,MyModel,RE,
                    LM_MultivariateGaussian<typename RE::E>,Landa> mcmc;
                TT<Ad,std::vector<MyData>,RE> tt;

                tt.run<E>(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,nAdapt,pTjump,slogL_max,ndts_max,seed,eviName,eviNameLog0,startTime,timeOpt,maxSimFileSize,does_stdout,state_file );

                flog.close();

              }
            else
              {
                Adaptive_parameterized<AP> landa(aps,apsPar,gainMoment);

                if (does_stdout)
                  std::cout<<landa;
                Metropolis_Hastings_mcmc<Adaptive_parameterized<AP>,
                    std::vector<MyData>,MyModel,RE,LM_MultivariateGaussian<typename RE::E>,Landa> mcmc;
                TT<Adaptive_parameterized<AP>,std::vector<MyData>,RE> tt;

                tt.run<E>(mcmc,LMLik,DLik,m,d,landa,aBeta,maxTime,samples,nskip,nAdapt,pTjump,slogL_max,ndts_max,seed,eviName,eviNameLog0,startTime,timeOpt,maxSimFileSize,does_stdout,state_file );


                flog.close();
              }


          }



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



//inline
//std::string operator+(const std::string& one,const  std::string& two)
//{
//  return one+two;

//}








#endif // COMMANDMANAGER

