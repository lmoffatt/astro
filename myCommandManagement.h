#ifndef MYCOMMANDMANAGEMENT_H
#define MYCOMMANDMANAGEMENT_H

#include <string>
#include <map>
#include <initializer_list>
#include <sstream>
#include <fstream>


#include "myTuples.h"

template<typename T>
class C
{

};

template<typename... T>
class Cs
{

};

template<template<typename> class C>
class Co
{

};


template<template<typename> class Cls,typename T>
bool load_from_file(const std::string& file_name, T& x, std::ostream& logstream)
{
  std::string filename=file_name;
  std::ifstream fi;
  fi.open(filename.c_str(),std::ios_base::in);
  if (!fi)
    {
      fi.close();
      filename=file_name+".txt";
      fi.open(filename.c_str(),std::ios_base::in);
      if (!fi)
        {
          fi.close();
          filename=Cls<T>::name()+"_"+file_name+".txt";
          fi.open(filename.c_str(),std::ios_base::in);
          if (!fi)
            {
              fi.close();
              logstream<<"cannot load "<<file_name<<" file "<<filename<<" not found\n";
              return false;
            }
          else
            logstream<<"file "<<filename<<" opened successfully\n";
        }
      else
        logstream<<"file "<<filename<<" opened successfully\n";

    }
  else
    logstream<<"file "<<filename<<" opened successfully\n";
  if (Cls<T>::read(fi,x,logstream))
    {
      logstream<<Cls<T>::name()<<" "<<file_name<<" loaded successfully \n";
      fi.close();
      return true;
    }
  else
    {
      logstream<<Cls<T>::name()<<" "<<file_name<<" not loaded successfully \n";
      fi.close();
      return false;
    }

}

template<typename Cm>
class myBaseCommand
{
public:
  virtual void run(Cm* cm,
                   const std::string& idResult,
                   const std::map<std::string,std::string>& args,
                   std::ostream& logstream)=0;
};



template <class Cm,template<typename>class Cls,class F,class R, typename... Argsout>
R run_impl_1(C<R>,Cm* cm,Co<Cls>,
             const F& f,
             std::ostream& log,
             std::map<std::string, std::string> in,
             Argsout...argsOut)
{
  return f(argsOut...);

}


template <class Cm,template<typename>class Cls,class F,class R, typename T,typename... ArgsIn,typename... Argsout>
R run_impl_1(C<R>,Cm* cm,Co<Cls>,
             const F& f,
             std::ostream& log,
             std::pair<std::string,T> a,
             std::pair<ArgsIn,std::string>... argsIn,
             const std::map<std::string, std::string>& in,
             Argsout...argsOut)
{
  auto it=in.find(a.first);
  if (it==in.end())
    {
      log<<a.first<<"=  [not listed, use default], ";
    }
  else
    {
      if (cm->template get<T>(it->second,a.second))
        {
          log<<it->first<<"="<<it->second<<" [variable], ";
        }
      else  if (cm->template load<T>(it->second,a.second))
        {
          log<<it->first<<"="<<it->second<<" [in file], ";
        }
      else
        {
          std::stringstream ss(it->second);
          if (ss>>a.second)
            {
              log<<it->first<<"="+it->second<<", ";
            }
          else
            {
              log<<it->first<<"="<<it->second<<" [error, use default],";

            }
        }
    }
  return run_impl_1(C<R>(),cm,Co<Cls>(),f,log,argsIn...,in,argsOut...,a.second);
}


template <class Cm,template<typename>class Cls,class F,class R, typename T,int... Is,typename... ArgsIn>
R run_impl_0
(C<R> ,Cm* cm,Co<Cls>,const F& f,std::ostream& ls,
 index_tuple<Is...>,
 std::tuple<std::pair<ArgsIn,std::string>...> argsIn,
 const std::map<std::string, std::string>& in)
{
  return run_impl_1(C<R>(),cm,Co<Cls>(),f,ls,std::get<Is...>(argsIn,in));
}




/*!
 * \brief wrapper class around a function
 *  \tparam Cm Command Manager class
 * \tparam F function class
 * \tparam R return type
 * \tparam Args ordered list of function arguments
 */
template <class Cm,template<typename>class Cls,class F,class R, typename... Args>
class myCommand: public myBaseCommand<Cm>
{


  /*!
   * \brief run run command on line input
   * \param cm_ the commandManager
   * \param line arguments: arg1=value1, arg2=value2 value could be text value or id
   * \param logstream logs the arguments actually used with their origin
   *
   * \return object of running the command and log text
   *
   * the arguments could be in any order. If an argument is absent, its default value is used instead. The log strean receives all the arguments actually used with their origin
   */
  R runit(Cm* cm_,Co<Cls>(),std::map<std::string, std::string> m, std::ostream& logstream)
  {

    return run_impl_0(cm_,f_,Co<Cls>(),logstream,typename make_indexes<Args...>::type(),args_,m);
  }

  myCommand(const F& f,std::pair<std::string,Args>... args):
    f_(f),
    args_{args...}{}

private:
  const F& f_;
  std::tuple<std::pair<std::string,Args>...> args_;

  static std::map<std::string, std::string>
  parse_args(const std::string& line)
  {
    std::map<std::string, std::string> out;
    // name0=value0, name1=value1,

    std::stringstream ss(line);
    std::string s;
    while (std::getline(ss,s,','))
      {
        auto npos=s.find('=');
        out[s.substr(0,npos)]=s.substr(npos+1);
      }
    return out;
  }

  // myBaseCommand interface
public:
  virtual void run(Cm *cm,
                   const std::string& idResult,
                   std::map<std::string, std::string> args,
                   std::ostream& logstream) override
  {

    R o=runit(cm,Co<Cls>(),args,logstream);
    cm->template push_back<R>(idResult,o);
  }
};







template <class Cm,template<typename>class Cls,class F, typename... Argsout>
void run_void_impl_1(Cm* cm,
                     Co<Cls>,
                     const F& f,
                     std::ostream& log,
                     std::map<std::string, std::string> in,
                     Argsout...argsOut)
{
  f(argsOut...);

}


template <class Cm,template<typename>class Cls,class F,typename T,typename... ArgsIn,typename... Argsout>
void run_void_impl_1(Cm* cm,
                     Co<Cls>,
                     const F& f,
                     std::ostream& log,
                     const std::map<std::string, std::string>& in,
                     std::pair<T, std::string> a,
                     std::pair<ArgsIn,std::string>... argsIn,
                     Argsout...argsOut)
{
  auto it=in.find(a.second);
  if (it==in.end())
    {
      log<<a.second<<"=  [not listed, use default], ";
    }
  else
    {
      std::stringstream ss(it->second);
      if (Cls<T>::read(ss,a.first,log))
        {
          log<<it->first<<"="+it->second<<", ";
        }
      else if (cm->template get<T>(it->second,a.first))
        {
          log<<it->first<<"="<<it->second<<" [variable], ";
        }
      else  if (load_from_file<Cls,T>(it->second,a.first,log))
        {
          log<<it->first<<"="<<it->second<<" [in file], ";
        }
      else
        {

              log<<it->first<<"="<<it->second<<" [error, use default],";

        }
    }
  run_void_impl_1(cm,Co<Cls>(),f,log,in,argsIn...,argsOut...,a.first );
}


template <class Cm,template<typename>class Cls,class F,typename T,typename... ArgsIn>
void run_void_impl_1(Cm* cm,
                     Co<Cls>,
                     const F& f,
                     std::ostream& log,
                     const std::map<std::string, std::string>& in,
                     std::pair<T, std::string> a,
                     std::pair<ArgsIn,std::string>... argsIn)
{
  auto it=in.find(a.second);
  if (it==in.end())
    {
      log<<a.second<<"=  [not listed, use default], ";
    }
  else
    {
      if (cm->template get<T>(it->second,a.first))
        {
          log<<it->first<<"="<<it->second<<" [variable], ";
        }
      else  if (load_from_file<Cls,T>(it->second,a.first,log))
        {
          log<<it->first<<"="<<it->second<<" [in file], ";
        }
      else
        {
          std::stringstream ss(it->second);
          if (ss>>a.second)
            {
              log<<it->first<<"="+it->second<<", ";
            }
          else
            {
              log<<it->first<<"="<<it->second<<" [error, use default],";

            }
        }
    }
  run_void_impl_1(cm,Co<Cls>(),f,log,in,argsIn...,a.first );
}



template <class Cm,template <typename T>class Cls,class F, std::size_t... Is,typename... ArgsIn>
void run_void_impl_0
(Cm* cm,
 Co<Cls>,
 const F& f,
 std::ostream& ls,
 const std::map<std::string, std::string>& in,
 std::index_sequence<Is...>,
 std::tuple<std::pair<ArgsIn,std::string>...> argsIn
 )
{
  run_void_impl_1(cm,Co<Cls>(),f,ls,in,std::get<Is>(argsIn)...);
}







/*!
 * \brief wrapper class around a function
 *  \tparam Cm Command Manager class
 * \tparam F function class
 * \tparam R return type
 * \tparam Args ordered list of function arguments
 */
template <class Cm,template<typename>class Cls,class F,typename... Args>
class myCommand<Cm,Cls,F,void,Args...>: public myBaseCommand<Cm>
{
public:

  /*!
   * \brief run run command on line input
   * \param cm_ the commandManager
   * \param line arguments: arg1=value1, arg2=value2 value could be text value or id
   * \param logstream logs the arguments actually used with their origin
   *
   * \return object of running the command and log text
   *
   * the arguments could be in any order. If an argument is absent, its default value is used instead. The log strean receives all the arguments actually used with their origin
   */
  void runit(Cm* cm_,const std::map<std::string, std::string>& m, std::ostream& logstream)
  {
    run_void_impl_0(cm_,
                    Co<Cls>(),
                    f_,
                    logstream,
                    m,
                    std::index_sequence_for<Args...>(),
                    args_);
  }

  myCommand(const F& f,std::pair<Args,std::string>... args):
    f_(f),
    args_{args...}{}

private:
  const F& f_;
  std::tuple<std::pair<Args,std::string>...> args_;


  // myBaseCommand interface
public:
  virtual void  run(Cm* cm,
                    const std::string& idResult,
                    const std::map<std::string,std::string>& args,
                    std::ostream& logstream)
  {

    runit(cm,args,logstream);
  }
};






template <class Cm,template <typename>class Cls,class R,class F, typename... Args>
myCommand<Cm,Cls,F,R,Args...>* make_Command
(C<Cm>,C<R>,Co<Cls>,const F& f,std::pair<Args,std::string>... args)
{
  return new myCommand<Cm,Cls,F,R,Args...>(f,args...);
}



template<typename T, typename...Ts>
bool get_impl_1(const std::string& id, T& x)
{
  return false;
}


template<typename T, typename K, typename...Ts>
bool get_impl_1(const std::string& id, T& x,
                std::map<std::string,K>,
                std::map<std::string,Ts>... m)
{
  return get_impl_1(id,x,m...);
}

template<typename T, typename...Ts>
bool get_impl_1(const std::string& id, T& x,
                std::map<std::string,T> mymap,
                std::map<std::string,Ts>... m)
{
  auto it=mymap.find(id);
  if (it!=mymap.end())
    {
      x=it->second;
      return true;
    }
  else return false;
}

template<typename T, std::size_t...Is, typename...Ts>
bool get_impl_0(const std::string& id, T& x,
                std::index_sequence<Is...>,
                std::tuple<std::map<std::string,Ts>...> m)
{
  return get_impl_1(id,x,std::get<Is>(m)...);
}


template<typename T, typename...Ts>
bool get_map(const std::string& id, T& x, std::tuple<std::map<std::string,Ts>...> m)
{
  return get_impl_0(id,x,std::index_sequence_for<Ts...>(),m);
}




/*!
 * \tparam Cls
 * Cls<T>::name()
 * Cls<T>::read(std::istream&,T& x,std::ostream& logstream)
 */
template<template<typename> class Cls,class...Ts>
class myCommandManager
{
public:
  template<typename T>
  using myCls=Cls<T>;

  static std::string ClassName(){return "CommandManager";}

  typedef myBaseCommand<myCommandManager> Command;
  std::tuple<std::string,std::string, std::map<std::string,std::string>>
  line_to_CommandName_Arg_Result(const std::string& line)
  {
    auto n0=line.find_first_of('(');
    auto n1=line.find_first_of(')');
    auto n2=line.find_first_of('>',n1);
    std::string commandName=line.substr(0,n0);
    std::string argsList=line.substr(n0+1,n1-n0-1);
    std::string resultName=line.substr(n2+1);
    std::map<std::string, std::string> out;
    // name0=value0, name1=value1,

    std::stringstream ss(argsList);
    std::string s;
    while (std::getline(ss,s,','))
      {
        auto n0=s.find_first_not_of(" ");
        auto npos=s.find('=');
        out[s.substr(n0,npos-n0)]=s.substr(npos+1);
      }


    return {commandName,resultName,out};
  }


  void execute(std::string line, std::ostream& logstream)
  {
    if (!line.empty())
      {
        std::stringstream ss(line);
        auto o=line_to_CommandName_Arg_Result(line);
        if (!std::get<0>(o).empty())
          {
            Command* c=command(std::get<0>(o));
            if (c!=0)
              {
                logstream<<line<<"\n";
                c->run(this,std::get<1>(o),std::get<2>(o),logstream);
              }
            else
              logstream<<"error in"<<line<<" "<<std::get<0>(o)<<" is not a command";
          }
      }
  }


  //CommandBase* Command(std::string commandName);

  template<typename T>
  bool get(const std::string& id, T& x)
  {
    return get_map(id,x,data_);
  }




  template <typename T>
  void push_back(const std::string& id, T* x)
  {
    auto m=std::get<typename std::map<std::string,T*>>(data_);
    m[id]=x;
  }
  void push_command(const std::string& id, Command * cmd)
  {
    cmds_[id]=cmd;
  }

  Command* command(const std::string& id)
  {
    auto it=cmds_.find(id);
    if (it!=cmds_.end())
      return it->second;
    else
      return nullptr;

  }


private:
  std::tuple<std::map <std::string, Ts*>...> data_;

  std::map<std::string,myBaseCommand<myCommandManager>*> cmds_;

};




#endif // MYCOMMANDMANAGEMENT_H
