#ifndef BASECLASS
#define BASECLASS

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>

#include <ctime>

std::istream& safeGetline(std::istream& is, std::string& t);

inline
std::string leadingZero(int i)
{
  if (i==0)
    return "00";
  else if (i<10)
    return "0"+std::to_string(i);
  else return std::to_string(i);
}

inline
std::string time_now()
{
  auto tc=std::chrono::system_clock::now();
  std::time_t rawtime=std::chrono::system_clock::to_time_t(tc);

  auto tcount=(tc.time_since_epoch().count()/1000)%1000000;


  struct std::tm * t;
  time (&rawtime);
  t = localtime (&rawtime);
  return leadingZero(t->tm_hour)+leadingZero(t->tm_min)+leadingZero(t->tm_sec)+"s"+std::to_string(tcount);

}



template<typename T>
void writeField(std::ostream& s
                , const std::string& fieldName
                , const T& value)
{
  s<<fieldName<<"\t"<<value<<"\n";
}


inline
void writeTable(std::ostream& s,
                const std::string& tableTitle,
                const std::string& xtitle,
                const std::vector<double>& x,
                const std::vector<std::vector<std::string>> ytitles,
                const std::vector<std::vector< std::vector<double>>>& ys,
                std::size_t i0=0, std::size_t iend=0)
{

  s<<std::scientific;
  s.precision(10);
  s<<tableTitle<<"\n";
  s<<xtitle<<"\t";
  for (std::size_t j=0; j<ytitles.size(); ++j)
    {
      for (auto& yt:ytitles[j])
        s<<yt<<"\t";
      s<<"\t";
    }
  s<<"\n";
  if (iend==0)
    iend=x.size();
  for (std::size_t i=i0; i<iend; ++i)
    {
      s<<x[i]<<"\t";
      for (std::size_t j=0; j<ys.size(); ++j)
        {
          if (ys[j].size()==1)
            s<<ys[j][0][i]<<"\t";
          else
            for (std::size_t k =0; k<ys[j][i].size(); ++k)
              {

                s<<ys[j][i][k]<<"\t";
              }
          s<<"\t";
        }
      s<<"\n";
    }
  s<<"\n"<<"\\\\"<<std::string(40,'-')<<"\\\\\n";
  s<<"\\\\"<<std::string(40,'-')<<"\\\\\n";
  s<<"\n";
}





bool readValue(std::string& line,
               std::istream&,
               double& val);

bool readValue(std::string& line,
               std::istream&,
               std::size_t& val);
bool readValue(std::string& line,
               std::istream&,
               std::string& val);



template<typename T>
bool readValue(std::string& line,
               std::istream& s,
               std::vector<T>& val)
{
  safeGetline(s,line);
  val.clear();
  T x;
  while (readValue(line,s,x))
    val.push_back(x);
  return !val.empty();
}


template<typename T>
bool readValue(std::string& line,
               std::istream& s,
               std::vector<std::vector<T>>& matrix)
{
  matrix.clear();
  std::vector<T> v;
  while (readValue(line,s,v))
    {
      matrix.push_back(v);
    }
  return !matrix.empty();
}

template<typename T>
bool readField(std::string& line,
               std::istream& s
               , const std::string& fieldName
               , T& value)
{
  while (line.empty()&& safeGetline(s,line)) {}
  std::stringstream ss(line);
  std::string fname;
  ss>>fname;

  if (fieldName==fname)
    {
      safeGetline(ss,line);
      return readValue(line,s,value);
    }
  else return false;
}





template<typename T>
void writeField(std::ostream& s
                , const std::string& fieldName
                , const std::vector<T>& value)
{
  s<<fieldName<<"\n";
  for (const T& x:value)
    s<<x<<"\t";
  s<<"\n";
}

template<typename T>
bool readField(std::string& line,
               std::istream& s
               , const std::string& fieldName
               ,  std::vector<T>& value)
{
  while (line.empty()&&safeGetline(s,line)) {}
  std::stringstream ss(line);
  std::string fname;
  if ((ss>>fname)&& fname==fieldName)
    {

      // safeGetline(s,line);
      return readValue(line,s,value);
    }
  else return false;
}

template<typename T>
void writeField(std::ostream& s
                , const std::string& fieldName
                , const std::vector<std::vector<T>>& matrix)
{
  s<<fieldName<<"\n";
  for (const std::vector<T>& v:matrix)
    {
      for (const T& x:v)
        s<<x<<"\t";
      s<<"\n";
    }
  s<<"\n";
}





template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v)
{
  for (T x:v)
    s<<x<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<std::vector<T>>& m)
{
  for (const std::vector<T>& v:m)
    s<<v;
  s<<"\n";
  return s;
}



template<typename T>
std::istream& operator>>(std::istream& s, std::vector<T>& v)
{
  std::string line;
  safeGetline(s,line);
  while (v.empty()  &&line.empty()&& s.good())
    safeGetline(s,line);

  T x;
  std::stringstream ss(line);
  while (ss>>x)
    v.push_back(x);
  return s;
}

template<typename T>
std::istream& operator>>(std::istream& s, std::vector<std::vector<T>>& m)
{
  std::vector<T> v;
  while((s>>v)&& !v.empty())
    {
      m.push_back(v);
      v.clear();
    }
  return s;
}






class BaseClass
{
public:
  std::string id()const
  {
    return id_;
  }

  virtual std::string myClass()const =0;
  ~BaseClass(){}

protected:
  void setId(std::string id)
  {
    id_=id;
  }

private:
  std::string id_;
};


class BaseObject: public BaseClass
{
public:


  virtual ~BaseObject(){}

  virtual BaseObject* create()const=0;

  std::ostream& write(std::ostream& s)const
  {
    s<<myClass()<<"\t"<<id()<<"\n";
    writeBody(s);
    return s;
  }



  std::string getSaveName(std::string name)
  {

    std::string filename=myClass()+"_"+name+"_"+time_now();
    setId(filename);
    return filename;
  }

  std::string save(std::string name)
  {
    std::ofstream f;
    std::string filename=getSaveName(name);
    f.open(filename.c_str());
    write(f);
    f.close();
    return filename;
  }

  void store(std::string name)
  {
    std::ofstream f;
    f.open(name.c_str());
    write(f);
    f.close();
  }




  std::string load(const std::string& name)
  {
    std::string filename=name;
    std::ifstream fi;
    fi.open(filename.c_str(),std::ios_base::in);
    if (!fi)
      {
        fi.close();
        filename=name+".txt";
        fi.open(filename.c_str(),std::ios_base::in);
        if (!fi)
          {
            fi.close();
            filename=myClass()+"_"+name+".txt";
            fi.open(filename.c_str(),std::ios_base::in);
            if (!fi)
              {
                fi.close();

                std::cerr<<"cannot load "<<name<<" file "<<filename<<" not found\n";
                return "";
              }
            else
              std::cout<<"file "<<filename<<" opened successfully\n";
          }
        else
          std::cout<<"file "<<filename<<" opened successfully\n";

      }
    else
      std::cout<<"file "<<filename<<" opened successfully\n";
    std::string line;
    if (read(line,fi))
      {
        std::cout<<myClass()<<" "<<id()<<" loaded successfully \n";
        fi.close();
        return line;
      }
    else
      {
        std::cerr<<myClass()<<" "<<id()<<" could not be loaded successfully \n";
        fi.close();
        return line;
      }

  }






  virtual std::ostream& writeBody(std::ostream& s)const=0;

  virtual void clear()=0;


  virtual std::ostream& extract(std::ostream& s,
                                const std::string& ,
                                const std::string& )const{return s;}




  bool read(std::string &line, std::istream &s)
  {
    std::string idv;
    if (readField(line,s,myClass(),idv))
      {
        setId(idv);
        safeGetline(s,line);
        if (readBody(line,s))
          {
            update();
            return true;
          }
        else
          return false;
      }
    else
      return false;

  }

  virtual bool readBody(std::string& line,std::istream& s)=0;


protected:
  virtual void update()=0;

};


class CommandManager;

class BaseAgent: public BaseObject
{
public:
  void setCommandManager(CommandManager* cm){cm_=cm;}
  CommandManager* getCommandManager()const{return cm_;}

  ~BaseAgent(){}

protected:
  CommandManager* cm_=nullptr;

};


inline
std::ostream& operator<<(std::ostream& s
                         , const BaseObject& object)
{
  object.write(s);
  return s;
}

inline
bool readValue(std::string& line,
               std::istream& s,
               BaseObject& val)
{
  return val.read(line,s);

}


inline
bool readValue(std::string& line,
               std::istream&,
               double& val)
{
  std::stringstream ss(line);
  bool o=bool(ss>>val);
  if (!o)  //check for inf nan and -nan values
    {
      auto i0=line.find_first_not_of(" \t");
      if (i0==line.npos)
        {
          safeGetline(ss,line);
          return false;
        }
      auto ie=line.find_first_of(" \t",i0);
      std::string sv=line.substr(i0,ie-i0);
      if (sv=="nan")
        {
          val=std::numeric_limits<double>::quiet_NaN();
          line=line.substr(ie);
          return true;
        }
      else if (sv=="-nan")
        {
          val=-std::numeric_limits<double>::quiet_NaN();
          line=line.substr(ie);
          return true;
        }
      else if (sv=="inf")
        {
          val=-std::numeric_limits<double>::infinity();
          line=line.substr(ie);
          return true;
        }
      else if (sv=="-inf")
        {
          val=-std::numeric_limits<double>::infinity();
          line=line.substr(ie);
          return true;
        }
      else
        return false;
    }
  else
    {
      safeGetline(ss,line);
      return true;
    }
}


inline
bool readValue(std::string& line,
               std::istream&,
               std::string& val){
  std::stringstream ss(line);
  bool o=bool(ss>>val);
  safeGetline(ss,line);
  return o;
}

inline
bool readValue(std::string& line,
               std::istream&,
               std::size_t& val){
  std::stringstream ss(line);
  bool o=bool(ss>>val);
  safeGetline(ss,line);
  return o;

}

#endif // BASECLASS

