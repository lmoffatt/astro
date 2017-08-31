#ifndef MYINPUTSERIALIZER_H
#define MYINPUTSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "myTuples.h"

inline
std::istream &safeGetline(std::istream &is, std::string &t)
{
  is.clear();
  std::getline(is,t);
  auto it=t.find('\r');
  if (it!=t.npos)
    t.erase(it);
  return is;
}

template<typename T>
auto operator>>(std::istream& is, T& v)
->decltype(v.read(std::declval<std::string&>(),is,std::declval<std::ostream&>()),is)
{
  std::string s;
  if (!v.read(s,is, std::cerr))
    {
      is.setstate(std::ios::failbit);

    }
  return is;

}

template<typename T>
auto operator>>(std::istream& is, T*& v)
->decltype(v->read(std::declval<std::string&>(),is,std::declval<std::ostream&>() ),is)
{
  std::string s;
  if (v!=nullptr)
    {
      v->read(s,is,std::cerr);
    }
  else
    {
      is.setstate(std::ios::failbit);
    }
  return is;

}




template<typename T
         ,typename std::enable_if<!std::is_pointer<T>::value,int>::type = 0>
std::istream& operator>>(std::istream& is, std::vector<T>& v)
{
  std::string line;
  char c;
  is>>c;
  if (c=='[')
    {
      v.clear();
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              is.putback(c);
              T e;
              if(is>>e)
                 v.push_back(e);
              else
                 break;
            }
        }
      return is;
    }
  else
    {
      safeGetline(is,line);
      while (v.empty()  &&line.empty()&& is.good())
        safeGetline(is,line);

      T x;
      std::stringstream ss(line);
      while (ss>>x)
        v.push_back(x);
      return is;
    }
}
template<typename T>

std::istream& operator>>(std::istream& is, std::vector<T*>& v)
{
  std::string line;
  char c;
  is>>c;
  if (c=='[')
    {
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              T* e=new T{};
              is.putback(c);
              if(is>>*e)
                 v.push_back(e);
              else
                 break;
            }
        }
      return is;
    }
  else
    {
      safeGetline(is,line);
      while (v.empty()  &&line.empty()&& is.good())
        safeGetline(is,line);

      T* x=new T;
      std::stringstream ss(line);
      while (ss>>*x)
        {
          v.push_back(x);
          x=new T;
        }
      delete x;
      return is;
    }
}


template<typename T>
std::istream& operator>>(std::istream& is, std::vector<std::vector<T>>& m)
{
  char c;
  is>>c;
  if (c=='[')
    {
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              is.putback(c);
              std::vector<T> e;
              if(is>>e)
                 m.push_back(e);
              else
                 break;
            }
        }
      return is;
    }
  else
    {
      is.putback(c);
      std::vector<T> v;
      while((is>>v)&& !v.empty())
        {
          m.push_back(v);
          v.clear();
        }
      return is;
    }
}




template<typename K,typename T>
std::istream& operator>>(std::istream& is, std::map<K,T>& v)
{ char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              is.putback(c);
              std::pair<K,T> e;
              if(is>>e)
                 v.insert(e);
              else
                 break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
      return is;
    }
}

template<typename T1, typename T2>
std::istream& operator>>(std::istream& os,std::pair<T1,T2>& other)
{
  char ch;
  os>>other.first>>ch>>other.second;
  return os;
}

template<typename K,typename T>
std::istream& operator>>(std::istream& is,  std::multimap<K,T>& v)
{
  char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              std::pair<K,T> e;
              is.putback(c);
              if(is>>e)
                 v.insert(e);
              else
                 break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
    }
}
template<typename T>
std::istream& operator>>(std::istream& is, std::multiset<T>& v)
{
  char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              T e;
              is.putback(c);
              if(is>>e)
                 v.insert(e);
              else
                 break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
      return is;
    }
  }


template<typename Last>
std::istream& get_impl(std::istream& is, Last& last)
{
  is>>last;
  return is;
}

template<typename First, typename ... Rest>
std::istream& get_impl(std::istream& is,  First& first, Rest&...rest)
{
  get_impl(is,first);
  get_impl(is,rest...);
  return is;
}

template< int ... Indexes, typename ... Args>
std::istream& get_helper(std::istream& is, index_tuple<Indexes...>, std::tuple<Args...>& tup)
{
  get_impl( is, std::get<Indexes>(tup)...);
  return is;
}


template<typename ...Args>
std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu)
{
  return get_helper(is,typename make_indexes<Args...>::type(),
                    tu);
}



#endif // MYINPUTSERIALIZER_H
