#ifndef MYOUTPUTSERIALIZER_H
#define MYOUTPUTSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other)
{
  os<<other.first<<","<<other.second;
  return os;
}

inline
std::ostream& operator<<(
    std::ostream& s,const std::vector< std::vector< double> >& matrix)
{
    s<<"\n";
    for (std::size_t i=0; i<matrix.size();++i)
    {
        for (std::size_t j=0; j<matrix[i].size();j++)
            s<<matrix[i][j]<<"\t";
        s<<"\n";
    }
    return s;
}

inline
std::ostream& operator<<(
    std::ostream& s,const std::vector<  double> & aVector){
    for (std::size_t j=0; j<aVector.size();j++)
        s<<aVector[j]<<"\t";
    s<<"\n";
return s;
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

template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::map<K,T>& v)
{
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::multimap<K,T>& v)
{
  for (auto& it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::set<T>& v)
{
  for (auto& it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::multiset<T>& v)
{
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}




#endif // MYOUTPUTSERIALIZER_H
