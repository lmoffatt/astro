#ifndef MYINPUTSERIALIZER_H
#define MYINPUTSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>


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
std::istream& operator>>(std::istream& is, std::vector<T>& v)
{
  std::string line;
  if (is.peek()=='[')
    {
      std::vector<T> o;
      char ch;
      while ((is>>ch)&&(ch!='[')){}
      while (ch!=']')
        {
          std::string s;
          while ((is.get(ch))&&((ch!=']')))
            {
              s.push_back(ch);
            }
          std::stringstream ss(s);
          T e;
          while (ss>>e) o.push_back(e);

        }
      v=o;
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
std::istream& operator>>(std::istream& is, std::vector<std::vector<T>>& m)
{
  char ch;
  is>>ch;
  if (ch=='[')
    {
      std::vector<std::vector<T>> mo;
      while ((is>>ch)&&(ch!='[')){}
      while (ch!=']')
        {
          std::vector<T> o;
          while (ch!=']')
            {

              std::string s;
              while ((is.get(ch))&&((ch!=']')))
                {
                  s.push_back(ch);
                }
              std::stringstream ss(s);
              T e;
              while (ss>>e) o.push_back(e);
            }
         mo.push_back(o);
         while ((is>>ch)&&(ch!='[')&&(ch!=']')){}

        }
      m=mo;
      return is;

    }
  else
    {
      is.putback(ch);
      std::vector<T> v;
      while((is>>v)&& !v.empty())
        {
          m.push_back(v);
          v.clear();
        }
      return is;
    }
}








#endif // MYINPUTSERIALIZER_H
