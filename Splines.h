#ifndef SPLINES
#define SPLINES

#include "MatrixInverse.h"


#include <vector>
#include <map>
#include <cmath>
#include <iostream>



class Spline
{

  Spline(std::vector<double> x, std::vector<double> y);


  double eval( double x)const
  {
    std::size_t i; double t,s;
    if (get_index(x,i,t,s))
      {
        --i;
        return i_eval(i,t,s);
      }
    else if (i==0)
      return lower_default_;
    else
      return upper_default_;
  }

  double deval( double x)const
  {
    std::size_t i; double t,s;
    if (get_index(x,i,t,s))
      {
        --i;
        return i_deval(i,t,s);
      }
    else if (i==0)
      return lower_default_;
    else
      return upper_default_;
  }
  double d2eval( double x)const
  {
    std::size_t i; double t,s;
    if (get_index(x,i,t,s))
      {
        --i;
        return i_d2eval(i,t);
      }
    else if (i==0)
      return lower_default_;
    else
      return upper_default_;
  }

  std::vector<double> eval( const std::vector<double>& x)const
  {
    std::vector<double> o(x.size());
    for (std::size_t i=0; i<o.size(); ++i)
      o[i]=eval(x[i]);
    return o;
  }
  std::vector<double> deval( const std::vector<double>& x)const
  {
    std::vector<double> o(x.size());
    for (std::size_t i=0; i<o.size(); ++i)
      o[i]=deval(x[i]);
    return o;
  }
  std::vector<double> d2eval( const std::vector<double>& x)const
  {
    std::vector<double> o(x.size());
    for (std::size_t i=0; i<o.size(); ++i)
      o[i]=d2eval(x[i]);
    return o;
  }


private:
  bool get_index(double x,std::size_t& i, double& t, double& s)const
  {
    auto it=x_map_.lower_bound(x);
    if (it!=x_map_.end())
      {
        i=it->second;
        if (i>0)
          {
            t=(x-x_[i-1])/dx_[i-1];
            s=1-t;
            return true;
          }
        else if (x==x_[i])
          {
            i=i+1;
            t=0;
            s=1;
            return true;
          }
        else return false;
      }
    else
      return false;
  }
  double i_eval(std::size_t i, double t, double s)const
  {
    return s*y_[i]+t*y_[i+1]+t*s*(a_[i]*s+b_[i]*t);
  }
  double i_deval(std::size_t i, double t, double s)const
  {
    return (y_[i+1]-y_[i])/dx_[i]+(s-t)*(a_[i]*s+b_[i]*t)/dx_[i]+t*s*(b_[i]-a_[i])/dx_[i];

  }
  double i_d2eval(std::size_t i, double t)const
  {
    return 2*(b_[i]-2*a_[i]+(a_[i]-b_[i])*3*t)/dx_[i]/dx_[i];

  }

  std::vector<std::vector<double>> getIndexMatix(const std::vector<double> x)
  {
    auto n=x.size();
    std::vector<std::vector<double >> o(n,std::vector<double>(n,0.0));

    for (std::size_t i=1; i<n-1; ++i)
      {
        o[i][i-1]=1.0/(x[i]-x[i-1]);
        o[i][i+1]=1.0/(x[i+1]-x[i]);
        o[i][i]=2.0*(o[i][i-1]+o[i][i+1]);
      }
    o[0][0]=1.0/std::pow(x[1]-x[0],2);
    o[0][2]=-1.0/std::pow(x[2]-x[1],2);
    o[0][1]=o[0][0]+o[0][2];
    o[n-1][n-1]=-1.0/std::pow(x[n-1]-x[n-2],2);
    o[n-1][n-3]=1.0/std::pow(x[n-2]-x[n-3],2);
    o[n-1][n-2]=o[n-1][n-1]+o[n-1][n-3];
    return o;
  }

  std::vector<double> getYvector(const std::vector<double>& x, const std::vector<double>& y)
  {
    auto n=x.size();
    if (y.size()!=x.size())
      {
        std::cerr<<" unequal vector sizes \n";
        return {};
      }
    std::vector<double> o(n);
    for (std::size_t i=1;i<n-1; ++i)
      {
        o[i]=3.0*((y[i]-y[i-1])/std::pow(x[i]-x[i-1],2)+(y[i+1]-y[i])/std::pow(x[i+1]-x[i],2));
      }
    o[0]=2*((y[1]-y[0])/std::pow(x[1]-x[0],3)-(y[2]-y[1])/std::pow(x[2]-x[1],3));
    o[n-1]=2*((y[n-1]-y[n-2])/std::pow(x[n-1]-x[n-2],3)-(y[n]-y[n-1])/std::pow(x[n]-x[n-1],3));
  }


  std::vector<double> getKs(const std::vector<double>& x, const std::vector<double>& y)
  {
    auto ys=getYvector(x,y);
    auto A=getIndexMatix(x);
    auto Ainv=inv(A);
    return mult(Ainv,ys);

  }



  std::map<double,std::size_t> x_map_;
  std::size_t n_knots_;

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> k_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> dx_;
  double lower_default_;
  double upper_default_;
};



class MSpline
{
public:
  MSpline(const std::vector<double>& x, const std::vector<std::vector<double>>& y):
    x_map_()
  ,n_knots_(x.size()-1)
  ,n_cols_(y[0].size())
  , x_(x)
  , dx_(std::vector<double>(n_knots_))
  , yM_(y)
  , kM_(getKs(x,y))
  , aM_(std::vector<std::vector<double>>(n_knots_,std::vector<double>(n_cols_)))
  , bM_(std::vector<std::vector<double>>(n_knots_,std::vector<double>(n_cols_)))

  {

    for (std::size_t i=0; i<n_knots_; ++i)

      {
        dx_[i]=x_[i+1]-x_[i];
        for (std::size_t j=0; j<n_cols_; ++j)
          {
            aM_[i][j]= kM_[i][j]*dx_[i]
                      -yM_[i+1][j]+yM_[i][j];
            bM_[i][j]=-kM_[i+1][j]*dx_[i]
                      +yM_[i+1][j]-yM_[i][j];
          }
      }
    for (std::size_t i=0; i<x.size(); ++i)
      x_map_[x[i]]=i;
  }


  std::vector<double> eval( double x)const
  {
    std::vector<double> o(n_cols_);
    std::size_t i; double t,s;
    if (get_index(x,i,t,s))
      {
        --i;
        for (std::size_t j=0; j<n_cols_; ++j)
          o[j]= i_eval(i,j,t,s);
        return o;
      }
    else if (i==0)
      return std::vector<double>(n_cols_,lower_default_);
    else
      return std::vector<double>(n_cols_,upper_default_);
  }

  std::vector<double> deval( double x)const
  {
    std::vector<double> o(n_cols_);
    std::size_t i; double t,s;
    if (get_index(x,i,t,s))
      {
        --i;
        for (std::size_t j=0; j<n_cols_; ++j)
          o[j]= i_deval(i,j,t,s);
      }
    else if (i==0)
      return std::vector<double>(n_cols_,lower_default_);
    else
      return std::vector<double>(n_cols_,upper_default_);
  }
  std::vector<double> d2eval( double x)const
  {
    std::vector<double> o(n_cols_);
    std::size_t i; double t,s;
    if (get_index(x,i,t,s))
      {
        --i;
        for (std::size_t j=0; j<n_cols_; ++j)
          o[j]= i_d2eval(i,j,t);
      }
    else if (i==0)
      return std::vector<double>(n_cols_,lower_default_);
    else
      return std::vector<double>(n_cols_,upper_default_);
  }

  std::vector<std::vector<double>> eval(const std::vector<double> x)
  {
    std::vector<std::vector<double>> o(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      o[i]=eval(x[i]);
    return o;
  }


  std::vector<double> dEval( const std::vector<double>& x)const;
  std::vector<double> d2Eval( const std::vector<double>& x)const;


private:
  double i_eval(std::size_t i, std::size_t j,double t, double s)const
  {
    return s*yM_[i][j]+t*yM_[i+1][j]+t*s*(aM_[i][j]*s+bM_[i][j]*t);
  }
  double i_deval(std::size_t i, std::size_t j,double t, double s)const
  {
    return (yM_[i+1][j]-yM_[i][j])/dx_[i]+(s-t)*(aM_[i][j]*s+bM_[i][j]*t)/dx_[i]+t*s*(bM_[i][j]-aM_[i][j])/dx_[i];

  }
  double i_d2eval(std::size_t i, std::size_t j,double t)const
  {
    return 2*(bM_[i][j]-2*aM_[i][j]+(aM_[i][j]-bM_[i][j])*3*t)/dx_[i]/dx_[i];

  }

  bool get_index(double x,std::size_t& i, double& t, double& s)const
  {
    auto it=x_map_.lower_bound(x);
    if (it!=x_map_.end())
      {
        i=it->second;
        if (i>0)
          {
            t=(x-x_[i-1])/dx_[i-1];
            s=1-t;
            return true;
          }
        else if (x==x_[i])
          {
            i=i+1;
            t=0;
            s=1;
            return true;
          }
        else return false;
      }
    else
      return false;
  }

  std::vector<std::vector<double>> getIndexMatix(const std::vector<double> x)
  {
    auto n=x.size();
    std::vector<std::vector<double >> o(n,std::vector<double>(n,0.0));

    for (std::size_t i=1; i<n-1; ++i)
      {
        o[i][i-1]=1.0/(x[i]-x[i-1]);
        o[i][i+1]=1.0/(x[i+1]-x[i]);
        o[i][i]=2.0*(o[i][i-1]+o[i][i+1]);
      }
    o[0][0]= 1.0/std::pow(x[1]-x[0],2);
    o[0][2]=-1.0/std::pow(x[2]-x[1],2);
    o[0][1]=o[0][0]+o[0][2];
    o[n-1][n-1]=-1.0/std::pow(x[n-1]-x[n-2],2);
    o[n-1][n-3]= 1.0/std::pow(x[n-2]-x[n-3],2);
    o[n-1][n-2]=o[n-1][n-1]+o[n-1][n-3];
    return o;
  }

  std::vector<std::vector<double>>
  getYvector(const std::vector<double>& x, const std::vector<std::vector<double>>& y)
  {
    auto n=x.size();
    auto m=y[0].size();
    if (y.size()!=x.size())
      {
        std::cerr<<" unequal vector sizes \n";
        return {};
      }
    std::vector<std::vector<double>> o(n,std::vector<double>(m));
    for(std::size_t j=0; j<m; ++j)
      {
        o[0][j]=2.0*((y[1][j]-y[0][j])/std::pow(x[1]-x[0],3)
                    -(y[2][j]-y[1][j])/std::pow(x[2]-x[1],3));
        o[n-1][j]=2.0*((y[n-2][j]-y[n-3][j])/std::pow(x[n-2]-x[n-3],3)
                      -(y[n-1][j]-y[n-2][j])/std::pow(x[n-1]-x[n-2],3));
      }
    for (std::size_t i=1;i<n-1; ++i)
      for(std::size_t j=0; j<m; ++j)
        {
          o[i][j]=3.0*((y[i][j]-y[i-1][j])/std::pow(x[i]-x[i-1],2)
                      +(y[i+1][j]-y[i][j])/std::pow(x[i+1]-x[i],2));
        }
    return o;
  }


  std::vector<std::vector<double>> getKs(const std::vector<double>& x, const std::vector<std::vector<double>>& y)
  {
    auto ys=getYvector(x,y);
    auto A=getIndexMatix(x);
    auto Ainv=inv(A);

    return mult(Ainv,ys);

  }



  std::map<double,std::size_t> x_map_;
  std::size_t n_knots_;
  std::size_t n_cols_;
  std::vector<double> x_;
  std::vector<double> dx_;
  std::vector<std::vector<double>> yM_;
  std::vector<std::vector<double>> kM_;
  std::vector<std::vector<double>> aM_;
  std::vector<std::vector<double>> bM_;
  double lower_default_;
  double upper_default_;
};


#endif // SPLINES

