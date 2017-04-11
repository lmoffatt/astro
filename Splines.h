#ifndef SPLINES
#define SPLINES

#include "MatrixInverse.h"


#include <vector>
#include <map>
#include <cmath>
#include <iostream>
# include <limits>


inline
double bisect(double (*f)(double),double targetf,double xpos,double xneg, std::size_t maxIter)
{
  double ypos, yneg;
  ypos=f(xpos)-targetf;
  if (ypos==0)
    return xpos;
  else if(ypos<0)
    {
      std::swap(xpos,xneg);
      yneg=ypos;
    }
  yneg=f(xneg)-targetf;
  if (yneg==0)
    return xneg;
  else if (yneg>0)
    return std::numeric_limits<double>::quiet_NaN();
  else
    {
      std::size_t niter=0;
      double c=(xpos+xneg)/2;
      while (niter<maxIter)
        {
          double yc=f(c)-targetf;
          ++niter;
          if (yc==0)
            return c;
          else if (yc>0)
            {
              xpos=c;
            }
          else
            {
              xneg=c;
            }
          c=(xpos+xneg)/2;

        }
      return c;
    }

}
inline
double secante(double (*f)(double),double targetf,double xpos,double xneg, std::size_t maxIter)
{
  double ypos, yneg;
  ypos=f(xpos)-targetf;
  if (ypos==0)
    return xpos;
  else if(ypos<0)
    {
      std::swap(xpos,xneg);
      yneg=ypos;
      ypos=f(xpos)-targetf;
    }
  else
    {
      yneg=f(xneg)-targetf;
    }
  if (yneg==0)
    return xneg;
  else if (yneg>0)
    return std::numeric_limits<double>::quiet_NaN();
  else
    {
      std::size_t niter=0;
      double c=xpos-ypos*(xpos-xneg)/(ypos-yneg);
      while (niter<maxIter)
        {
          double yc=f(c)-targetf;
          ++niter;
          if (yc==0)
            return c;
          else if (yc>0)
            {
              xpos=c;
              ypos=yc;
            }
          else
            {
              xneg=c;
              yneg=yc;
            }
          c=xpos-ypos*(xpos-xneg)/(ypos-yneg);

        }
      return c;
    }
}


inline
double exp1_over_x(double x)
{
  if (std::abs(x)>1e-3)
    return std::expm1(x)/x;
  else
    return 1+x/2.0+x*x/6.0+x*x*x/24.0+x*x*x*x/120.0;
}



class inverse
{
public:
  inverse(double (*func)(double)
          ,std::vector<double> vals
          ,std::size_t maxIter)
    : f_(func),
      finv_(),
      maxIter_(maxIter)
  {
    for (double x:vals)
      finv_[f_(x)]=x;
  }

  double eval(double x)
  {
    auto it=finv_.upper_bound(x);
    if (it==finv_.begin())
      {
        return lowerbound_;
      }
    else if (it==finv_.end())
      return upperbound_;
    else
      {
        double max=it->second;
        --it;
        double min=it->second;
        return secante(f_,x,min,max,maxIter_);
      }
  }

private:
  double (*f_)(double);
  std::map<double,double> finv_;
  std::size_t maxIter_;
  double lowerbound_;
  double upperbound_;
};




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
    return o;
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
  MSpline(const std::vector<double>& x, const std::vector<std::vector<double>>& y, bool NotAKnot):
    x_map_()
  ,n_knots_(x.size()-1)
  ,n_cols_(y[0].size())
  , x_(x)
  , dx_(std::vector<double>(n_knots_))
  , yM_(y)
  , kM_(getKs(x,y,NotAKnot))
  , aM_(std::vector<std::vector<double>>(n_knots_,std::vector<double>(n_cols_)))
  , bM_(std::vector<std::vector<double>>(n_knots_,std::vector<double>(n_cols_)))
  ,notAKnot_(NotAKnot)
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

  std::vector<std::vector<double>> getIndexMatix(const std::vector<double> x, bool notAKnot)
  {
    auto n=x.size();
    std::vector<std::vector<double >> o(n,std::vector<double>(n,0.0));

    for (std::size_t i=1; i<n-1; ++i)
      {
        o[i][i-1]=1.0/(x[i]-x[i-1]);
        o[i][i+1]=1.0/(x[i+1]-x[i]);
        o[i][i]=2.0*(o[i][i-1]+o[i][i+1]);
      }
    if (notAKnot)
      {
        o[0][0]= 1.0/std::pow(x[1]-x[0],2);
        o[0][2]=-1.0/std::pow(x[2]-x[1],2);
        o[0][1]=o[0][0]+o[0][2];
        o[n-1][n-1]=-1.0/std::pow(x[n-1]-x[n-2],2);
        o[n-1][n-3]= 1.0/std::pow(x[n-2]-x[n-3],2);
        o[n-1][n-2]=o[n-1][n-1]+o[n-1][n-3];
      }
    else
      {
        o[0][0]= 2.0/(x[1]-x[0]);
        o[0][1]=1.0/(x[1]-x[0]);
        o[n-1][n-1]=2.0/(x[n-1]-x[n-2]);
        o[n-1][n-2]=1.0/(x[n-1]-x[n-2]);

      }
    return o;
  }

  std::vector<std::vector<double>>
  getYvector(const std::vector<double>& x
             , const std::vector<std::vector<double>>& y
             ,bool notAKnot)
  {
    auto n=x.size();
    auto m=y[0].size();
    if (y.size()!=x.size())
      {
        std::cerr<<" unequal vector sizes \n";
        return {};
      }
    std::vector<std::vector<double>> o(n,std::vector<double>(m));
    if (notAKnot)
      for(std::size_t j=0; j<m; ++j)
        {
          o[0][j]=2.0*((y[1][j]-y[0][j])/std::pow(x[1]-x[0],3)
              -(y[2][j]-y[1][j])/std::pow(x[2]-x[1],3));
          o[n-1][j]=2.0*((y[n-2][j]-y[n-3][j])/std::pow(x[n-2]-x[n-3],3)
              -(y[n-1][j]-y[n-2][j])/std::pow(x[n-1]-x[n-2],3));
        }
    else
      for(std::size_t j=0; j<m; ++j)

        {
          o[0][j]=3.0*(y[1][j]-y[0][j])/std::pow(x[1]-x[0],2);
          o[n-1][j]=3.0*(y[n-1][j]-y[n-2][j])/std::pow(x[n-1]-x[n-2],2);

        }
    for (std::size_t i=1;i<n-1; ++i)
      for(std::size_t j=0; j<m; ++j)
        {
          o[i][j]=3.0*((y[i][j]-y[i-1][j])/std::pow(x[i]-x[i-1],2)
              +(y[i+1][j]-y[i][j])/std::pow(x[i+1]-x[i],2));
        }
    return o;
  }


  std::vector<std::vector<double>> getKs(const std::vector<double>& x, const std::vector<std::vector<double>>& y,bool notAKnot)
  {
    auto ys=getYvector(x,y,notAKnot);
    auto A=getIndexMatix(x,notAKnot);
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
  bool notAKnot_;
};



class MInterpol
{
public:
  MInterpol(const std::vector<double>& x, const std::vector<std::vector<double>>& y):
    x_map_()
  ,n_knots_(x.size()-1)
  ,n_cols_(y[0].size())
  , x_(x)
  , yM_(y)
  {

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

  std::vector<std::vector<double>> eval(const std::vector<double> x)
  {
    std::vector<std::vector<double>> o(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      o[i]=eval(x[i]);
    return o;
  }


private:
  double i_eval(std::size_t i, std::size_t j,double t, double s)const
  {
    return s*yM_[i][j]+t*yM_[i+1][j];
  }

  bool get_index(double x,std::size_t& i, double& t, double& s)const
  {
    auto it=x_map_.lower_bound(x);
    if (it!=x_map_.end())
      {
        i=it->second;
        if (i>0)
          {
            t=(x-x_[i-1])/(x_[i]-x_[i-1]);
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






  std::map<double,std::size_t> x_map_;
  std::size_t n_knots_;
  std::size_t n_cols_;
  std::vector<double> x_;
  std::vector<std::vector<double>> yM_;
  double lower_default_;
  double upper_default_;
};

class MQSpline
{
public:
  MQSpline(const std::vector<double>& x, const std::vector<std::vector<double>>& y):
    x_map_()
  ,n_knots_(x.size()-1)
  ,n_cols_(y[0].size())
  , x_(x)
  , yM_(y)
  , mM_(getMs(x,y))
  {
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


  std::vector<std::vector<double>> eval(const std::vector<double> x)
  {
    std::vector<std::vector<double>> o(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      o[i]=eval(x[i]);
    return o;
  }



private:
  double i_eval(std::size_t i, std::size_t j,double t, double s)const
  {
    return s*yM_[i][j]+t*yM_[i+1][j]+t*s*(mM_[i][j]);
  }

  bool get_index(double x,std::size_t& i, double& t, double& s)const
  {
    auto it=x_map_.lower_bound(x);
    if (it!=x_map_.end())
      {
        i=it->second;
        if (i>0)
          {
            t=(x-x_[i-1])/(x_[i]-x_[i-1]);
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


  std::vector<std::vector<double>> getMs(const std::vector<double>& x, const std::vector<std::vector<double>>& y)
  {

    std::vector<std::vector<double>> o(x.size(),std::vector<double> (y[0].size()));
    double v1=1.0/(x_[1]-x_[0]);
    double v2=1.0/(x_[2]-x_[1]);

    for (std::size_t j=0; j<n_cols_; ++j)
      o[0][j]=y[1][j]-v1/(v1+v2)*y[0][j]-v2/(v1+v2)*y[2][j];
    for (std::size_t i=1; i<n_knots_; i++)
      {
        v1=v2;
        v2=1.0/(x_[i+1]-x_[i]);
        for (std::size_t j=0; j<n_cols_; ++j)
          o[i][j]=((v1+v2)*y[i][j]-v1*(y[i-1][j]+o[i-1][j]))/v2-y[i+1][j];

      }
    return o;

  }



  std::map<double,std::size_t> x_map_;
  std::size_t n_knots_;
  std::size_t n_cols_;
  std::vector<double> x_;
  std::vector<std::vector<double>> yM_;
  std::vector<std::vector<double>> mM_;
  double lower_default_;
  double upper_default_;
};



/// los primeros dos intervalos DEBEN ser iguales
/// el vector sumy contiene la suma de los ys de cada intervalo
class MIExpSpline
{
public:
  MIExpSpline(const std::vector<double>& x, const std::vector<std::vector<double>>& sumy):
   inv_(getInv())
   , x_map_()
  ,n_knots_(x.size()-1)
  ,n_cols_(sumy[0].size())
  , x_(x)
  , YM_(sumy)
  , yM_(std::vector<std::vector<double>>(n_knots_,std::vector<double>(n_cols_)))
  , bM_(std::vector<std::vector<double>>(n_knots_,std::vector<double>(n_cols_,0.0)))
  {
    for (std::size_t i=0; i<x.size(); ++i)
      x_map_[x[i]]=i;
    for (std::size_t j=0; j<n_cols_; ++j)
      {

        double b=log((YM_[2][j]-YM_[1][j])/(YM_[1][j]-YM_[0][j]))/(x_[1]-x_[0]);
        bM_[0][j]=b;
        double y=-bM_[0][j]*(YM_[1][j]-YM_[0][j])/std::expm1(bM_[0][j]*(x_[0]-x_[1]));
        yM_[0][j]=y;
        bM_[1][j]=b;
        double y1=yM_[0][j]*exp(bM_[1][j]*(x_[2]-x_[1]));
        yM_[1][j]=y1;
      }
    for (std::size_t i=2; i<n_knots_; ++i)
      for (std::size_t j=0; j<n_cols_; ++j)
        {
          double b=inv_.eval((YM_[i+1][j]-YM_[i][j])/yM_[i-1][j]/(x_[i+1]-x_[i]))/(x_[i+1]-x_[i]);
          bM_[i][j]=b;
          double y=yM_[i-1][j]*std::exp(b*(x_[i+1]-x_[i]));
          yM_[i][j]=y;

        }
  }




  std::vector<double> eval( double x)const
  {
    std::vector<double> o(n_cols_);
    std::size_t i=0; double t;
    if (get_index(x,i,t))
      {
        for (std::size_t j=0; j<n_cols_; ++j)
          o[j]= i_eval(i,j,t);
        return o;
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



private:
  inverse inv_;
  static inverse getInv()
  {
    std::vector<double> v{-1e9,-1e7,-1e5,-1e4,-1e3,-1e2,-10 ,-7, -5, -4, -3 ,-2,-1.5,-1,0.7,-0.5,-0.4,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9,1,1.5,2,3,5,7,9,15,20,30,50,100};


    return inverse(exp1_over_x,v,8);
  }

  double i_eval(std::size_t i, std::size_t j,double t)const
  {
    //double y= yM_[i-1][j];
    //double b=bM_[i-1][j];
    //double e=std::expm1(bM_[i-1][j]*t);
    //double yb=y/b;
    return YM_[i][j]+yM_[i-1][j]/bM_[i-1][j]*std::expm1(bM_[i-1][j]*t);

  }

  bool get_index(double x,std::size_t& i, double& t)const
  {
    auto it=x_map_.upper_bound(x);
    if ((it!=x_map_.end()))
      {
        i=it->second;
        if (it!=x_map_.begin())
         {
           t=(x-x_[i]);
        return true;
          }
        else return false;
      }
    else return false;
  }



  std::map<double,std::size_t> x_map_;
  std::size_t n_knots_;
  std::size_t n_cols_;
  std::vector<double> x_;
  std::vector<std::vector<double>> YM_;
  std::vector<std::vector<double>> yM_;
  std::vector<std::vector<double>> bM_;
  double lower_default_;
  double upper_default_;
};




#endif // SPLINES

