#include "Splines.h"

Spline::Spline(std::vector<double> x, std::vector<double> y):
  x_map_()
,n_knots_(x.size()-1)

, x_(x)
, y_(y)
, k_(getKs(x,y))
, a_(std::vector<double>(n_knots_))
, b_(std::vector<double>(n_knots_))
, dx_(std::vector<double>(n_knots_))

{
  for (std::size_t i=0; i<n_knots_; ++i)
    {
      dx_[i]=x_[i+1]-x_[i];
      a_[i]=k_[i]*dx_[i]-y_[i+1]-y_[i];
      b_[i]=-k_[i+1]*dx_[i]+y_[i+1]-y_[i];
    }
  for (std::size_t i=0; i<x.size(); ++i)
    x_map_[x[i]]=i;
}




inverse MIExpSpline::inv_=getInv();

