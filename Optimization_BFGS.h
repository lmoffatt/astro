#ifndef OPTIMIZATION_BFGS_H
#define OPTIMIZATION_BFGS_H


#include "Matrix.h"

#include <chrono>


#include "Optimization.h"

#ifndef PI____
#define PI____
const double PI  =3.141592653589793238463;
#endif


struct opt_max_point
{
  M_Matrix<double> b;
  M_Matrix<double> w;
  double fout;
  bool isValid=false;
};


struct opt_max_step:public opt_max_point
{
  opt_max_step(opt_max_point p): opt_max_point(p),G(),Hinv(){}
  opt_max_step(){}
  M_Matrix<double> G;
  M_Matrix<double> Hinv;
};


struct opt_max_run
{
  std::chrono::steady_clock::time_point startTime;
  std::size_t maxIter;
  double maxDur_in_min;
  double paramChangeMin;
  double PostLogLikChangeMin;
  double GradientNormPostMin;
};



struct opt_max_iter
{
  std::size_t i;
  double timeOpt;
  double timeIter;
  opt_max_step sample;
  opt_max_point candidate;
  double ParamChange()const
  {
    return maxAbs(sample.b-candidate.b);
  }
  double foutChange() const
  {
    return sample.fout-candidate.fout;
  }
  double normGradien()const
  {
    return maxAbs(sample.G);
  }
};
inline
std::ostream& operator<<(std::ostream& os, const opt_max_point& iter)
{
  os<<"parameters\t"<<iter.b;
  os<<"\tfout\t"<<iter.fout;
  os<<"\tisValid\t"<<iter.isValid;
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, const opt_max_step& iter)
{
  const opt_max_point& i=iter;
  os<<i;
  os<<"\nG\n"<<iter.G;
  os<<"\nHinv\n"<<iter.Hinv;
  return os;

}


inline
std::ostream& operator<<(std::ostream& os, const opt_max_iter& iter)
{
  os<<"iters\t"<<iter.i;
  os<<"\ttimeOpt\t"<<iter.timeOpt;
  os<<"\tsample\t"<<iter.sample;
  os<<"\tcandidate\t"<<iter.candidate;

  return os;
}

struct BFGS_optimal
{
  template<class F, class D>
  static opt_max_iter opt(const F& f,
                          const D& x,
                          const M_Matrix<double>& binit,
                          double alfaInit=1e-5,
                          std::size_t maxIter=100,
                          double maxTime=600)
  {

    opt_max_run r;
    r.startTime= std::chrono::steady_clock::now();
    r.maxIter=maxIter;
    r.maxDur_in_min=maxTime;
    r.paramChangeMin=1e-6*alfaInit;
    r.GradientNormPostMin=1e-6*alfaInit;
    r.PostLogLikChangeMin=1e-6;


    opt_max_iter iter;
    iter.i=0;
    iter.sample=new_point(f,x,binit);
    if (!iter.sample.isValid) return iter;
    iter.sample=init_H(f,x,iter.sample);
    iter.candidate=next_candidate(f,x,iter.sample,alfaInit);
    iter.sample=compute_H(f,x,iter.sample,iter.candidate);
    iter.candidate=next_candidate(f,x,iter.sample,alfaInit);
    iterAct(r,iter);

    while (!meetCovergenceCriteria(r,iter))
      {
        iter.sample=compute_H(f,x,iter.sample,iter.candidate);
        iter.candidate=next_candidate(f,x,iter.sample,alfaInit);
        iterAct(r,iter);
      }
    return iter;
  }
private:

  static void iterAct(const opt_max_run& r,opt_max_iter& iter)
  {
    ++iter.i;
    auto tnow=std::chrono::steady_clock::now();
    auto d=tnow-r.startTime;
    double t0=iter.timeOpt;
    iter.timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
    iter.timeIter=60*(iter.timeOpt-t0);

  }

  static bool meetCovergenceCriteria(const opt_max_run& r,const opt_max_iter& iter)
  {
    if (iter.sample.G.empty()|| iter.sample.b.empty()) return true;
    else
      {
        bool surpassDuration_=iter.timeOpt+iter.timeIter/60.0>=r.maxDur_in_min;
        bool surpassIter_=iter.i>=r.maxIter;
        bool  smallParamChange_=iter.ParamChange()<r.paramChangeMin;
        bool  smallPostLikChange=iter.foutChange()<r.PostLogLikChangeMin;
        bool  smallGradient_=iter.normGradien()<r.GradientNormPostMin;



        return surpassIter_||
            surpassDuration_||
            smallParamChange_||
            smallPostLikChange||
            smallGradient_||
            !iter.sample.isValid||
            !iter.candidate.isValid;

      }
  }



  template<class F, class D>
  static opt_max_point new_point(const F& f,
                                 const D& x,
                                 const M_Matrix<double>& b)
  {
    opt_max_point out;
    out.b=b;
    out.fout=f(x,b);
    out.isValid=!std::isnan(out.fout);
    return out;
  }

  template<class F, class Data>
  static opt_max_step init_H(const F& f,
                             const Data& x,
                             opt_max_point p)
  {
    opt_max_step out(p);
    out.G=compute_G(f,x,out.b,out.fout);
    out.Hinv=eye<double>(p.b.size());
    return out;
  }




  template<class F, class Data>
  static opt_max_step compute_H(const F& f,
                                const Data& x,
                                const opt_max_step& p0,
                                opt_max_point p)
  {
    opt_max_step out(p);
    out.G=compute_G(f,x,out.b,out.fout);
    if (out.G.empty()) return {};
    else
      {
        out.Hinv=p0.Hinv;

        M_Matrix<double> delta_x_=-out.b+p0.b;
        M_Matrix<double> delta_G=-out.G+p0.G;

        double s2=(multTransp(delta_x_,delta_G))[0];
        if (s2>0)
          {

            double sigma=pow(s2,0.5);
            M_Matrix<double> s=delta_x_/sigma;
            M_Matrix<double> y=delta_G/sigma;
            M_Matrix<double> ds=s-(y*p0.Hinv);
            M_Matrix<double> Hx=TranspMult(ds,s);

            out.Hinv+=Hx+Transpose(Hx)-TranspMult(s,s)*(multTransp(ds,y))[0] ;

            /**     last formula is from
    http://www.math.washington.edu/~burke/crs/408f/notes/nlp/direction.pdf

    */

          }
        return out;
      }
  }




  template<class F, class Data>
  static opt_max_point next_candidate(const F& f,
                                      const Data& x,
                                      const opt_max_step& sample,
                                      double& alpha_,
                                      std::size_t maxfevalLoop_=10,
                                      double Wolf_Condition_c1_=0.1,
                                      double Wolf_Condition_c2_=0.5)
  {
    if (sample.G.empty()) return {};
    auto d_=  -sample.G*sample.Hinv; //that determine the direction of search.
    double df0_=multTransp(d_,sample.G)[0];
    double alfa=0;
    double beta=INFINITY;
    double alphamin=1e-9;
    std::size_t ifevalLoop=0;
    std::size_t neval_=0;
    // this loop look after a value of a that satisfy Wolfe conditions
    double df1;
    opt_max_point candidate1;
    while (true)
      {
        M_Matrix<double> x01=sample.b;
        M_Matrix<double> x1;
        while(true)
          {
            x1=x01+d_*alpha_;
            candidate1=new_point(f,x,x1);
            neval_++;
            ifevalLoop++;
            if ((candidate1.isValid||(std::isnan(alpha_)))||(alpha_<alphamin))
              break;
            else
              alpha_/=10;
          }
        //Armijo rule
        if (((candidate1.fout>sample.fout+df0_*Wolf_Condition_c1_*alpha_)||
             (!candidate1.isValid))
            &&(ifevalLoop<maxfevalLoop_))
          {
            beta=alpha_;
            alpha_=0.5*(alfa+beta);
          }
        else
          {
            df1=compute_dG(f,x,candidate1.b,candidate1.fout,d_);
            neval_++;
            ifevalLoop++;
            //curvature conditions
            if ((df1<Wolf_Condition_c2_*df0_)&&ifevalLoop<maxfevalLoop_)
              {
                alfa=alpha_;
                if (beta==INFINITY)
                  {
                    alpha_=2*alfa;
                  }
                else
                  {
                    alpha_=0.5*(alfa+beta);
                  };
              }
            else
              break;
          }
      }

    /** if it is possible, apply a cubic interpolation phase                                                                              */

    double b1=df0_+df1-3*(-candidate1.fout+sample.fout)/alpha_;

    double alfa00=alpha_;
    if ((b1*b1-df0_*df1)>0)
      {
        double b2=pow((b1*b1-df0_*df1),0.5);
        if (std::abs(df1-df0_+2*b2)>1e-6)
          alpha_=alpha_-alpha_*(df1+b2-b1)/(df1-df0_+2*b2);
      }
    if (!std::isfinite(alpha_)) alpha_=alfa00;
    M_Matrix<double> x2=sample.b+d_*alpha_;
    opt_max_point candidate2=new_point(f,x,x2);
    neval_++;
    if ((candidate2.isValid)&&(candidate2.fout>candidate1.fout))
      {
        return candidate2;
      }
    else
      {
        return candidate1;
      }

  }






  template<class F, class Data>
  static M_Matrix<double> compute_G(const F& f,
                                    const Data& x,
                                    const M_Matrix<double>& b,
                                    double y_0,
                                    double delta=1e-7,
                                    double delta_div=10,
                                    double deltamin=1e-10)
  {
    auto npar=b.size();
    M_Matrix<double> out(1,npar);
    for (std::size_t j=0; j<npar; ++j)
      {
        double deltarun=delta;

        M_Matrix<double> p(b);
        p[j]+=deltarun;
        double y_i=f(x,p);
        while (std::isnan(y_i)&&deltarun>deltamin)
          {
            deltarun=deltarun/delta_div;
            p=b;
            p[j]+=deltarun;
            y_i=f(x,p);
          }
        if (std::isnan(y_i))
          return {};

        out(0,j)=(y_i-y_0)/deltarun;

      }
    return out;



  }


  template<class F, class D>
  static double compute_dG(const F& f,
                           const D& x,
                           const M_Matrix<double>& b,
                           double y_0,
                           const M_Matrix<double>& direction,
                           double delta=1e-5,
                           double delta_div=10,
                           double deltamin=1e-7)
  {
    double out;
    double deltarun=delta;
    M_Matrix<double> p(b);
    p+=direction*deltarun;
    double y_i=f(x,p);
    while (std::isnan(y_i)&&deltarun>deltamin)
      {
        deltarun=deltarun/delta_div;
        p=b;
        p+=direction*deltarun;
        y_i=f(x,p);
      }
    if (std::isnan(y_i))
      return {};

    out=(y_i-y_0)/deltarun;


    return out;

  }


};














#endif // OPTIMIZATION_BFGS_H

