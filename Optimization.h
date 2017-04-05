#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "Matrix.h"

#include <chrono>

#ifndef PI____
#define PI____
const double PI  =3.141592653589793238463;
#endif




template<typename T>
double Normal(const M_Matrix<T>& x,const M_Matrix<T>& m,const M_Matrix<T>& s_or_cov , bool is_Variance=false)
{
  if (s_or_cov.size()==m.size())
    {
      auto& s=s_or_cov;
      double logL=0;
      for (std::size_t i=0; i<m.size(); ++i)
        {
          double v;
          if (is_Variance)
            v=s[i];
          else
            v=sqr(s[i]);
          logL+=-sqr(x[i]-m[i])/v-std::log(2*PI*v);
        }
      logL/=2;
      return logL;
    }
  else
    {
      auto& cov=s_or_cov;

      auto covinv=inv(cov);
      auto cho_cov_=chol(cov,"upper");
      auto logDetCov_=log(diagProduct(cho_cov_));

      return -0.5*x.size()*log(PI)-logDetCov_-0.5*xTSigmaX(x.toVector_of_Cols()-m.toVector_of_Cols(),covinv);

    }
}

inline double Normal(double x,double m,double  s , bool is_Variance=false)
{
  double v;
  if (is_Variance)
    v=s;
  else
    v=sqr(s);

  double logL=-(sqr(x-m)/v+std::log(2*PI*v))/2;
  return logL;
}



template<typename T>
double Normal(const M_Matrix<T>& x,double m,double  s , bool is_Variance=false)
{
  double v;
  if (is_Variance)
    v=s;
  else
    v=sqr(s);

  double logL=0;
  for (std::size_t i=0; i<x.size(); ++i)
    {
      logL+=-sqr(x[i]-m)/v;
    }

  logL/=2;
  logL-=std::log(2*PI*v)*x.size()/2;

  return logL;
}

template<typename T>
double Normal(const std::vector<T>& x,double m,double  s , bool is_Variance=false)
{
  double v;
  if (is_Variance)
    v=s;
  else
    v=sqr(s);

  double logL=0;
  for (std::size_t i=0; i<x.size(); ++i)
    {
      logL+=-sqr(x[i]-m)/v;
    }

  logL/=2;
  logL-=std::log(2*PI*v)*x.size()/2;

  return logL;
}



template<typename T>
M_Matrix<T> NormalSample(const M_Matrix<T>& m,const M_Matrix<T>& s_or_cov , std::mt19937_64& mt)
{

  M_Matrix<T> x(m);
  if (s_or_cov.size()==m.size())
    {
      auto& s=s_or_cov;
      std::normal_distribution<> normal;
      for (std::size_t i=0; i<m.size(); ++i)
        {
          double z=normal(mt);
          x[i]+=s[i]*z;
        }
      return x;
    }
  else
    {
      auto& cov=s_or_cov;

      auto cho_cov_=chol(cov,"upper");
      auto z(m);
      std::normal_distribution<> normal;
      for (std::size_t i=0; i<m.size(); ++i)
        {
          z[i]=normal(mt);
        }

      M_Matrix<T> r=m+(z*cho_cov_);
      return r;

    }
}



template<typename T>
double Normal(const M_Matrix<T>& yfit,const M_Matrix<T>& y,std::size_t df=0)
{
  double SS=0;
  for (std::size_t i=0; i<y.size(); ++i)
    {
      SS+=sqr(yfit[i]-y[i]);
    }
  auto n=y.size();
  return -0.5*n*log(2*PI)-0.5*n*log(SS/(n-df))-0.5*(n-df);
}





template<typename T>
M_Matrix<T> Gradient_of_Normal(const M_Matrix<T>& beta,const M_Matrix<T>& m,const M_Matrix<T>& s_or_cov )
{
  if (s_or_cov.size()==m.size())
    {
      auto& s=s_or_cov;
      M_Matrix<T> G(1,m.size());
      for (std::size_t i=0; i<m.size(); ++i)
        {
          double v=sqr(s[i]);
          G[i]=-(beta[i]-m[i])/v;
        }
      return G;
    }
  else
    {
      auto& cov=s_or_cov;

      auto covinv=inv(cov);
      return -(beta.toVector_of_Cols()-m.toVector_of_Cols())*covinv;

    }
}


inline double  Gradient_of_Normal(double  beta,const double m,double s, bool is_variance=false)
{
  double v;
  if (!is_variance)
    double v=sqr(s);
  else
    v=s;
  double g=-(beta-m)/v;
  return g;

}



struct fit_point
{
  M_Matrix<double> b;
  M_Matrix<double> x;
  M_Matrix<double> y;
  M_Matrix<double> yfit;
  M_Matrix<double> w;

  double logPrior;
  M_Matrix<double> priorG;
  M_Matrix<double> priorH;

  double logL;
  bool isValid=false;
};


struct fit_step:public fit_point
{
  fit_step(fit_point p): fit_point(p),G(),H(){}
  fit_step(){}
  M_Matrix<double> G;
  M_Matrix<double> H;
};


struct fit_run
{
  std::chrono::steady_clock::time_point startTime;
  std::size_t maxIter;
  double maxDur_in_min;
  double maxLanda;
  double paramChangeMin;
  double PostLogLikChangeMin;
  double GradientNormPostMin;
};



struct fit_iter
{
  std::size_t i;
  double timeOpt;
  double timeIter;
  double landa;
  fit_step sample;
  fit_point candidate;
  double ParamChange()const
  {
    return maxAbs(sample.b-candidate.b);
  }
  double postLChange() const
  {
    return candidate.logL+candidate.logPrior-sample.logL-sample.logPrior;
  }
  double normGradien()const
  {
    return maxAbs(sample.G);
  }



};


struct golden_section_search
{

  template<class F>
  static double opt(const F& f,double a, double c,std::size_t itermax)
  {
    double goldenRatio = (1.0 + std::sqrt(5.0)) / 2.0;
    std::size_t iter=0;
    double fa=f(a); ++iter;
    double fc=f(c); ++iter;
    double b= c - (2.0 - goldenRatio) * (c - a);
    double  fb=f(b); ++iter;
    return step(f,a,fa,b,fb,c,fc,iter,itermax);
  }
  template<class F>
  static double step(const F& f,
                     double a, double fa,
                     double b, double fb,
                     double c, double fc,
                     std::size_t iter,
                     std::size_t itermax)
  {

    double goldenRatio = (1.0 + std::sqrt(5.0)) / 2.0;

    double x;
    if (b < c)
      x = b + (2.0 - goldenRatio) * (c - b);
    else
      x = b - (2.0 - goldenRatio) * (b - a);
    if (iter>itermax)
      return (a+c)/2;
    double  fx=f(x);
    if (fx < fb)
      {if (b < c)
          return step(f, b, fb,x,fx, c,fc, ++iter,itermax);
        else
          return step(f, a, fa,x, fx,b, fb,++iter,itermax);
      }
    else
      {
        if (b < c)

          return step(f, a, fa,b,fb, x,fx, ++iter,itermax);
        else
          return step(f, x, fx, b, fb,c, fc,++iter,itermax);

      }


  };

};


template<bool minimize>
struct wolf_conditions
{
  double c1=1E-4;
  double c2=0.9;

  struct termination
  {
    double alfa;
    M_Matrix<double> x;
    double f;
    M_Matrix<double> G;
    bool isValid;
    bool firstCondition;
    bool secondCondition;
    std::size_t nIter;

    bool isGood()const {return isValid&&firstCondition&&secondCondition;}
  };


  template<class F, class G>
  bool first(const F& f, const G& g, M_Matrix<double> x, M_Matrix<double> p, double alfa)const
  {
    if (minimize)
      return f(x+p*alfa)<=f(x)+c1*alfa*TranspMult(p,g(x));
    else
      return f(x+p*alfa)>=f(x)+c1*alfa*TranspMult(p,g(x));

  }

  bool first(double f_x,const M_Matrix<double>& g_x, const M_Matrix<double>& p,double alfa,double f_xa)const
  {
    if (!std::isfinite(f_xa))
      return false;
    else  if (minimize)
      return f_xa<=f_x+c1*alfa*TranspMult(p,g_x)[0];
    else return f_xa>=f_x+c1*alfa*TranspMult(p,g_x)[0];

  }

  template<class F, class G>
  bool second(const F& f, const G& g, M_Matrix<double> x, M_Matrix<double> p, double alfa)const
  {
    if (minimize)
      return TranspMult(p,g(x+p*alfa))[0]>=c2*TranspMult(p, g(x))[0];
    else
      return TranspMult(p,g(x+p*alfa))[0]<=c2*TranspMult(p, g(x))[0];

  }

  bool second(const M_Matrix<double>& g_x,const M_Matrix<double>& g_xa, const M_Matrix<double>& p)const
  {
    if (!all(g_xa,[](double x){return std::isfinite(x);}))
      return false;
    else if (minimize)
      return TranspMult(p,g_xa)[0]>=TranspMult(p, g_x)[0]*c2;
    else
      return TranspMult(p,g_xa)[0]<=TranspMult(p, g_x)[0]*c2;
  }

  template<class F, class G>
  bool third(const F& f, const G& g, M_Matrix<double> x, M_Matrix<double> p, double alfa)const
  {
    return std::abs(TranspMult(p,g(x+p*alfa))[0])<=c2*std::abs(TranspMult(p, g(x))[0]);
  }

  bool third(const M_Matrix<double>& g_x,const M_Matrix<double>& p,const M_Matrix<double>& g_xa)const
  {
    return std::abs(TranspMult(p,g_xa)[0])<=std::abs(TranspMult(p, g_x)[0])*c2;
  }





  template<class F, class G>
  termination
  opt(const F& f,
      const G& g,
      const M_Matrix<double>& x,
      double f_x,
      const M_Matrix<double>& g_x,
      M_Matrix<double> p,
      std::size_t maxEvals,
      double a_max=1.0,
      double a_min=0)const
  {
    std::size_t nEvals=0;
    double next_alfa=1;
    M_Matrix<double> next_xa=x+p*next_alfa;
    double next_f_xa=f(next_xa); ++nEvals;
    return step(f,g,x,f_x,g_x,p,next_alfa,next_xa,next_f_xa, a_min,a_max, nEvals,maxEvals);
  }

  template<class F, class G>
  termination
  second_step(const F& f,
              const G& g,
              const M_Matrix<double>& x,
              double f_x,
              const M_Matrix<double>& g_x,
              M_Matrix<double> p,
              double alfa,
              const M_Matrix<double>& xa,
              double f_xa,
              const M_Matrix<double>& g_xa,
              double a_max, std::size_t nEvals, std::size_t maxEvals
              )const

  {


    if ( second(g_x,g_xa,p))
      return termination{alfa,xa,f_xa,g_xa,true,true,true,nEvals};
    else if (nEvals>=maxEvals)
      {
        return termination{alfa,xa,f_xa,g_xa,true,true,false,nEvals};
      }
    else
      {
        double next_a_min=alfa;
        double next_alfa;
        if (std::isfinite(a_max))
          next_alfa=0.5*(next_a_min+a_max);
        else
          next_alfa=2*alfa;
        M_Matrix<double> next_xa=x+p*next_alfa;
        double next_f_xa=f(next_xa);++nEvals;
        return step(f,g,x,f_x,g_x,p,next_alfa,next_xa,next_f_xa,next_a_min,a_max,nEvals,maxEvals);

      }

  };

  template<class F, class G>
  termination
  step(const F& f,
       const G& g,
       const M_Matrix<double>& x,
       double f_x,
       const M_Matrix<double>& g_x,
       M_Matrix<double> p,
       double alfa,
       const M_Matrix<double>& xa,
       double f_xa,
       double a_min,
       double a_max
       , std::size_t nEvals, std::size_t maxEvals
       )const

  {


    if (first(f_x,g_x,p,alfa,f_xa))
      {
        M_Matrix<double> next_g_xa=g(xa);++nEvals;
        return second_step(f,g,x,f_x,g_x,p,alfa,xa,f_xa,next_g_xa, a_max,nEvals,maxEvals);
      }
    else if (nEvals>=maxEvals)
      {
        return {alfa,xa,f_xa,{},std::isfinite(f_xa),false,false,nEvals};
      }
    else
      {
        double next_a_max=alfa;
        double next_alfa=0.5*(a_min+next_a_max);
        M_Matrix<double> next_xa=x+p*next_alfa;
        double next_f_xa=f(next_xa);++nEvals;
        return step(f,g,x,f_x,g_x,p,next_alfa,next_xa,next_f_xa,                     a_min,next_a_max,nEvals,maxEvals);
      }


  };


};




struct LevenbergMarquardt_fit
{
  template<class F>
  static fit_iter opt(const F& f,
                      const M_Matrix<double>& x,
                      const M_Matrix<double>&y,
                      const M_Matrix<double>& bmean,
                      const M_Matrix<double>& bsd,
                      std::mt19937_64& mt,
                      double landa=1e3,
                      std::size_t maxIter=100,
                      double maxTime=60)
  {

    fit_run r;
    r.startTime= std::chrono::steady_clock::now();
    r.maxIter=maxIter;
    r.maxDur_in_min=maxTime;

    fit_iter iter;
    iter.landa=landa;
    iter.i=0;
    while(!iter.sample.isValid)
      {
        auto b0= NormalSample(bmean,bsd,mt);
        iter.sample=new_step(f,b0,x,y,bmean,bsd);
      }
    iter.candidate=next_candidate(f,bmean,bsd,iter.sample,iter.landa);

    while (!meetCovergenceCriteria(r,iter))
      {
        iter.sample=compute_H(f,iter.candidate);
        iter.candidate=next_candidate(f,bmean,bsd,iter.sample,iter.landa);
        iterAct(r,iter);
      }
    return iter;
  }
private:

  static void iterAct(const fit_run& r,fit_iter& iter)
  {
    ++iter.i;
    auto tnow=std::chrono::steady_clock::now();
    auto d=tnow-r.startTime;
    double t0=iter.timeOpt;
    iter.timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
    iter.timeIter=60*(iter.timeOpt-t0);

  }

  static bool meetCovergenceCriteria(const fit_run& r,const fit_iter& iter)
  {
    bool surpassDuration_=iter.timeOpt+iter.timeIter/60.0>=r.maxDur_in_min;
    bool surpassIter_=iter.i>=r.maxIter;
    bool  surpassLanda_=iter.landa>=r.maxLanda;
    bool  smallParamChange_=iter.ParamChange()<r.paramChangeMin;
    bool  smallPostLikChange=iter.postLChange()<r.PostLogLikChangeMin;
    bool  smallGradient_=iter.normGradien()<r.GradientNormPostMin;



    return surpassIter_||
        surpassDuration_||
        smallParamChange_||
        smallPostLikChange||
        smallGradient_||
        surpassLanda_||
        !iter.sample.isValid||
        !iter.candidate.isValid;

  }


  template<class F>
  static fit_point new_point(const F& f,
                             const M_Matrix<double>& bmean,
                             const M_Matrix<double>& bsd,
                             const fit_step& sample,
                             double landa)
  {
    M_Matrix<double> lHinv=inv(diag_landa(sample.H+sample.priorH,landa));
    M_Matrix<double> G=sample.G+sample.priorG;

    auto d=-G*lHinv;

    auto bnew=sample.b+d;
    return new_point(f,bmean,bsd,bnew,sample.x,sample.y);
  }

  template<class F>
  static fit_point new_point(const F& f,
                             const M_Matrix<double>& bmean,
                             const M_Matrix<double>& bsd,
                             const M_Matrix<double>& b,
                             const M_Matrix<double>& x,
                             const M_Matrix<double>&y)
  {
    fit_point out;
    out.b=b;
    out.x=x;
    out.y=y;
    out.logPrior=Normal(b,bmean,bsd);
    out.priorG=Gradient_of_Normal(b,bmean,bsd);
    if (bsd.size()==bmean.size())
      out.priorH=diag(bsd.apply([](double x){return std::pow(x,-2);}));
    else
      out.priorH=inv(bsd);
    out.yfit=f(b,x);
    out.logL=Normal(y,out.yfit);
    return out;
  }


  template<class F>
  static fit_step compute_H(const F& f,
                            fit_point p)
  {
    fit_step out(std::move(p));
    auto J=compute_J(f,out.b,out.x,out.yfit);
    auto epsilon=(out.y-out.yfit);
    out.G=J*epsilon.toVector_of_Rows();
    out.H=TranspMult(J,J);
    return out;
  }




  template<class F>
  static fit_point next_candidate(const F& f,
                                  const M_Matrix<double>& bmean,
                                  const M_Matrix<double>& bsd,
                                  const fit_step& sample,
                                  double& landa1,
                                  double v=3,
                                  double landaMax=1e9)
  {
    fit_point candidate,candidate0;

    double landa0;
    candidate=new_point(f,bmean,bsd,sample,landa1);
    std::size_t ifevalLoop=0;
    std::size_t maxFevalLoop=200;
    /// no es mejor
    if ((candidate.logL+candidate.logPrior<=sample.logL+sample.logPrior)
        ||!candidate.isValid)
      {
        /// mientras no sea mejor
        while(((candidate.logL+candidate.logPrior<=sample.logL+sample.logPrior)
               ||!candidate.isValid )&&(ifevalLoop<maxFevalLoop))

          {
            /// si me paso freno...
            if (landa1*v>=landaMax) break;
            landa0=landa1;
            landa1=landa0*v;
            candidate0=std::move(candidate);
            candidate=new_point(f,bmean,bsd,sample,landa1);
            ifevalLoop++;
            //   std::cerr<<landa_<<" ";
          }

        /// de aca sali porque cosegui algo mejor o porque el landa es muy grande
        /// o porque hice muchas evaluaciones
      }
    else
      {
        /// aqui ya tengo un candidato
        landa0=landa1;
        landa1=landa0/v;
        candidate0=candidate;

        candidate=new_point(f,bmean,bsd,sample,landa1);
        ifevalLoop++;
        /// mientras pueda mejorar
        while(candidate.isValid
              &&(candidate.logL+candidate.logPrior>candidate0.logL+candidate0.logPrior)
              &&(landa1>0.1))
          {
            landa0=landa1;
            landa1=landa0/v;

            candidate0=candidate;
            candidate=new_point(f,bmean,bsd,sample,landa1);
            ifevalLoop++;
          }

        /// si me pase, recupero el anterior
        if((candidate.logL+candidate.logPrior<=sample.logL+sample.logPrior)
           ||!candidate.isValid)
          {
            landa1=landa0;
            candidate=candidate0;
          }
      }
    return candidate;

  }






  template<class F>
  static M_Matrix<double> compute_J(const F& f,
                                    const M_Matrix<double>& b,
                                    const M_Matrix<double>& x,
                                    const M_Matrix<double>& y_0,
                                    double delta=1e-5,
                                    double delta_div=10,
                                    double deltamin=1e-7)
  {
    auto n=y_0.size();
    auto npar=b.size();
    M_Matrix<double> out(n,npar,0.0);
    for (std::size_t j=0; j<npar; ++j)
      {
        double deltarun=delta;

        M_Matrix<double> p(b);
        p[j]+=deltarun;
        M_Matrix<double> y_i=f(p,x);
        while ((isnan(y_i)||(y_i.empty()))&&deltarun>deltamin)
          {
            deltarun=deltarun/delta_div;
            p=b;
            p[j]+=deltarun;
            y_i=f(p,x);
          }
        if (isnan(y_i)||y_i.empty())
          return {};

        for (std::size_t i=0; i<n; ++i)
          out(i,j)=(y_i[i]-y_0[i])/deltarun;

      }
    return out;



  }



};




template<bool minimize>
struct Newton_fit_fossil
{
  struct fit_point
  {
    fit_point(const M_Matrix<double>& xs,double logLs ): x(xs),logL(logLs),isValid(std::isfinite(logLs)){}
    fit_point(){}

    M_Matrix<double> x;
    double logL;
    bool isValid=false;
  };


  struct fit_step:public fit_point
  {
    fit_step(){}
    fit_step(fit_point p): fit_point(p),G(),H(){}
    M_Matrix<double> G;
    M_Matrix<double> H;
    double landa;
    M_Matrix<double> Hlinv;
    bool validCov=false;
    M_Matrix<double> d;
  };


  struct fit_iter
  {
    std::size_t i;
    double timeOpt;
    double timeIter;
    fit_step sample;
    double landa;
    double alfa;
    fit_point candidate;
    double ParamChange()const
    {
      return maxAbs(sample.x-candidate.x);
    }
    double postLChange() const
    {
      return candidate.logL-sample.logL;
    }
    double normGradien()const
    {
      return maxAbs(sample.G);
    }


  };

  struct fit_run
  {
    std::chrono::steady_clock::time_point startTime;
    std::size_t maxIter;
    std::size_t maxEvalLoop;
    double maxDur_in_min;
    double maxLanda=1e10;
    double minLanda=1e-3;
    double paramChangeMin=1e-7;
    double PostLogLikChangeMin=1e-7;
    double GradientNormPostMin=1e-5;

    double vlanda=3;
    double min_alfa=0.3;

  };

  template<class L, class G,class H, class Inv>
  static fit_iter opt(const L& logL,
                      const G& g,
                      const H& h,
                      const Inv& inverse,
                      const M_Matrix<double>& x,
                      double landa=1e-3,
                      double vlanda=3,
                      std::size_t maxIter=1000,
                      double maxTime=60)
  {
    fit_run r;
    r.startTime= std::chrono::steady_clock::now();
    r.maxIter=maxIter;
    r.maxDur_in_min=maxTime;
    r.maxEvalLoop=20;
    r.vlanda=vlanda;
    r.min_alfa=1.0/vlanda;

    fit_iter iter;
    iter.i=0;
    double alfa=1;
    iter.landa=landa;

    iter.sample=new_step(logL,g,h,inverse,x,iter.landa,vlanda);
    iter.candidate=next_candidate(logL,g,iter.sample,alfa,r.maxEvalLoop);
    iter.alfa=alfa;
    if ((iter.alfa==1.0)&&(iter.sample.landa==iter.landa))
      iter.landa=std::max(r.minLanda,iter.landa/std::pow(r.vlanda,3));
    else if ((iter.alfa<r.min_alfa)&&(iter.landa<r.maxLanda))
      iter.landa=std::max(iter.landa,1.0)*r.vlanda;
    alfa=1;

    while (!meetCovergenceCriteria(r,iter))
      {
        iter.sample=compute_H(logL,g,h,inverse,iter.candidate,iter.landa,r.vlanda);
        iter.candidate=next_candidate(logL,g,iter.sample,alfa,r.maxEvalLoop);
        iter.alfa=alfa;
        if ((iter.alfa==1.0)&&(iter.sample.landa==iter.landa))
          iter.landa=std::max(r.minLanda,iter.landa/std::pow(r.vlanda,3));
        else if ((iter.alfa<r.min_alfa)&&(iter.sample.landa<r.maxLanda))
          iter.landa=std::max(iter.sample.landa,1.0)*r.vlanda;

        alfa=1;
        iterAct(r,iter);
      }
    return iter;
  }
private:
  template<class L, class G, class H, class Hinv>
  static fit_step new_step(const L& l, const G& g, const H& h, const Hinv& hinv, const M_Matrix<double>& x, double landa, double v)
  {
    double logL=l(x);
    fit_point p(x,logL);

    return compute_H(l,g,h,hinv,p,landa,v);

  }


  template<class L, class G, class H, class Hinv>
  static fit_step compute_H(const L& l, const G& g, const H& h, const Hinv& hinv, fit_point p, double landa, double v)
  {
    fit_step out(std::move(p));
    if (out.isValid)
      {
        out.G=g(out.x);
        out.H=h(out.x);
        out.landa=landa;
        auto Hl=out.H+diag((diag(out.H)).apply([](double x){return -std::abs(x);}))*out.landa;
        out.Hlinv=hinv(Hl);
        out.validCov=all(diag(out.Hlinv), [](double x){return x<0;});

        while (!out.validCov)
          {
            out.landa*=v;
            Hl=out.H+diag((diag(out.H)).apply([](double x){return -std::abs(x);}))*out.landa;
            out.Hlinv=hinv(Hl);
            out.validCov=all(diag(out.Hlinv), [](double x){return x<0;});
          }

        out.d=-out.G*out.Hlinv;
      }

    return out;

  }


  static void iterAct(const fit_run& r,fit_iter& iter)
  {
    ++iter.i;
    auto tnow=std::chrono::steady_clock::now();
    auto d=tnow-r.startTime;
    double t0=iter.timeOpt;
    iter.timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
    iter.timeIter=60*(iter.timeOpt-t0);

  }

  static bool meetCovergenceCriteria(const fit_run& r,const fit_iter& iter)
  {
    bool surpassDuration_=iter.timeOpt+iter.timeIter/60.0>=r.maxDur_in_min;
    bool surpassIter_=iter.i>=r.maxIter;
    bool  smallGradient_=iter.normGradien()<r.GradientNormPostMin;
    bool surpass_MaxLanda=iter.sample.landa>r.maxLanda;


    return surpassIter_||
        surpassDuration_||
        surpass_MaxLanda||
        smallGradient_||
        !iter.sample.isValid||
        !iter.candidate.isValid;

  }


  template<class L>
  static fit_point candidate_point(const L& logL,
                                   const fit_step& sample,
                                   double alfa)
  {

    auto dl=sample.d*alfa;

    auto bnew=sample.x+dl;
    return new_point(logL,bnew);
  }

  template<class L>
  static fit_point new_point(const L& logL,
                             const M_Matrix<double>& b)
  {
    fit_point out;
    out.x=b;
    out.logL=logL(b);
    if (std::isfinite(out.logL))
      out.isValid=true;
    else
      out.isValid=false;
    return out;
  }






  template<class L, class G>
  static fit_point next_candidate(const L& logL,
                                  const G& g,
                                  const fit_step& sample,
                                  double& alfa,
                                  std::size_t maxEvals)
  {
    wolf_conditions<minimize> w;
    M_Matrix<double> g_x=g(sample.x);

    auto run=w.opt(logL,g,sample.x,sample.logL,g_x,sample.d,maxEvals);
    fit_point candidate(run.x,run.f);
    alfa=run.alfa;
    return candidate;

  }






  template<class F>
  static M_Matrix<double> compute_J(const F& f,
                                    const M_Matrix<double>& b,
                                    const M_Matrix<double>& x,
                                    const M_Matrix<double>& y_0,
                                    double delta=1e-5,
                                    double delta_div=10,
                                    double deltamin=1e-7)
  {
    auto n=y_0.size();
    auto npar=b.size();
    M_Matrix<double> out(n,npar,0.0);
    for (std::size_t j=0; j<npar; ++j)
      {
        double deltarun=delta;

        M_Matrix<double> p(b);
        p[j]+=deltarun;
        M_Matrix<double> y_i=f(p,x);
        while ((isnan(y_i)||(y_i.empty()))&&deltarun>deltamin)
          {
            deltarun=deltarun/delta_div;
            p=b;
            p[j]+=deltarun;
            y_i=f(p,x);
          }
        if (isnan(y_i)||y_i.empty())
          return {};

        for (std::size_t i=0; i<n; ++i)
          out(i,j)=(y_i[i]-y_0[i])/deltarun;

      }
    return out;



  }



};



template<bool minimize, bool verbose>
struct Newton_fit
{
  struct fit_point
  {
    fit_point(const M_Matrix<double>& xs,double logLs ): x(xs),logL(logLs),isValid(std::isfinite(logLs)){}
    fit_point(){}

    M_Matrix<double> x;
    double logL;
    bool isValid=false;
  };


  struct fit_step:public fit_point
  {
    fit_step(){}
    fit_step(fit_point p): fit_point(p),G(),H(){}
    M_Matrix<double> G;
    M_Matrix<double> H;
  };


  class fit_iter
  {
  public:


    std::size_t i;
    double timeOpt;
    double timeIter;
    fit_step sample;
    double landa;
    M_Matrix<double> Hlinv;
    bool validCov=false;
    M_Matrix<double> d;
    double alfa;
    fit_point candidate;
    double ParamChange()const
    {
      return maxAbs(sample.x-candidate.x);
    }
    double postLChange() const
    {
      return candidate.logL-sample.logL;
    }
    double normGradien()const
    {
      return maxAbs(sample.G);
    }


  };

  struct fit_run
  {
    std::chrono::steady_clock::time_point startTime;
    std::size_t maxIter;
    std::size_t maxEvalLoop;
    double maxDur_in_min;
    double maxLanda=1e10;
    double minLanda=1e-3;
    double paramChangeMin=1e-7;
    double PostLogLikChangeMin=1e-7;
    double GradientNormPostMin=1e-5;

    double vlanda=3;
    double min_alfa=0.3;

  };

  struct landa_H
  {
    double landa;
    M_Matrix<double> Hlinv;
    bool validCov;
  };

  static bool isValidCov(const M_Matrix<double>& HlandaInv)
  {
    if (HlandaInv.empty()) return false;
    else
      {

        for (std::size_t i=0; i<HlandaInv.nrows(); ++i)
          {
            if (minimize)
              {
                if(HlandaInv(i,i)<=0) return false;
              }
            else
              {
                if(HlandaInv(i,i)>=0) return false;

              }
          }
        return true;

      }
  }

  template<class Hinv>
  static landa_H compute_landaH(const Hinv& hinv,const M_Matrix<double> H, double landa0, double v,double landamax)
  {

    auto Hl=H;
    if (minimize)
      {
        for (std::size_t i=0; i<H.nrows(); ++i)
          Hl(i,i)+=std::abs(H(i,i))*landa0;
      }
    else
      {
        for (std::size_t i=0; i<H.nrows(); ++i)
          Hl(i,i)-=std::abs(H(i,i))*landa0;
      }
    auto Hlinv=hinv(Hl);

    if (isValidCov(Hlinv))
      return {landa0,Hlinv,true};
    else if (landa0*v>=landamax)
      return {landa0,Hlinv,false};
    else
      return compute_landaH(hinv,H,landa0*v,v,landamax);
  }


  template<class L, class G, class H, class Hinv>
  static fit_iter step(const L& logL, const G& g, const H& h, const Hinv& hinv, fit_run r,
                       fit_iter iter)
  {
    if (meetCovergenceCriteria(r,iter))
      return iter;
    else
      {
        auto i=iter.i+1;

        double newLanda;
        fit_step sample(iter.candidate);
        sample.G=g(sample.x);
        sample.H=h(sample.x);

        if (iter.alfa<r.min_alfa)
          newLanda=iter.landa*r.vlanda;
        else if (iter.alfa>0.9)
          newLanda=iter.landa/r.vlanda;
        else
          newLanda=iter.landa;


        return candidate_point(logL,g,h,hinv,r,sample,newLanda,i);
      }



  }

  template<class L, class G, class H, class Hinv>
  static fit_iter candidate_point(const L& logL,
                                  const G& g,
                                  const H& h,
                                  const Hinv& hinv,
                                  fit_run r,
                                  fit_step sample,
                                  double landa,
                                  std::size_t i)
  {
    landa_H lH=compute_landaH(hinv,sample.H,landa,r.vlanda,r.maxLanda);
    auto newd=-sample.G*lH.Hlinv;
    wolf_conditions<minimize> w;
    typename wolf_conditions<minimize>::termination run=w.opt(logL,g,sample.x,sample.logL,sample.G,newd,r.maxEvalLoop);

    if ((!run.isGood())&&lH.landa<r.maxLanda)
      return candidate_point(logL,g,h,hinv,r,sample,lH.landa*r.vlanda,i);
    else
      {

        fit_iter iter;
        iter.i=i;
        auto tnow=std::chrono::steady_clock::now();
        auto dur=tnow-r.startTime;
        double t0=iter.timeOpt;
        iter.timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(dur).count()/60.0;
        iter.timeIter=60*(iter.timeOpt-t0);
        iter.sample=sample;
        iter.landa=lH.landa;
        iter.Hlinv=lH.Hlinv;
        iter.d=newd;
        iter.alfa=run.alfa;
        fit_point candidate(run.x,run.f);
        iter.candidate=std::move(candidate);
        if (verbose) std::cerr<<iter;
        return step(logL,g,h,hinv,r,iter);
      }
  }




  template<class L, class G,class H, class Inv>
  static fit_iter opt(const L& logL,
                      const G& g,
                      const H& h,
                      const Inv& inverse,
                      const M_Matrix<double>& x,
                      double landa=1e-3,
                      double vlanda=3,
                      std::size_t maxIter=1000,
                      double maxTime=60)
  {
    fit_run r;
    r.startTime= std::chrono::steady_clock::now();
    r.maxIter=maxIter;
    r.maxDur_in_min=maxTime;
    r.maxEvalLoop=20;
    r.vlanda=vlanda;
    r.min_alfa=1.0/vlanda;

    fit_step sample(fit_point(x,logL(x)));
    sample.G=g(x);
    sample.H=h(x);

    return  candidate_point(logL,g,h,inverse,r,sample,landa,0);

  }

  friend
  std::ostream& operator<<(std::ostream &os,const fit_iter& iter)
  {
    os<<iter.i<<" "<<iter.timeOpt<<" "<<iter.timeIter<<" "<<iter.sample.logL<<" "<<iter.landa<<"\n";
    return os;
  }

private:




  static bool meetCovergenceCriteria(const fit_run& r,const fit_iter& iter)
  {
    bool surpassDuration_=iter.timeOpt+iter.timeIter/60.0>=r.maxDur_in_min;
    bool surpassIter_=iter.i>=r.maxIter;
    bool  smallGradient_=iter.normGradien()<r.GradientNormPostMin;
    bool surpass_MaxLanda=iter.landa>r.maxLanda;


    return surpassIter_||
        surpassDuration_||
        surpass_MaxLanda||
        smallGradient_||
        !iter.sample.isValid||
        !iter.candidate.isValid;

  }

};


template <bool verbose>
inline std::ostream& operator<<(std::ostream& os, const typename Newton_fit<false,verbose>::fit_iter& it)
{
  os<<"iteration="<<it.i<<" timeOpt="<<it.timeOpt<<"timeIter="<<it.timeIter<<"\n";
  //os<<"landa="<<it.landa<<"\n";
  os<<"sample logL="<<it.sample.logL<<" candidate diff="<<it.candidate.logL-it.sample.logL<<"\n";
  os<<"Parameter change="<<it.ParamChange()<<" postLChange="<<it.postLChange();
  os<<" normGradient="<<it.normGradien()<<"\n";


  return os;
}



#endif // OPTIMIZATION_H
