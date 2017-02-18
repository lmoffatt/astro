#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "Matrix.h"

#include <chrono>

#ifndef PI____
#define PI____
const double PI  =3.141592653589793238463;
#endif




template<typename T>
double Normal(const M_Matrix<T>& x,const M_Matrix<T>& m,const M_Matrix<T>& s_or_cov )
{
  if (s_or_cov.size()==m.size())
    {
      auto& s=s_or_cov;
      double logL=0;
      for (std::size_t i=0; i<m.size(); ++i)
        {
          double v=sqr(s[i]);
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






#endif // OPTIMIZATION_H
