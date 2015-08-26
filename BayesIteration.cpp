#include <vector>
#include <fstream>
#include <cmath>
#include <map>
#include "BayesIteration.h"
#include "MatrixInverse.h"

//std::vector<double> ABC_Freq_obs::getDataWeigth()
//{
////  std::vector<double> w=getDataStandardError();
////  for (std::size_t i=0; i<w.size(); i++)
////    {
////      w[i]=1.0/w[i]/w[i];
////    }
////  std::vector<double> wW(w);
////  return wW;
//}


ABC_Multinomial_Model::~ABC_Multinomial_Model(){}

BayesIteration::~BayesIteration(){
}



Parameters BayesIteration::Posterior()const
{
  return posterior_;
}

Parameters BayesIteration::Prior(std::size_t n)const
{
  return priors_[n];

}


BayesIteration::BayesIteration(const ABC_Multinomial_Model* m,
                               Parameters prior,
                               const ABC_Freq_obs* d,
                               const std::string& filename):
  m_(m),
  priors_(1,prior),
  data_(1,d),
  posterior_(),
  numSeeds_(20),
  filename_(filename)
{
  //getPosterior();
}


BayesIteration& BayesIteration::addNewData(ABC_Freq_obs* d)
{
  data_.push_back(d);
  priors_.push_back(posterior_);
  getPosterior();
}



BayesIteration::BayesIteration(const BayesIteration& other):
  m_(other.m_),
  data_(other.data_),
  priors_(other.priors_),
  posterior_(other.posterior_)
{}


/*   void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);

    BayesIteration();

    ~BayesIteration(){}
*/
// void reset(const SimParameters& sp,const Treatment& tr);






BayesIteration::BayesIteration(){}




std::vector<std::vector<double>> BayesIteration::p_exp(const Parameters& parameters)const
{

  std::vector<std::vector<double>> simulatedResults=m_->p_exp(parameters);
  std::vector<double> param=parameters.trMeans();

  //simulatedResults.insert(simulatedResults.end(),param.begin(),param.end());
  return simulatedResults;
}


const ABC_Freq_obs& BayesIteration::getData()const
{
//  std::vector<double>data =data_.back()->getData();
//  std::vector<double> param=priors_.back().pMeans();

//  data.insert(data.end(),param.begin(),param.end());
//  std::vector<double> d(data);
//  return d;


}

//std::vector<double> BayesIteration::getDataStandardError()const
//{

//  std::vector<double>se =data_.back()->getDataStandardError();
//  std::vector<double> param=priors_.back().pStds();

//  se.insert(se.end(),param.begin(),param.end());
//  std::vector<double> o(se);
//  return o;

//}

BayesIteration& BayesIteration::getPosterior(const Parameters& startingPoint)
{

  Parameters p=priors_.back();

  std::size_t numIterations=1000;

  LevenbergMarquardtMultinomial LM(this,
                        startingPoint,
                        numIterations);
  LM.optimize();
  std::ofstream f;

  f.open(filename_.c_str(),std::ios_base::app);

  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------The prior is-----------------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";

  put(f,p);



  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Start from the following point------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  put(f,startingPoint);


  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Result of Levenberg Marquardt------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<LM;
  f<<"SS \t"<<LM.LogLik()<<"\n";
  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
  f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
  f<<"logLikelihood \t"<<LM.LogLik()<<"\n";





  f<<"chi2Distance to seed\t"<<startingPoint.chi2Distance(LM.OptimParameters())<<"\n";
  f<<"dBDistance to seed\t"<<dbDistance(startingPoint,LM.OptimParameters())<<"\n";
  f<<"chi2Distance to prior\t"<<p.chi2Distance(LM.OptimParameters())<<"\n";
  f<<"dBDistance to prior\t"<<dbDistance(p,LM.OptimParameters())<<std::endl;
  put(f,LM.OptimParameters());
  f<<"SS \t"<<LM.LogLik()<<"\n";


  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Result of Hessian calculation---\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  Parameters h=getHessian(LM.OptimParameters(),1e-5);

  double logDetPost=log(det(h.getCovariance()));
  f<<"SS \t"<<LM.LogLik()<<"\n";

  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Evidence Hess\t"<<LM.getEvidence()-0.5*LM.logDetPostCov()+0.5*logDetPost<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrCov Hess\t"<<logDetPost<<"\n";

  h.write(f);


  Parameters corrMax=getEvidence(LM.OptimParameters(),1000);
  corrMax.write(f);

  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Result of Posterior LikelihoodSampling---\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  corrMax.write(f);
  logDetPost=log(det(corrMax.getCovariance()));
  double SSc=SumWeighedSquare(corrMax);
  f<<"SS \t"<<LM.LogLik()<<"\n";
  f<<"SS corr\t"<<SSc<<"\n";

  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Evidence corr SS\t"<<-0.5*LM.LogLik()+0.5*logDetPost<<"\n";
  f<<"Evidence corr SSc\t"<<-0.5*SSc+0.5*logDetPost<<"\n";

  f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
  f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrCovCorr \t"<<logDetPost<<"\n";
  f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";


  put(f,corrMax);



  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"chi2Distance to seed\t"<<startingPoint.chi2Distance(LM.OptimParameters())<<"\n";
  f<<"dBDistance to seed\t"<<dbDistance(startingPoint,LM.OptimParameters())<<"\n";
  f<<"chi2Distance to prior\t"<<p.chi2Distance(LM.OptimParameters())<<"\n";
  f<<"dBDistance to prior\t"<<dbDistance(p,LM.OptimParameters())<<std::endl;


  f.close();
  return *this;
}

BayesIteration& BayesIteration::getPosterior(const Parameters& startingPoint,
                                             double factor,
                                             std::size_t numSeeds,
                                             double probParChange)
{
//  std::vector<double> data=getData();
//  std::vector<double> w=getDataWeigth();

  Parameters p=priors_.back();

  std::size_t numIterations=1000;

  std::ofstream f;

  f.open(filename_.c_str(),std::ios_base::app);

  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------The prior is-----------------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";

  put(f,p);

  for (std::size_t i=0;i<numSeeds;i++)
    {

      Parameters seed;
      seed=startingPoint.randomSample(p,factor,probParChange);
      LevenbergMarquardtMultinomial LM(this,
                            seed,
                            numIterations);
      LM.optimize();






      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Start from the following point------------------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      put(f,seed);





      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Result of Levenberg Marquardt------------------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      if (LM.LogLik()<200)
        f<<LM;
      f<<"Suma de cuadrados"<<LM.LogLik();
      f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
      f<<"Evidence \t"<<LM.getEvidence()<<"\n";
      f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
      f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
      f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
      f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
      f<<"logLikelihood \t"<<LM.LogLik()<<"\n";

      f<<"chi2Distance to seed\t"<<startingPoint.chi2Distance(LM.OptimParameters())<<"\n";
      f<<"dBDistance to seed\t"<<dbDistance(startingPoint,LM.OptimParameters())<<"\n";
      f<<"chi2Distance to prior\t"<<p.chi2Distance(LM.OptimParameters())<<"\n";
      f<<"dBDistance to prior\t"<<dbDistance(p,LM.OptimParameters())<<"\n";

      put(f,LM.OptimParameters());
      f<<"SS \t"<<LM.LogLik()<<"\n";


    }

  f.close();
  return *this;

}



BayesIteration& BayesIteration::getPosterior()
{
  std::vector<LevenbergMarquardtMultinomial> LMs;
  std::vector<Parameters> Ps;


  Parameters p=priors_.back();


  std::size_t factor=2;
  std::size_t numIterations=30;
  std::size_t numSeeds=10;
  std::map<double,Parameters> seeds=getRandomParameters(numSeeds*factor,1);

  std::map<double,Parameters> friuts;

  LevenbergMarquardtMultinomial LM(this,
                        p,
                        numIterations);
  LM.optimize();
  LMs.push_back(LM);
  Ps.push_back(LM.OptimParameters());

  friuts[LM.LogLik()]=LM.OptimParameters();
  std::ofstream f;

  f.open(filename_.c_str(),std::ios_base::app);

  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------The prior is-----------------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";

  put(f,p);



  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Start from the center-------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<LM;
  put(f,Ps.back());
  f.close();


  std::map<double,Parameters>::iterator it=seeds.begin();

  for (std::size_t i=0; i<numSeeds; i++)
    {
      Parameters initParam=(*it).second;
      double ss=(*it).first;
      ++it;
      LevenbergMarquardtMultinomial LM(this,
                            initParam,
                            numIterations);
      LM.optimize();
      friuts[LM.LogLik()]=LM.OptimParameters();

      LMs.push_back(LM);
      Ps.push_back(LM.OptimParameters());
      f.open(filename_.c_str(),std::ios_base::app);
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Start from the perisphery-------------------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<LM;
      put(f,Ps.back());
      f.close();

    }
  double errorFactor=0.2;
  double errorShrinkFactor=3;

  std::size_t numCycle=3;
  for (size_t iCycle=0;
       iCycle<numCycle;
       iCycle++)
    {
      seeds=friuts;
      it=seeds.begin();
      size_t j=0;
      size_t jfactor=5;
      Parameters initParam=(*it).second;
      errorFactor/=errorShrinkFactor;



      for (std::size_t i=0; i<numSeeds; i++)
        {
          ++j;
          if(j==jfactor)
            {
              ++it;
              initParam=(*it).second;
              j=0;
            }
          initParam=initParam.randomSample(errorFactor);

          LevenbergMarquardtMultinomial LM(this,
                                initParam,
                                numIterations);
          LM.optimize();
          friuts[LM.LogLik()]=LM.OptimParameters();

          LMs.push_back(LM);
          Ps.push_back(LM.OptimParameters());
          f.open(filename_.c_str(),std::ios_base::app);
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";
          f<<"----------Start from the perisphery-------------------\n";
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";
          f<<LM;
          put(f,Ps.back());
          f.close();

        }


    }

  return *this;
}


std::ostream& BayesIteration::put(std::ostream& s,const Parameters& parameters)const
{
 // m_->put(s,parameters);
  return s;
}


void BayesIteration::setFilename(const std::string filename)
{
  filename_=filename;
}

double BayesIteration::SumWeighedSquare(const Parameters& p)
{
//  std::vector<double> d=this->getData();
//  std::vector<double> w=this->getDataWeigth();

//  std::vector<double> y=this->p_exp(p);
//  double SSW=0;
//  for (std::size_t i=0;i<d.size();i++)
//    {
//      SSW+=pow(y[i]-d[i],2)*w[i];
//    };
//  return SSW;
}

Parameters BayesIteration::getHessian(const Parameters& MAP,double eps)
{
  std::vector<std::vector<double> > Hessian(MAP.size(),std::vector<double> (MAP.size(),0));
  double ss=0.5*SumWeighedSquare(MAP);

  for (std::size_t i=0; i<MAP.size(); i++)
    {
      Parameters fip(MAP);
      fip[i]=fip[i]+eps;
      Parameters fin(MAP);
      fin[i]=fin[i]-eps;
      double ssip=0.5*SumWeighedSquare(fip);
      double ssin=0.5*SumWeighedSquare(fin);

      Hessian[i][i]=(ssip+ssin-2.0*ss)/eps/eps;
      for (std::size_t j=i+1; j<MAP.size(); j++)
        {
          Parameters fipjp(MAP);
          fipjp[i]=fipjp[i]+eps;
          Parameters fipjn(fipjp);
          fipjp[j]=fipjp[j]+eps;
          fipjn[j]=fipjn[j]-eps;

          Parameters finjp(MAP);
          finjp[i]=finjp[i]-eps;
          Parameters finjn(finjp);
          finjp[j]=finjp[j]+eps;
          finjn[j]=finjn[j]-eps;

          double ssipjp=0.5*SumWeighedSquare(fipjp);
          double ssinjp=0.5*SumWeighedSquare(finjp);
          double ssipjn=0.5*SumWeighedSquare(fipjn);
          double ssinjn=0.5*SumWeighedSquare(finjn);

          Hessian[i][j]=(ssipjp+ssinjn-ssipjn-ssinjp)/eps/eps/4.0;
          Hessian[j][i]=Hessian[i][j];
        }

    }
  Parameters result(MAP);
  result.setCovariance(inv(Hessian));
  return result;
}


Parameters BayesIteration::getHessianInterpol(const Parameters& MAP, double mindSS, double maxdSS)
{
  std::vector<std::vector<double> > Hessian(MAP.size(),std::vector<double> (MAP.size(),0));
  double ss=0.5*SumWeighedSquare(MAP);

  for (std::size_t i=0; i<MAP.size(); i++)
    {
      double ep=1e-3;
      double h;

      Parameters fip(MAP);
      Parameters fin(MAP);
      double ssip;
      double ssin;

      fip[i]=MAP[i]+ep;
      fin[i]=MAP[i]-ep;
      ssip=0.5*SumWeighedSquare(fip);
      ssin=0.5*SumWeighedSquare(fin);
      h=(ssip+ssin-2.0*ss);
      while (h<mindSS)
        {

        }




      Hessian[i][i]=h;
      for (std::size_t j=i+1; j<MAP.size(); j++)
        {
          Parameters fipjp(MAP);
          fipjp[i]=fipjp[i]+ep;
          Parameters fipjn(fipjp);
          fipjp[j]=fipjp[j]+ep;
          fipjn[j]=fipjn[j]-ep;

          Parameters finjp(MAP);
          finjp[i]=finjp[i]-ep;
          Parameters finjn(finjp);
          finjp[j]=finjp[j]+ep;
          finjn[j]=finjn[j]-ep;

          double ssipjp=0.5*SumWeighedSquare(fipjp);
          double ssinjp=0.5*SumWeighedSquare(finjp);
          double ssipjn=0.5*SumWeighedSquare(fipjn);
          double ssinjn=0.5*SumWeighedSquare(finjn);

          Hessian[i][j]=(ssipjp+ssinjn-ssipjn-ssinjp)/ep/ep/4.0;
          Hessian[j][i]=Hessian[i][j];
        }

    }
  Parameters result(MAP);
  result.setCovariance(inv(Hessian));
  return result;
}



Parameters BayesIteration::getEvidence(const Parameters& maximumPostLik, std::size_t num)
{
  double ss0=SumWeighedSquare(maximumPostLik);
  std::vector<Parameters> parvec;
  std::vector<double> parvecVal;
  std::vector<double> m=maximumPostLik.trMeans();
  std::vector<std::vector<double> > covinv=inv(maximumPostLik.getCovariance());
  double sumw=0;
  double sumP=0;

  for (std::size_t i=0; i<num; i++)
    {
      Parameters p;
      double factor=1.0;
      double  dss;
      double sampLogLik=0;
      p=maximumPostLik.randomSample(factor);
      dss=SumWeighedSquare(p)-ss0;
      std::vector<double> d=m;
      for (std::size_t i0=0;i0<p.size(); i0++ )
        {
          d[i0]=p[i0]-m[i0];
        }
      for (std::size_t i0=0;i0<p.size(); i0++ )
        {
          sampLogLik+=d[i0]*d[i0]*covinv[i0][i0];
          for (std::size_t j0=i0+1;j0<p.size(); j0++ )
            {
              sampLogLik+=2*d[i0]*d[j0]*covinv[i0][j0];
            }
        }

      double pv=exp(-0.5*(dss-sampLogLik));
      double plik=0;
      while((pv!=pv)||(pv<1e-2))
        {
          factor/=sqrt(2);
          for (std::size_t i0=0;i0<p.size(); i0++ )
            {
              p[i0]=m[i0]+factor*d[i0];
            }
          dss=SumWeighedSquare(p)-ss0;
          pv=exp(-0.5*(dss-sampLogLik*factor*factor));

        }


      parvec.push_back(p);
      parvecVal.push_back(pv);
      sumw+=pv;
    }

  // now calculate expected mean and covariance
  Parameters result(maximumPostLik);
  std::vector<double> mout(m.size(),0);
  std::vector<std::vector<double> > covout(m.size(),std::vector<double>(m.size(),0));

  for(std::size_t im=0; im<parvec.size(); ++im)
    {
      for (std::size_t i=0; i<parvec[im].size(); i++)
        mout[i]+=parvec[im][i]*parvecVal[im]/sumw;
    }
  for(std::size_t im=0; im<parvec.size(); ++im)
    {
      for (std::size_t i=0; i<parvec[im].size(); i++)
        for (std::size_t j=0; j<parvec[im].size(); j++)
          covout[i][j]+=(parvec[im][i]-mout[i])*(parvec[im][j]-mout[j])*
              parvecVal[im]/sumw;
    }
  result.settMeans(mout);
  result.setCovariance(covout);

  return result;




}




std::map<double,Parameters> BayesIteration::getRandomParameters(const Parameters& per,
                                                                std::size_t num,
                                                                double factor)
{
  std::map<double,Parameters> myMap;
  for (std::size_t i=0; i<num; i++)
    {
      Parameters p=per.randomSample(factor);
      double ss=SumWeighedSquare(p);

      std::cout<<ss<<"\n";
      //    std::cout<<p;
      if (!(ss!=ss))
        myMap[ss]=p;
    }


  for (std::map<double,Parameters>::iterator it=myMap.begin();it!=myMap.end();++it)
    {
      double ss=(*it).first;
      std::cout<<(*it).first<<"\n";

    }


  return myMap;
}

std::map<double,Parameters> BayesIteration::getRandomParameters(std::size_t num,
                                                                double factor)
{
  std::map<double,Parameters> myMap;
  for (std::size_t i=0; i<num; i++)
    {
      Parameters p=priors_.back().randomSample(factor);
      double ss=SumWeighedSquare(p);

      std::cout<<ss<<"\n";
      //    std::cout<<p;
      if (!(ss!=ss))
        myMap[ss]=p;
    }


  for (std::map<double,Parameters>::iterator it=myMap.begin();it!=myMap.end();++it)
    {
      double ss=(*it).first;
      std::cout<<(*it).first<<"\n";

    }


  return myMap;
}


