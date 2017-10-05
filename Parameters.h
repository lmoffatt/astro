#ifndef PARAMETERS
#define PARAMETERS
#include<map>
#include <vector>

#include <string>
#include <cmath>
#include <random>

#include "BaseClass.h"


enum TRANSFORM {LINEAR,LOG,LOGRATIO};

#ifndef PI____
#define PI____
const double PI  =3.141592653589793238463;
#endif

class Transformation
{
public:
  virtual std::string myClass()const=0;
  virtual double eval(double x)const=0;
  virtual double inverse(double x)const=0;
  virtual TRANSFORM myEnum()const=0;

  virtual ~Transformation(){}


};



class Log10Tranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;
  virtual TRANSFORM myEnum()const
  {
    return LOG;
  }

  static std::string ClassName();
};

class Log10RatioTranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;
  virtual TRANSFORM myEnum()const
  {
    return LOGRATIO;
  }

  static std::string ClassName();
};

class LinearTranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;
  virtual TRANSFORM myEnum()const
  {
    return LINEAR;
  }

  static std::string ClassName();
};




class Parameters: public BaseObject
{
public:



  double model()const;

  void setModel(double mod);

  double get(const std::string& name)const;
  //void push_back(const std::string& name, double val,std::string comment);
  void push_back(const std::string& name, double val);

  void push_back_dB(const std::string& name,
                    TRANSFORM tranformation,
                    double meanValue,
                    const std::string &unit,
                    double dBError,
                    const std::vector<double>& Correlations={},
                    const std::string& comment="");

  void push_back_dB(const std::string& name,
                    const std::string &tranformation,
                    double meanValue,
                    const std::string &unit,
                    double dBError,
                    const std::vector<double>& Correlations={},
                    const std::string& comment="");


  void push_back_Hyper_dB(const std::string& name,
                    const std::string &tranformation,
                    double meanValue,
                    const std::string &unit,
                    double dBErrorMean,
                    double mean_log_std_dev,
                    double std_log_std_dev,
                    const std::string& comment="");



  double mean(const std::string& name)const;
  double mean(std::size_t i) const;

  /// pMean("someParameter")=log10(mean("someParameter")
  double tMean(const std::string& name)const;

  double pStd(const std::string& name)const;

  Transformation* getTransform(const std::string& name)const;

  Transformation* getTransform(std::size_t i)const;


  bool hasCovariance()const
  {
    return !cov_.empty();
  }

  /// returns the standard deviation of the logartithm of the parameter in deciBels
  double dBStd(const std::string& name)const;
  double dBStd(std::size_t i)const;

  bool setMeans(const std::vector<std::string> names,const std::vector<double> values);
  bool settMeans(const std::vector<std::string> names,const std::vector<double> log10values);
  bool setpStd(const std::vector<std::string> names,const std::vector<double> values);

  void settMeans(const std::vector<double> log10values);


  std::vector<double> means(const std::vector<std::string> names)const;

  std::vector<double> trMeans()const;

  const std::vector<double> &pStds()const;



  std::vector<std::string> names()const
  {
    return names_;
  }



  std::vector<double> trMeans(const std::vector<std::string> names)const;






  bool setpMean(const std::string& name, double value);

  bool setpStd(const std::string& name, double value);

  bool hasName(const std::string& name)const;

  std::string indexToName(std::size_t i)const;

  std::size_t nameIndex(const std::string& name)const;

  std::vector<std::string> commonNames(const Parameters& other)const;



  ///returns the log of the mean
  const double& operator[](std::size_t)const;
  double& operator[](std::size_t);


  //applies the commonNames mean values to *this
  Parameters& applyParameters(const Parameters& other);

  Parameters& scaleError(double factor);

  void setCovariance(const std::vector< std::vector <double> >& cov);

  const std::vector< std::vector <double> >& getCovariance()const;

  const std::vector< std::vector <double> >& getInvCovariance()const;

  M_Matrix<double> getHessian()const
  {
       M_Matrix<double> out(size(),size(),M_Matrix<double>::DIAGONAL);
       for (std::size_t i =0; i<size(); ++i)
         out[i]=-std::pow(std_of_tr_[i],-2);
       return out;
  }

  M_Matrix<M_Matrix<double>>  getHyperHessian()const
  {
    M_Matrix<M_Matrix<double>>
        out(2,2,M_Matrix<M_Matrix<double>>::SYMMETRIC);

    M_Matrix<double> covinvMean(size(),size(),M_Matrix<double>::DIAGONAL);
    for (std::size_t i=0; i<size(); ++i)
      covinvMean[i]=-std::pow(std_of_tr_[i],-2);
    M_Matrix<double> covinvStd(size(),size(),M_Matrix<double>::DIAGONAL);
    for (std::size_t i=0; i<size(); ++i)
      covinvStd[i]=-std::pow(std_of_random_effects_std_[i],-2);
    out(0,0)=covinvMean;
    out(1,1)=covinvStd;
    return out;

  }

  Parameters randomSample(std::mt19937_64 &mt, double factor=1)const;

  Parameters randomHiperSample(std::mt19937_64 &mt, double factor=1)const;

  M_Matrix<M_Matrix<double>> hyperParameters()const
  {
    M_Matrix<M_Matrix<double>> out(1,2);
    out[0]=M_Matrix<double>(1,size(),mean_of_tr_);
    out[1]=M_Matrix<double>(1,size(),std_of_tr_);
    return out;
  }




  Parameters randomSample(std::mt19937_64 &mt, Parameters prior, double factor)const;


  std::vector<double> randomSampleValues(std::mt19937_64 &mt, Parameters prior, double factor)const;


  Parameters toParameters(const std::vector<double>& o)const
  {
    return Parameters(id(),model_,name_to_i_,names_,o,trans_,unit_);
  }
  Parameters toHyperParameters(const M_Matrix<M_Matrix<double>>& o)const
  {
    auto std_of_tr=o[1].apply([](double x){return std::exp(x);}).toVector();

    return Parameters(id(),model_,name_to_i_,names_,o[0].toVector(),trans_,unit_,std_of_tr);
  }


  double logProb(const Parameters &sample) const;

  M_Matrix<double> logHiperProb(const Parameters &sample)const
  {
    M_Matrix<double> out(1,2);
    out[0]=logProb(sample);
    std::size_t n=sample.size();
    out[1]=0;
    for (std::size_t i=0; i<n; ++i)
      {
        double chi   =sqr(std::log(sample.std_of_tr_[i])-mean_of_random_effects_std_[i])
            /sqr(std_of_random_effects_std_[i]);
        double loglik=-0.5*(chi+std::log(2*PI))
            -std::log(std_of_random_effects_std_[i]);
        out[1]+=loglik;
      }

    return out;

  }

  Parameters(const Parameters& other);
  Parameters();
  ~Parameters(){}

  Parameters& operator=(const Parameters& other);

  void friend swap(Parameters& one, Parameters& other);

  //  friend std::ostream& operator<<(std::ostream& s, const Parameters& p);
  //  friend std::istream& operator>>(std::istream& s, Parameters& p);

  double chi2Distance(const Parameters& one,const Parameters& other)const;

  double chi2Distance(const Parameters &other) const;
  double logDetCov()const;

  unsigned size() const;
  std::string unit(const std::string &name) const;

  const std::vector<double>& mean_of_log_std() const
  {
     return mean_of_random_effects_std_;
  }
  static Transformation* Tr(TRANSFORM tr);

  static TRANSFORM toTr(const std::string& trs);

protected:
  Parameters(const std::string id,
             double model
             ,std::map<std::string, std::size_t> nametoi,
             std::vector<std::string> names,
             std::vector<double> meanoftr,
             std::vector<TRANSFORM> trans,
             std::vector<std::string> unit):
    model_(model)
  ,name_to_i_(nametoi)
  ,names_(names)
  ,mean_of_tr_(meanoftr)
  ,trans_(trans)
  ,unit_(unit)
  ,mean_of_random_effects_std_{}
  ,std_of_random_effects_std_{}
  ,corr_{}
  ,comment_{}
  ,cov_{}
  ,cov_inv_{}
  ,cho_{}
  ,logDetCov_(){
    setId(id);
  }

  Parameters(const std::string id,
             double model
             ,std::map<std::string, std::size_t> nametoi,
             std::vector<std::string> names,
             std::vector<double> meanoftr,
             std::vector<TRANSFORM> trans,
             std::vector<std::string> unit,
             std::vector<double> stdoftr):
    model_(model)
  ,name_to_i_(nametoi)
  ,names_(names)
  ,mean_of_tr_(meanoftr)
  ,trans_(trans)
  ,unit_(unit)
  ,std_of_tr_{stdoftr}
  ,mean_of_random_effects_std_{}
  ,std_of_random_effects_std_{}
  ,corr_{}
  ,comment_{}
  ,cov_{}
  ,cov_inv_{}
  ,cho_{}
  ,logDetCov_(){
    setId(id);
    update();
  }



private:
  std::size_t index(const std::string& name)const;


  std::vector<std::vector<double> >
  buildCovariance(const std::vector<double>& std_of_tr,
                  const std::vector<std::vector<double>> correlations);

  std::vector<std::vector<double> >
  buildCorrelations(const std::vector<double>& std_of_tr,
                    const std::vector<std::vector<double>> covariance);

  std::vector<double> getStdDev(const std::vector<std::vector<double>>& covariance);

  void update() override;

  static std::map<TRANSFORM,Transformation*> tr_map_;

  /// this variables


  double model_;

  std::map<std::string, std::size_t> name_to_i_;

  std::vector<std::string> names_;

  std::vector<double> mean_of_tr_;
  std::vector<TRANSFORM> trans_;
  std::vector<std::string> unit_;
  std::vector<double> std_of_tr_;

  std::vector<double> mean_of_random_effects_std_;
  std::vector<double> std_of_random_effects_std_;


  std::vector<std::vector<double>> corr_;

  std::vector<std::string> comment_;



  std::vector< std::vector <double> > cov_;
  std::vector< std::vector <double> > cov_inv_;

  std::vector< std::vector <double> > cho_;
  double logDetCov_;





  // BaseClass interface
public:
  static std::string ClassName(){return "parameters";}
  virtual std::string myClass() const override {return ClassName();}

  // BaseObject interface
public:
  virtual Parameters *create() const override
  {
    return new Parameters;
  }
  virtual std::ostream &writeBody(std::ostream &s) const override;
  virtual void clear() override;
  virtual bool readBody(std::string &line, std::istream &s, std::ostream& logs) override;


  std::ostream& writeHeaderDataFrame(std::ostream& os)const
    {
      for (std::size_t i=0; i<size(); ++i)
        {
          os<<"Value...paramName.."<<indexToName(i)<<"\t";
          os<<"TrasnfValue...transfParm.."<<Tr(trans_[i])->myClass()<<indexToName(i);
          if (i+1<size()) os<<"\t";

        }
          return os;
  }

  std::ostream& writeRowDataFrame(std::ostream& os)const
  {
    for (std::size_t i=0; i<size(); ++i)
      {
        os<<mean(i)<<"\t";
        os<<(*this)[i];
        if (i+1<size()) os<<"\t";
      }
    return os;
  }


  void writeValueInSingleRowHeader(std::ostream& os)const
  {
    for (std::size_t i=0; i<size(); ++i)
      os<<indexToName(i)<<"\t";
    for (std::size_t i=0; i<size(); ++i)
      os<<Tr(trans_[i])->myClass()<<indexToName(i)<<"\t";
  }


  void writeValueInSingleRow(std::ostream& os)const
  {
    for (std::size_t i=0; i<size(); ++i)
      os<<mean(i)<<"\t";
    for (std::size_t i=0; i<size(); ++i)
      os<<(*this)[i]<<"\t";
  }



  void writeDataFrameHeader(std::ostream& os)const
  {
    writeDataFrameHeaderRow(os);
    writeDataFrameHeaderRow(os,"_i");
    writeDataFrameHeaderRow(os,"_i","_j");

  }



  void writeDataFrameRow(std::ostream& os)const
  {
    os<<id()<<"\t";
    os<<model()<<"\t";
    os<<logDetCov()<<"\t";
  }

  void writeDataFrameRow(std::ostream& os,
                         std::size_t ipar) const
  {
    os<<indexToName(ipar)<<"\t";
    os<<unit_[ipar]<<"\t";
    os<<Tr(trans_[ipar])->myClass()<<"\t";
    os<<comment_[ipar]<<"\t";
    os<<mean(ipar)<<"\t";
    os<<(*this)[ipar]<<"\t";
    os<<pStds()[ipar]<<"\t";
  }
  void writeDataFrameRow(std::ostream& os,
                         std::size_t ipar,
                         std::size_t jpar) const
  {
    os<<cov_[ipar][jpar]<<"\t";
    os<<corr_[ipar][jpar]<<"\t";
    os<<cov_inv_[ipar][jpar]<<"\t";
    os<<cho_[ipar][jpar]<<"\t";

  }
  void writeDataFrameHeaderRow(std::ostream& os) const
  {
    os<<"id"<<"\t";
    os<<"model"<<"\t";
    os<<"logDetCov"<<"\t";
  }

  void writeDataFrameHeaderRow(std::ostream& os,
                               const std::string s) const
  {
    os<<"parameterName"<<s<<"\t";
    os<<"unit"<<s<<"\t";
    os<<"transformation"<<s<<"\t";
    os<<"comment"<<s<<"\t";
    os<<"value"<<s<<"\t";
    os<<"transValue"<<s<<"\t";
    os<<"stddev"<<s<<"\t";
  }
  void writeDataFrameHeaderRow(std::ostream& os,
                               std::string ipar,
                               std::string jpar) const
  {
    os<<"cov"<<ipar<<jpar<<"\t";
    os<<"corr"<<ipar<<jpar<<"\t";
    os<<"cov_inv"<<ipar<<jpar<<"\t";
    os<<"cholesky"<<ipar<<jpar<<"\t";

  }

  void writeDataFrame(std::ostream& os) const
  {
    writeDataFrameHeader(os);
    os<<"\n";

    for (std::size_t ipar=0; ipar<size(); ++ipar)
      for (std::size_t jpar=0; jpar<size(); ++jpar)
        {
          writeDataFrameRow(os);
          writeDataFrameRow(os,ipar);
          writeDataFrameRow(os,ipar,jpar);
          os<<"\n";
        }
    os<<"\n";

  }


};





//std::ostream& operator<<(std::ostream& s, const Parameters& p);
//std::istream& operator>>(std::istream& s, Parameters& p);

double dbDistance(const Parameters& one,const Parameters& other);

bool areTheSame(const Parameters& one, const Parameters& other);


double randNormal(std::mt19937_64 &mt, double mean, double stddev);
double randNormal(std::mt19937_64 &mt);





#endif // PARAMETERS

