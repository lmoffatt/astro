#ifndef PARAMETERS
#define PARAMETERS
#include<map>
#include <vector>

#include <string>

class Parameters
{
public:
  double get(const std::string& name)const
  {
    return mean(name);
  }
  void push_back(const std::string& name, double val,std::string comment);
  void push_back(const std::string& name, double val);


  void read(std::string &line, std::istream &s);

  std::ostream& write(std::ostream& s)const;

  double mean(const std::string& name)const;

  /// pMean("someParameter")=log10(mean("someParameter")
  double pMean(const std::string& name)const;

  double mean_ratio(const std::string& name)const;



  /// returns the standard deviation of the logartithm of the parameter
  double pStd(const std::string& name)const;

  /// returns the standard deviation of the logartithm of the parameter in deciBels
  double dBStd(const std::string& name)const;


  double residual(const std::string& name, double value)const;
  double pResidual(const std::string& name, double log10value)const;

  std::string mode()const;
  Parameters& setMode(const std::string& mode_);

  std::vector<double> residuals(const std::vector<std::string> names,const std::vector<double> log10Values);

  std::vector<double> residuals(const Parameters& sample)const;



  bool setMeans(const std::vector<std::string> names,const std::vector<double> values);
  bool setpMeans(const std::vector<std::string> names,const std::vector<double> log10values);
  bool setpStd(const std::vector<std::string> names,const std::vector<double> values);

  void setpMeans(const std::vector<double> log10values);


  std::vector<double> means(const std::vector<std::string> names)const;

  std::vector<double> pMeans()const;

  std::vector<double> pStds()const;

  std::vector<std::string> names()const;



  std::vector<double> pMeans(const std::vector<std::string> names)const;



  void push_back_dB(const std::string& name,
                    double meanValue,
                    const std::string &unit,
                    double dBStdValue,
                    const std::string& comment="");

  void push_back_1S(const std::string& name
                    , double minValue_p34
                    , const std::string &unit
                    , double maxValue_p68
                    , const std::string &comment="");

  void push_back_2S(const std::string& name
                    , double minValue_p02
                    , const std::string& unit
                    , double maxValue_p98
                    , const std::string& comment="");

  void push_back_3S(const std::string& name
                    ,double minValue_p001
                    ,const std::string& unit
                    , double maxValue_p999
                    , const std::string& comment="");


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

  std::vector< std::vector <double> > getCovariance()const;


  Parameters randomSample(double factor=1.0)const;




  Parameters randomSample(Parameters prior,double factor,double probIncludeParameter)const;
  Parameters randomSample(Parameters prior,double factor)const;



  Parameters(const Parameters& other);
  Parameters();
  ~Parameters(){}

  Parameters& operator=(const Parameters& other);

  void friend swap(Parameters& one, Parameters& other);

//  friend std::ostream& operator<<(std::ostream& s, const Parameters& p);
//  friend std::istream& operator>>(std::istream& s, Parameters& p);

  double chi2Distance(const Parameters& other)const;



  unsigned size() const;
  std::string unit(const std::string &name) const;
private:
  std::map<std::string, std::size_t> name_;
  std::vector<double> pMean_;
  std::vector<double> pStd_;  // not in dB
  std::vector< std::vector <double> > cov_;
  std::vector< std::vector <double> > cho_;

  std::string mode_;

  std::vector<std::string> unit_;
  std::vector<std::string> comment_;


};

//std::ostream& operator<<(std::ostream& s, const Parameters& p);
//std::istream& operator>>(std::istream& s, Parameters& p);

double dbDistance(const Parameters& one,const Parameters& other);

bool areTheSame(const Parameters& one, const Parameters& other);


double randNormal(double mean,double stddev);
double randNormal();


#endif // PARAMETERS

