#include "Parameters.h"
#include <limits>
#include <cmath>
#include <cstdlib>
#include <istream>
#include <sstream>
#include <iostream>
#include "MatrixInverse.h"

#include <random>



std::string Log10Tranformation::myClass()const
{
  return ClassName();
}
double Log10Tranformation::eval(double x)const
{
  return std::log10(x);
}
double Log10Tranformation::inverse(double x)const
{
  return std::pow(10.0,x);
}

std::string Log10Tranformation::ClassName()
{
  return "Log10";
}


std::string Log10RatioTranformation::myClass()const
{
  return ClassName();
}
double Log10RatioTranformation::eval(double x)const
{
  return std::log10(x/(1.0-x));
}
double Log10RatioTranformation::inverse(double x)const
{
  double r= std::pow(10.0,x);
  return r/(r+1.0);
}

std::string Log10RatioTranformation::ClassName()
{
  return "Log10Ratio";
}



std::string LinearTranformation::myClass()const
{
  return ClassName();
}
double LinearTranformation::eval(double x)const
{
  return x;
}
double LinearTranformation::inverse(double x)const
{
  return x;
}

std::string LinearTranformation::ClassName()
{
  return "Linear";
}














Parameters::Parameters(const Parameters& other):
  id_(other.id_)
  ,model_(other.model_)
,name_to_i_(other.name_to_i_)
, names_(other.names_)
,mean_of_tr_(other.mean_of_tr_)
,trans_(other.trans_)
,unit_(other.unit_)
,std_of_tr_(other.std_of_tr_)
,corr_(other.corr_)
,comment_(other.comment_)
,cov_(other.cov_)
,cov_inv_(other.cov_inv_)
,cho_(other.cho_)
,logDetCov_(other.logDetCov_){}

bool Parameters::readBody(std::string &line, std::istream &s)
{
  std::string name;
  std::stringstream ss(line);
  ss.str(line);
  ss.clear();
  ss>>name;
  while (name.empty()
         &&safeGetline(s,line))
    {
      ss.clear();
      ss.str(line);
      ss>>name;
    }
  if (name=="model")
    {
      double mod;
      ss>>mod;
      setModel(mod);
      safeGetline(s,line);
    }
  while (true)
    {
      double val, std_in_dB;
      std::string transformation,unit, db,comment;
      std::vector<double> corrCoef;
      ss.str(line);
      name.clear();
      ss.clear();
      ss>>name;
      while (name.empty()
             &&safeGetline(s,line))
        {
          ss.str(line);
          ss.clear();
          ss>>name;
        }
      if (name.empty()||name.find("--",0)!=name.npos)
        break;
      else
        {

          ss>>transformation>>val;
          char ch;
          while ((ss>>ch)&&(ch!='[')) {}
          if (ch=='[')
            {
              while ((ss>>ch)&&(ch!=']'))
                unit.push_back(ch);
            }
          else unit="";
          if (ss>>std_in_dB)
            {
              while ((ss>>ch)&&(ch!='[')) {}
              if (ch=='[')
                while ((ss>>ch)&&(ch!=']'))
                  db.push_back(ch);
              else
                db="";
              double corrcoef;
              while (ss>>corrcoef)
                corrCoef.push_back(corrcoef);
              ss.clear();
              while ((ss>>ch)&&(ch!='/')) {}
              if ((ss>>ch)&&(ch=='/'))
                std::getline(ss,comment);
              else
                comment="";
            }
          else
            {
              std_in_dB=0;
              db="";
              comment="";
            }
          push_back_dB(name,transformation,val,unit,std_in_dB,corrCoef,comment);
          line.clear();
        }
    }
  line.clear();
  return true;
}


Parameters& Parameters::operator=(const Parameters& other)
{
  if (this!=&other)
    {
      Parameters tmp(other);
      swap(*this,tmp);
    }
  return *this;
}

void swap(Parameters& one, Parameters& other)
{
  std::swap(one.id_,other.id_);
  std::swap(one.model_,other.model_);
  std::swap(one.name_to_i_,other.name_to_i_);
  std::swap(one.names_,other.names_);

  std::swap(one.mean_of_tr_,other.mean_of_tr_);
  std::swap(one.std_of_tr_,other.std_of_tr_);

  std::swap(one.trans_,other.trans_);
  std::swap(one.corr_,other.corr_);
  std::swap(one.logDetCov_,other.logDetCov_);
  std::swap(one.cov_inv_,other.cov_inv_);

  std::swap(one.cov_,other.cov_);
  std::swap(one.cho_,other.cho_);
  std::swap(one.unit_,other.unit_);
  std::swap(one.comment_,other.comment_);



}



double Parameters::model() const
{
  return model_;
}

void Parameters::setModel(double mod)
{
  model_=mod;
}

double Parameters::get(const std::string &name) const
{
  return mean(name);
}

void Parameters::push_back(const std::string &name, double val)
{
  push_back_dB(name,"<LINEAR>",val,"",0);
}



std::ostream &Parameters::writeBody(std::ostream &s) const
{
  s<<"model\t"<<model()<<"\n\n";
  for (std::size_t i=0; i<size(); i++)
    {
      s<<names_[i]<<"\t";
      s<<"<"<<Tr(trans_[i])->myClass()<<">\t";
      s<<Tr(trans_[i])->inverse(mean_of_tr_[i])<<"\t";
      s<<"["+unit_[i]+"]"<<"\t"<<std_of_tr_[i]*10<<"[dB]\t";
      if (!corr_.empty())
        {
          for (std::size_t j=0; j<corr_[i].size(); ++j)
            s<<corr_[i][j]<<"\t";
        }
      s<<"//"<<comment_[i]<<"\n";
    }

  return s;
}

void Parameters::clear()
{
  name_to_i_.clear();

    names_.clear();

    mean_of_tr_.clear();
    trans_.clear();
    unit_.clear();
    std_of_tr_.clear();

    corr_.clear();

    comment_.clear();



    cov_.clear();
    cov_inv_.clear();

    cho_.clear();
   logDetCov_=0;
}

double Parameters::mean(const std::string& name)const
{
  std::size_t i=index(name);
  if (i!=std::string::npos)
    {
      double m=Tr(trans_[i])->inverse(mean_of_tr_[i]);
      return m;
    }
  else
    return std::numeric_limits<double>::quiet_NaN();
}

double Parameters::mean(std::size_t i)const
{
  if (i!=std::string::npos)
    {
      double m=Tr(trans_[i])->inverse(mean_of_tr_[i]);
      return m;
    }
  else
    return std::numeric_limits<double>::quiet_NaN();
}


double Parameters::tMean(const std::string& name)const
{
  std::size_t i=index(name);
  if (i!=std::string::npos)
    {
      double m=mean_of_tr_[i];
      return m;
    }
  else
    return std::numeric_limits<double>::quiet_NaN();

}

double Parameters::pStd(const std::string &name) const
{
  std::size_t i=index(name);
  if (i!=std::string::npos)
    {
      double s=std_of_tr_[i];
      return s;
    }
  else
    return std::numeric_limits<double>::quiet_NaN();

}

Transformation* Parameters::getTransform(const std::string &name) const
{
  std::size_t i=index(name);
  if (i!=std::string::npos)
    {
      return Tr(trans_[i]);
    }
  else
    return Tr(LINEAR);

}

Transformation* Parameters::getTransform(std::size_t i) const
{
  if (i!=std::string::npos)
    {
      return Tr(trans_[i]);
    }
  else
    return Tr(LINEAR);

}


std::string Parameters::unit(const std::string& name)const
{
  std::map<std::string,size_t>::const_iterator it=name_to_i_.find(name);
  if((it!=name_to_i_.end())&& !std_of_tr_.empty())
    return unit_[(*it).second];
  else
    return "";

}

Transformation *Parameters::Tr(TRANSFORM tr)
{
  return tr_map_[tr];
}

TRANSFORM Parameters::toTr(const std::string &trs)
{
  if ((trs=="<LOG>")||(trs=="<Log10>"))
    return LOG;
  else if ((trs=="<LOGRATIO>")||(trs=="<Log10Ratio>"))
    return LOGRATIO;
  else if (trs=="<LINEAR>")
    return LINEAR;
  else return LINEAR;
}

std::size_t Parameters::index(const std::string &name) const
{
  std::map<std::string,std::size_t>::const_iterator it=name_to_i_.find(name);
  if (it!=name_to_i_.end())
    return it->second;
  else
    return std::string::npos;

}

std::map<TRANSFORM,Transformation*> Parameters::tr_map_=
{{LINEAR,new LinearTranformation()},
 {LOG, new Log10Tranformation()},
 {LOGRATIO,new Log10RatioTranformation()}
};





std::vector<std::vector<double> > Parameters::buildCovariance(const std::vector<double> &std_of_tr, const std::vector<std::vector<double> > correlations)
{
  std::size_t n=std_of_tr.size();
  if (n!=correlations.size())
    return {};
  else
    {
      std::vector<std::vector<double>> o(n,std::vector<double>(n));
      for (std::size_t i=0; i<n; ++i)
        {
          o[i][i]=std_of_tr[i]*std_of_tr[i];
          for (std::size_t j=0; j<i ; ++j)
            {
              double s2_ij=correlations[i][j]*std_of_tr[i]*std_of_tr[j];
              o[i][j]=s2_ij;
              o[j][i]=s2_ij;
            }
        }
      return o;
    }
}

std::vector<std::vector<double> > Parameters::buildCorrelations(const std::vector<double> &std_of_tr, const std::vector<std::vector<double> > covariance)
{
  std::size_t n=std_of_tr.size();
  if (n!=covariance.size())
    return {};
  else
    {
      std::vector<std::vector<double>> o(n);
      for (std::size_t i=0; i<n; ++i)
        {
          o[i]=std::vector<double>(i+1);
          o[i][i]=1;
          for (std::size_t j=0; j<i ; ++j)
            {
              double r_ij=covariance[i][j]/std_of_tr[i]/std_of_tr[j];
              o[i][j]=r_ij;
            }
        }
      return o;
    }
}

std::vector<double> Parameters::getStdDev(const std::vector<std::vector<double> > &covariance)
{
  std::size_t n=covariance.size();

  std::vector<double> o(n);
  for (std::size_t i=0; i<n; ++i)
    {
      o[i]=std::sqrt(covariance[i][i]);
    }
  return o;
}

void Parameters::update()
{
  if (corr_[1].size()==0)
    {
      cov_=std::vector<std::vector<double>>(size(),std::vector<double>(size(),0));
      for (std::size_t i=0; i<size(); ++i)
        cov_[i][i]=std_of_tr_[i]*std_of_tr_[i];

      cov_inv_=inv(cov_);

      cho_=chol(cov_);
      for (std::size_t i=0; i< size();++i)
        logDetCov_+=2*log(std_of_tr_[i]);
    }
  else
    {
      cov_=buildCovariance(std_of_tr_,corr_);
      cov_inv_=inv(cov_);
      cho_=chol(cov_);
      for (std::size_t i=0; i< size();++i)
        logDetCov_+=log(cho_[i][i]);
    }

}


/// returns the standard deviation in dB (deciBel)
double Parameters::dBStd(const std::string& name)const
{
  std::map<std::string,size_t>::const_iterator it=name_to_i_.find(name);
  if((it!=name_to_i_.end())&& !std_of_tr_.empty())
    return std_of_tr_[(*it).second]*10;
  else
    return std::numeric_limits<double>::quiet_NaN();

}


/// returns the standard deviation in dB (deciBel)
double Parameters::dBStd(std::size_t i)const
{
   return std_of_tr_[i]*10;
}



std::size_t Parameters::nameIndex(const std::string& name)const
{
  return (*name_to_i_.find(name)).second;

}



bool Parameters::setMeans(const std::vector<std::string> names,const std::vector<double> values)
{
  for (std::size_t i=0;i<names.size();i++)
    {
      if (hasName(names[i]))
        mean_of_tr_[nameIndex(names[i])]=getTransform(names[i])->eval(values[i]);
      else return false;
    }
  return true;

}
bool Parameters::settMeans(const std::vector<std::string> names,const std::vector<double> log10values)
{
  for (std::size_t i=0;i<names.size();i++)
    {
      if (hasName(names[i]))
        mean_of_tr_[nameIndex(names[i])]=log10values[i];
      else return false;
    }
  return true;

}

void Parameters::settMeans(const std::vector<double> trasf_values)
{
  mean_of_tr_=trasf_values;
}


bool Parameters::setpStd(const std::vector<std::string> names,const std::vector<double> values)
{
  for (std::size_t i=0;i<names.size();i++)
    {
      if (hasName(names[i]))
        std_of_tr_[nameIndex(names[i])]=values[i];
      else return false;
    }
  return true;

}


std::vector<double> Parameters::means(const std::vector<std::string> names)const
{
  std::vector<double> result(names.size());
  for (std::size_t i=0;i<names.size();i++)
    {
      result[i]=mean(names[i]);
    }
  return result;


}

std::vector<double> Parameters::trMeans(const std::vector<std::string> names)const
{
  std::vector<double> result(names.size());
  for (std::size_t i=0;i<names.size();i++)
    {
      result[i]=tMean(names[i]);
    }
  return result;
}


void Parameters::push_back_dB(const std::string& name,
                              const std::string &tranformation,
                              double meanValue,
                              const std::string& unit,
                              double dBError,
                              const std::vector<double> &Correlations,
                              const std::string& comment)
{
  push_back_dB(name,toTr(tranformation),meanValue,unit,dBError,Correlations,comment);
}


void Parameters::push_back_dB(const std::string& name,
                              TRANSFORM tranformation,
                              double meanValue,
                              const std::string& unit,
                              double dBError,
                              const std::vector<double>& Correlations,
                              const std::string& comment)
{
  name_to_i_[name]=mean_of_tr_.size();
  names_.push_back(name);
  TRANSFORM tr=tranformation;
  trans_.push_back(tr);
  double t=Tr(tr)->eval(meanValue);
  mean_of_tr_.push_back(t);
  std_of_tr_.push_back(dBError/10);
  corr_.push_back(Correlations);
  unit_.push_back(unit);
  comment_.push_back(comment);
}


bool Parameters::setpMean(const std::string& name, double value)
{
  std::map<std::string,std::size_t>::iterator it=name_to_i_.find(name);
  if(it!=name_to_i_.end())
    {
      mean_of_tr_[(*it).second]=value;
      return true;
    }
  else
    return false;

}

bool Parameters::setpStd(const std::string& name, double value)
{
  std::map<std::string,std::size_t>::iterator it=name_to_i_.find(name);
  if(it!=name_to_i_.end())
    {
      std_of_tr_[(*it).second]=value;
      return true;
    }
  else
    return false;

}

bool Parameters::hasName(const std::string& name)const{
  return name_to_i_.find(name)!=name_to_i_.end();
}


std::string Parameters::indexToName(std::size_t i)const
{
  if (i<names_.size())
    return names_[i];
  else
    return "";

}

std::vector<std::string> Parameters::commonNames(const Parameters& other)const
{
  std::vector<std::string>  myCommonNames;
  for (std::map<std::string,std::size_t>::const_iterator it=name_to_i_.begin();
       it!=name_to_i_.end();
       ++it)
    {
      if (other.hasName(it->first))
        myCommonNames.push_back(it->first);

    }
  return myCommonNames;

}



Parameters Parameters::randomSample(std::mt19937_64& mt,double factor)const
{


  Parameters sample(*this);
  if (cho_.empty())
    {
      for (std::size_t i=0; i<mean_of_tr_.size();i++)

        {

          sample.mean_of_tr_[i]=randNormal(mt,mean_of_tr_[i],std_of_tr_[i]*factor);
          sample.std_of_tr_[i]=0;
        }
    }
  else
    {
      std::vector<double> z(mean_of_tr_.size());
      for (std::size_t i=0; i<mean_of_tr_.size();i++)
        {
          z[i]=randNormal(mt)*factor;
        }
      for (std::size_t i=0; i<mean_of_tr_.size();i++)
        {

          sample.mean_of_tr_[i]=mean_of_tr_[i];
          sample.std_of_tr_[i]=0;
          for (std::size_t j=0; j<i+1;j++)
            {
              sample.mean_of_tr_[i]+=cho_[i][j]*z[j];
            }
        }
    }
  return sample;

}


Parameters Parameters::randomSample(std::mt19937_64& mt,Parameters prior,double factor)const

{
  Parameters sample(*this);
  if (prior.cho_.empty())
    {
      for (std::size_t i=0; i<mean_of_tr_.size();i++)

        {

          sample.mean_of_tr_[i]=randNormal(mt,mean_of_tr_[i]
                                           ,prior.std_of_tr_[i]*factor);
          sample.std_of_tr_[i]=0;
        }
    }
  else
    {
      std::vector<double> z(mean_of_tr_.size());
      for (std::size_t i=0; i<mean_of_tr_.size();i++)
        {
          z[i]=randNormal(mt)*factor;
        }
      for (std::size_t i=0; i<mean_of_tr_.size();i++)
        {

          sample.mean_of_tr_[i]=mean_of_tr_[i];
          sample.std_of_tr_[i]=0;
          for (std::size_t j=0; j<i+1;j++)
            {
              sample.mean_of_tr_[i]+=prior.cho_[i][j]*z[j];
            }
        }
    }
  Parameters s(sample);
  return s;

}

std::vector<double> Parameters::randomSampleValues(std::mt19937_64 &mt, Parameters prior, double factor) const

{
  std::vector<double> o(size());
  if (prior.cho_.empty())
    {
      for (std::size_t i=0; i<mean_of_tr_.size();i++)

        {

          o[i]=randNormal(mt,mean_of_tr_[i]
                                           ,prior.std_of_tr_[i]*factor);
        }
    }
  else
    {
      std::vector<double> z(mean_of_tr_.size());
      for (std::size_t i=0; i<mean_of_tr_.size();i++)
        {
          z[i]=randNormal(mt)*factor;
        }
      for (std::size_t i=0; i<mean_of_tr_.size();i++)
        {

          o[i]=mean_of_tr_[i];
          for (std::size_t j=0; j<i+1;j++)
            {
              o[i]+=prior.cho_[i][j]*z[j];
            }
        }
    }
  return o;
}

double Parameters::logProb(const Parameters& sample)const
{
  std::size_t n=sample.size();
   std::vector<double> d(n);
  for (std::size_t i=0; i<n; ++i)
    d[i]=sample.mean_of_tr_[i]-this->mean_of_tr_[i];

  double z=0;
  for (std::size_t i=0; i<n; ++i)
    for (std::size_t j=0; j<n; ++j)
      z+=d[i]*this->cov_inv_[i][j]*d[j];

  return -0.5*(logDetCov_+z+n*log(2*PI));

}













unsigned Parameters::size()const
{
  return name_to_i_.size();
}




double randNormal(std::mt19937_64& mt,double mean,double stddev)
{
  std::normal_distribution<> nd(mean,stddev);
  return nd(mt);
}
double randNormal(std::mt19937_64& mt)
{

  std::normal_distribution<> nd(0,1);
  return nd(mt);
}

Parameters& Parameters::applyParameters(const Parameters& other)
{
  std::vector<std::string> cnames=commonNames(other);
  settMeans(cnames,other.trMeans(cnames));
  return *this;
}

Parameters& Parameters::scaleError(double factor)
{
  for (std::size_t i=0; i<std_of_tr_.size(); i++)
    {
      std_of_tr_[i]*=factor;
    }
  return *this;
}

void Parameters::setCovariance(const std::vector< std::vector <double> >& cov)
{
  if (cov.size()==size())
    {
      cov_=cov;
      std_of_tr_=getStdDev(cov);
      corr_=buildCorrelations(std_of_tr_,cov);
      cho_=chol(cov);
      cov_inv_=inv(cov);
      logDetCov_=0;
      for (std::size_t i=0; i<size(); ++i)
        logDetCov_+=log(cho_[i][i]);
      logDetCov_*=2;
    }
}



const std::vector<std::vector<double> > &Parameters::getCovariance()const
{
  return cov_;
}

const std::vector<std::vector<double> > &Parameters::getInvCovariance() const
{
  return cov_inv_;
}






std::vector<double> Parameters::trMeans()const
{
  return mean_of_tr_;
}





const std::vector<double>& Parameters::pStds()const
{
  return std_of_tr_;
}


const double& Parameters::operator[](std::size_t i)const
{
  return mean_of_tr_[i];
}
double& Parameters::operator[](std::size_t i){
  return mean_of_tr_[i];
}





Parameters::Parameters(){}


bool areTheSame(const Parameters& one, const Parameters& other)
{

  if (one.size()!=other.size())
    return false;
  for (std::size_t i=0; i<one.size();i++)
    if (one.indexToName(i)!=other.indexToName(i))
      return false;
  return true;
}


double dbDistance(const Parameters& one,const Parameters& other)
{
  double result=0;

  for (std::size_t i=0;i<one.size();i++)
    {
      result+=pow(other.tMean(other.indexToName(i))-one.tMean(one.indexToName(i)),2);
    }
  result=sqrt(result/one.size())*10;
  return result;
}


double Parameters::chi2Distance(const Parameters &one, const Parameters &other)const
{
  double result=0;

  for (std::size_t i=0;i<size();i++)
    {
       double e=one.trMeans()[i]-other.trMeans()[i];
       double s=std_of_tr_[i];
       result+=std::pow(e/s,2)/2.0;
    }
  // result=sqrt(result/one.size())*10;
  return result;

}

double Parameters::chi2Distance(const Parameters &other) const
{
  return chi2Distance(*this,other);

}

double Parameters::logDetCov() const
{
  return logDetCov_;
}

