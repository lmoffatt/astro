#ifndef PARAMETERST_H
#define PARAMETERST_H

#include "Matrix.h"
#include "myTuples.h"
#include "DirectSum.h"




template< typename...Fields>
class ParametersT
{
public:
  typedef TypeCo<Fields...> fields;
  typedef std::tuple<typename Fields::field_value_type...> tuData;
 // typedef ParametersT<Fields...> field_value_type;


private:
  tuData data;

public:
  ParametersT( tuData d):
    data{d}{}


  template<typename... F2>
  ParametersT(ParametersT<F2...> p): data(p.data){}



  ParametersT():
   data{}{}


  template <typename FieldName>
  decltype(auto) get()
  {
    constexpr std::size_t I=index_of<FieldName>(fields());
    return std::get<I>(data);
  }


  template <typename FieldName>
  decltype(auto) get()const
  {
    constexpr std::size_t I=index_of<FieldName>(fields());
    return std::get<I>(data);
  }

  template <class Op,std::size_t...Is>
  decltype(auto)
   apply(std::index_sequence<Is...>)const
  {
      return Op::apply(std::get<Is>(data)...);
  }

  template <class Op>
  decltype(auto)
  apply() const
  {
    return apply(std::index_sequence_for<Fields...>());
  }


  template <std::size_t...Is>
  std::size_t size(std::index_sequence<Is...>)const
  {
      return sum_Is(mySize(std::get<Is>(data))...);
  }

  std::size_t size()const
  {
    return size(std::index_sequence_for<Fields...>());
  }


  template <typename F>
  constexpr const char* fieldNames(int)
  {
    return F::ClassName();
  }
  template <typename F,typename...Fs>
  constexpr const char* fieldNames(long)
  {
    return F::ClassName()+","+fieldNames<Fs...>(0);
  }


  constexpr const char* ClassName(){
    return "{"+fieldNames<Fields...>(0)+"}";
  };








  ParametersT<Fields...>& operator+=(const ParametersT<Fields...> other)
  {
    data+=other.data;
    return *this;
  }


  template< typename...F0>
  ParametersT<Fields...,F0...>
  cat(ParametersT<F0...> p1)
  {
    return ParametersT<Fields...,F0...>(std::tuple_cat(data, p1.data));
  }




};
template< typename...Fields>
ParametersT<Fields...> operator+(const ParametersT<Fields...>& one, const ParametersT<Fields...>& other)
{
   ParametersT<Fields...> out(one);
   out+=other;
   return out;
}

template< typename...Fields>
ParametersT<Fields...> operator+(ParametersT<Fields...> one, const ParametersT<Fields...>& other)
{
   one+=other;
   return std::move(one);
}


template <typename Field, typename... Fields>
decltype(auto) get(const ParametersT<Fields...>& me )
{
  return  me.template get<Field>();
}
template <typename Field, typename... Fields>
ParametersT<Field> getP(const ParametersT<Fields...>& me )
{
  return  ParametersT<Field>({me.template get<Field>()});
}



template< typename P>
decltype(auto)
cats_impl(P p)
{
  return p;
}

template< typename P,typename...Ps>
decltype(auto)
cats_impl(P p,Ps...ps)
{
  return p.cat(cats_impl(ps...));
}


template< typename P,typename...Ps>
decltype(auto)
cats(P p,Ps...ps)
{
  return p.cat(cats_impl(ps...));
}



template< typename...Ps>
decltype(auto)
cats(Ps...ps)
{
  return cats_impl(ps...);
}






#endif // PARAMETERST_H
