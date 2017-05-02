#ifndef FUNCTIONT_H
#define FUNCTIONT_H
#include "ParametersT.h"

namespace fn
{


template<class Field,typename T>
struct field
{
  typedef ParametersT<Field> domain_type;
  typedef ParametersT<Field> codomain_type;
  template<class P>
  static
  codomain_type
  apply(const P& p)
  {
    return getP<Field>(p);
  }
};


template <class... fs>
struct bundle
{
private:
  template<typename...>
  struct fn_impl_0;

  template<typename codomain_field>
  struct fn_impl_0<codomain_field>
  {};

  template<typename codomain_field, class fn,class... fns>
  struct fn_impl_0<codomain_field, fn,fns...>
  {
    typedef
    typename
    std::conditional<std::is_same
    <typename fn::codomain_type, typename codomain_field::codomain_type>::value,
    fn,
    fn_impl_0<codomain_field,fs...>>
    ::type type;

  };
  template<typename codomain_field>
  struct fn_impl
  {
    typedef
    typename
    fn_impl_0<codomain_field,fs...>
    ::type type;
  };

public:
  typedef direct_sums<typename fs::domain_type...> domain_type;

//  typedef ParametersT<typename fs::codomain_type...> codomain_type;

  template <typename codomain_field>
  using fn=typename fn_impl<codomain_field>::type;

  template<class P>
  static
  decltype(auto)
  apply(const P& p)
  {
    return cats(fs::apply(p)...);
  }
};



template<class f, class g>
struct compose
{

  static_assert(std::is_same<typename g::codomain_type,typename f::domain_type>::value,"g::codomain_type!=typename f::domain_type" );
   typedef typename g::domain_type domain_type;
  typedef typename f::codomain_type codomain_type;

  template<class P>
  static
  codomain_type
  apply(const P& p)
  {
      return f::apply(g::apply(p));
  }

};

template <typename Field, class f>
struct assign
{
  typedef typename f::domain_type domain_type;
 typedef typename Field::codomain_type codomain_type;

 template<class P>
 static
 codomain_type
 apply(const P& p)
 {
     return f::apply(p);
 }
};


struct sum_imp
{
  template<typename T,typename S>
  static   decltype(auto)
  op(T x, S y){return x+y;}

  constexpr const char* ClassName(){return "sum";};


};

struct prod_imp
{
  template<typename T,typename S>
  static   decltype(auto)
  op(T x, S y){return x*y;}

  constexpr const char* ClassName(){return "prod";};


};



struct minus_imp
{
  template<typename T>
  static   decltype(auto)
  op(T x){return -x;}
  constexpr const char* ClassName(){return "minus";};

};

struct exp_imp
{
  template<typename T>
  static   decltype(auto)
  op(T x){return std::exp(x);}
  constexpr const char* ClassName(){return "exp ";};

};





template<class Op,class... fs>
struct NOp_imp
{
private:
  template<typename T>
  struct Nop_op;
  template<typename... Fields>
  struct Nop_op<ParametersT<Fields...>>
  {

  private:
    template<typename P,typename F>
    static
    decltype(auto)
    apply_imp(const P& p, TypeCo<F>, int)
    {
      return get<F>(p);
    }


    template<typename P,typename F,typename... Fs>
    static
    decltype(auto)
    apply_imp(const P& p, TypeCo<F,Fs...>,long)
    {
      return Op::apply(get<F>(p),apply_imp(p, TypeCo<Fs...>(),0));
    }


   public:
    template<class P>
      static
      decltype(auto)
      apply(const P& p)
      {
        return apply_imp(p, TypeCo<Fields...>(),0);
      }
  };

public:
  typedef  typename bundle<fs...>::codomain_type domain_type;
  struct result:field<result, decltype(Nop_op<domain_type>::apply(domain_type()))> {};
  typedef   typename result::codomain_type codomain_type;

 template<class P>
 static
 codomain_type
 apply(const P& p)
 {
     return codomain_type({Nop_op<domain_type>::apply(p)});
 }
};


template<class Op, typename Field>
struct result:field<result<Op,Field>,decltype(Op::apply(typename Field::field_value_type()))>
{
  constexpr const char* ClassName(){return Op::ClassName()+"("+Field::ClassName()+")";};

};


template<class Op,class f>
struct Op_imp
{
private:
  template<typename T>
  struct op_op;
  template<typename... Fields>
  struct op_op<ParametersT<Fields...>>
  {
    typedef ParametersT<result<Op,Fields>...> codomain_type;

  template<class P>
      static
      ParametersT<result<Op,Fields>...>
      apply(const P& p)
      {
        return {{Op::apply(get<Fields>(p))...}};
      }
  };

public:
  typedef  typename f::codomain_type domain_type;
  typedef   typename op_op<domain_type>::codomain_type codomain_type;

 template<class P>
 static
 codomain_type
 apply(const P& p)
 {
     return op_op<domain_type>::apply(p);
 }
};

template<class... fs>
using sum=compose<NOp_imp<sum_imp,fs...>,bundle<fs...>>;

template<class... fs>
using prod=compose<NOp_imp<prod_imp,fs...>,bundle<fs...>>;

template<class... fs>
using minus=compose<Op_imp<minus_imp,fs...>,bundle<fs...>>;




#endif // FUNCTIONT_H

};
