#ifndef DERIVATIVES_H
#define DERIVATIVES_H


#include "myTuples.h"
#include "Matrix.h"
#include "ParametersT.h"


// esta clase guarda las derivadas entre fields
template<typename F, typename X>
struct der_;

// esta clase indica las derivadas entre los valores de losfields
template<typename F, typename X>
struct derv_;

// esta clase indica las derivadas entre parametros
template<typename F, typename X>
struct Diff_;


// detalles de implementacion de la derivada de un parametro con otro parametro
template<typename...Ts>
struct Diff_impl_;

// detalles de implementacion de la derivada de un field con otro parametro
template<typename...Ts>
struct diff_impl_;



// reparte las derivadas segun sea param o fields
template<typename F, typename X>
struct D_impl;


//impletenta D<param, param>
template<typename F, typename X>
struct Diff_impl;

// implementa D<F,param>
template<typename F, typename X>
struct diff_impl;

// implementa D<F,X>
template<typename F, typename X>
struct der_impl
{
  typedef der_<F,X>  type;
};


template<typename F, typename X>
using Diff=typename D_impl<F,X>::type;



template<typename F, typename X, typename Y>
class Derivative2_;

template<typename...Ts>
class Derivative2_impl_;

template<typename...Ts>
class Derivative2_impl;

template<typename F, typename X, typename Y>
using Diff2=typename Derivative2_impl<F,X,Y>::type;


template<typename...Ts>
using Diff2_=typename Derivative2_impl_<Ts...>::type;


template<class...>
struct isTypeConteiner: std::false_type{};

template<template <class...> class TC, typename...Ts>
struct isTypeConteiner<TC<Ts...>>: std::true_type{};


template< class, class = void >
struct isField: std::false_type{};

template<typename T>
struct isField
    <T, void_t<typename T::field_value_type>> : std::true_type  {

};










template<typename Field_f, typename Field_x>
struct der_
{
  typedef typename derv_<typename Field_f::field_value_type,typename Field_x::field_value_type>::field_value_type field_value_type;

  constexpr const char* ClassName(){return "d_"+Field_f::ClassName()+"__d_"+Field_x::ClassName();};
};




template<template<typename...> class TC,
         typename F,
         typename X,
         typename ...dFdXs>
struct diff_impl_<F,TC<X>, TC<dFdXs...>>
{
  typedef   TC<dFdXs...,Diff<F,X>> type;
};



template<template<typename...> class TC,
         typename F,
         typename X,
         typename... Xs,
         typename ...dFdXs>
struct diff_impl_<F,TC<X,Xs...>, TC<dFdXs...>>
{
  typedef
  typename
  diff_impl_<
  F,
  TC<Xs...>,
  TC<dFdXs...,Diff<F,X>>
  >::type type;
};


template<template<typename...> class TC,
         typename F, typename X,
         typename... Xs>
struct diff_impl<F,TC<X,Xs...>>
{
  typedef
  typename
  diff_impl_<
  F,
  TC<Xs...>,
  TC<Diff<F,X>>
  >::type type;
};


template<template<typename...> class TC,
         typename... dFdXs_s,
         typename... dFdXs>
struct Diff_impl_<TC<TC<dFdXs_s...>,dFdXs...>>
{
  typedef  TC<dFdXs...,dFdXs_s...> type;

};

template<template<typename...> class TC,
         typename F,
         typename... Xs,
         typename... dFdXs_s,
         typename... dFdXs>
struct Diff_impl_<TC<F>,TC<Xs...>, TC<TC<dFdXs_s...>,dFdXs...>>
{
  typedef
  typename
  Diff_impl_<
  TC<Diff<F,TC<Xs...>>,dFdXs...,dFdXs_s...>
  >::type type;

};

template<template<typename...> class TC,
         typename F,
         typename... Fs,
         typename... Xs,
         typename... dFdXs_s,
         typename... dFdXs>
struct Diff_impl_<TC<F,Fs...>,
    TC<Xs...>,
    TC<TC<dFdXs_s...>,dFdXs...>>
{
  typedef
  typename
  Diff_impl_<
  TC<Fs...>,
  TC<Xs...>,
  TC<Diff<F,TC<Xs...>>,dFdXs...,dFdXs_s...>
  >::type type;

};




template<template<typename...> class TC,
         typename F,
         typename... Fs,
         typename... Xs>
struct Diff_impl<TC<F,Fs...>,TC<Xs...>>
{
  typedef
  typename
  Diff_impl_<
  TC<Fs...>,
  TC<Xs...>,
  TC<Diff<F,TC<Xs...>>>
  >::type type;

};



template<>
struct derv_<double,double>
{
  typedef double field_value_type;
};

template<>
struct derv_<double,std::vector<double>>
{
  typedef std::vector<double> field_value_type;
};

template<>
struct derv_<std::vector<double>,double>
{
  typedef std::vector<double> field_value_type;
};

template<>
struct derv_<double,M_Matrix<double>>
{
  typedef M_Matrix<double> field_value_type;
};

template<>
struct derv_<M_Matrix<double>,double>
{
  typedef M_Matrix<double> field_value_type;
};

template<>
struct derv_<M_Matrix<double>,M_Matrix<double>>
{
  typedef M_Matrix<M_Matrix<double>> field_value_type;
};

template<>
struct derv_<std::vector<double>,std::vector<double>>
{
  typedef M_Matrix<double> field_value_type;
};

template<typename F, typename X>
struct Diff_impl
{
  typedef derv_<F,X> type;
};

template <class F, class X>
struct D_impl

{

  template <class F_, class X_>
  static
  std::enable_if_t<isField<X_>::value, typename der_impl<F_,X_>::type> D(TypeCo<F_>,TypeCo<X_>);

  template <class F_, class X_>
  static
  std::enable_if_t<!isField<X_>::value&&isField<F_>::value, typename diff_impl<F_,X_>::type> D(TypeCo<F_>,TypeCo<X_>);

  template <class F_, class X_>
  static
  std::enable_if_t<!isField<X_>::value&&!isField<F_>::value, typename Diff_impl<F_,X_>::type> D(TypeCo<F_>,TypeCo<X_>);


public:
  typedef decltype(D<F,X>(TypeCo<F>(),TypeCo<X>())) type;



};


template <class F, class... Xs>
struct D_impl<F,ParametersT<Xs...>>
{
    typedef Diff_impl<F,ParametersT<Xs...>> type;
};


//impletenta D<param, param>
template <class F, class... Xs>
struct Diff_impl<F,ParametersT<Xs...>>: public diff_impl<F,ParametersT<Xs...>>::type{};

//impletenta D<param, param>
template <class... Fs, class... Xs>
struct Diff_impl<ParametersT<Fs...>,ParametersT<Xs...>>: public Diff_impl_<ParametersT<Fs...>,ParametersT<Xs...>>::type{};


//---------------------

template<typename Field_f, typename Field_x, typename Field_y>
struct Derivative2_
{
  typedef typename Diff2<
  typename Field_f::field_value_type,
  typename Field_x::field_value_type,
  typename Field_y::field_value_type
  >::field_value_type field_value_type;

  constexpr const char* ClassName(){return "d2_"+Field_f::ClassName()+"__d_"+Field_x::ClassName()+"__d_"+Field_y::ClassName();};
};

template<typename Field_f, typename Field_x>
struct Derivative2_<Field_f,Field_x,Field_x>
{
  typedef typename Diff2<
  typename Field_f::field_value_type,
  typename Field_x::field_value_type,
  typename Field_x::field_value_type
  >::field_value_type field_value_type;

  constexpr const char* ClassName(){return "d2_"+Field_f::ClassName()+"__d_"+Field_x::ClassName()+"_2";};
};


template<typename F, typename X, typename Y>
struct Derivative2_impl<F,X,Y
    >
{
  typedef Derivative2_<F,X,Y> type;
};


template<typename...> class Derivative2_impl_2;

template<template<typename...> class TC,
         typename F, typename X,
         typename Y,
         typename...dFdXdYs>
struct Derivative2_impl_2<F,X,TC<Y>,TC<dFdXdYs...>>
{
  typedef TC<dFdXdYs...,Diff2<F,X,Y>> type;
};



template<template<typename...> class TC,
         typename F, typename X,
         typename Y,
         typename... Ys,
         typename...dFdXdYs>
struct Derivative2_impl_2<F,X,TC<Y,Ys...>,TC<dFdXdYs...>>
{
  typedef
  typename
  Derivative2_impl_2<
  F,
  X,
  TC<Ys...>,
  TC<dFdXdYs...,Diff2<F,X,Y>>
  >::type type;
};


template<template<typename...> class TC,
         typename F,
         typename X,
         typename Y,
         typename... Ys>
struct Derivative2_impl_<F,X,
    TC<Y,Ys...>>
{
  typedef
  typename
  Derivative2_impl_2<
  F,
  X,
  TC<Ys...>,
  TC<Diff2<F,X,Y>>
  >::type type;
};




template<typename...> class Derivative2_impl_1;
template<typename...> class Derivative2_impl_11;

template<template<typename...> class TC,
         typename...dFsdXdYs_s,
         typename...dFdXdYs>
struct Derivative2_impl_11<
    TC<TC<dFsdXdYs_s...>,
    dFdXdYs...> >
{
  typedef
  TC<dFdXdYs...,dFsdXdYs_s...>  type;

};



template<template<typename...> class TC,
         typename F,
         typename X,
         typename...Ys,
         typename...dFdXdYs_s,
         typename...dFdXdYs>
struct Derivative2_impl_1<F,TC<X>, TC<Ys...>,
    TC<TC<dFdXdYs_s...>,dFdXdYs...>>
{
  typedef
  typename
  Derivative2_impl_11<
  TC<Diff2_<F,X,TC<Ys...>>,dFdXdYs...,dFdXdYs_s...>
  >::type type;

};



template<template<typename...> class TC,
         typename F,
         typename X,
         typename... Xs,
         typename...Ys,
         typename...dFdXdYs_s,
         typename...dFdXdYs>
struct Derivative2_impl_1<F,TC<X,Xs...>, TC<Ys...>,
    TC<TC<dFdXdYs_s...>,dFdXdYs...>>
{
  typedef
  typename
  Derivative2_impl_1<
  F,
  TC<Xs...>,
  TC<Ys...>,
  TC<Diff2_<F,X,TC<Ys...>>,dFdXdYs...,dFdXdYs_s...>
  >::type type;

};


template<template<typename...> class TC,
         typename F,
         typename X,
         typename... Xs,
         typename...Ys>
struct Derivative2_impl<
    F,
    TC<X,Xs...>,
    TC<Ys...>>
{
  typedef
  typename
  Derivative2_impl_1<
  F,
  TC<Xs...>,
  TC<Ys...>,
  TC<Diff2_<F,X,TC<Ys...>>>
  >::type type;

};

template<typename...>
struct Derivative2_impl_0;
template<typename...>
struct Derivative2_impl_00;


template<template<typename...> class TC,
         typename...dFs_dXdYs,
         typename...dFdXdYs>
struct Derivative2_impl_00<
    TC<TC<dFs_dXdYs...>,
    dFdXdYs...> >
{
  typedef
  TC<dFdXdYs...,dFs_dXdYs...>  type;

};


template<template<typename...> class TC,
         typename F,
         typename... Xs,
         typename...Ys,
         typename...dFs_dXdYs,
         typename...dFdXdYs>
struct Derivative2_impl_0<TC<F>,TC<Xs...>, TC<Ys...>,
    TC<TC<dFs_dXdYs...>,dFdXdYs...> >
{
  typedef
  typename
  Derivative2_impl_0<
  TC<Diff2<F,TC<Xs...>,TC<Ys...>>,dFdXdYs...,dFs_dXdYs...>
  >::type type;

};


template<template<typename...> class TC,
         typename F,
         typename... Fs,
         typename... Xs,
         typename...Ys,
         typename...dFs_dXdYs,
         typename...dFdXdYs
         >
struct Derivative2_impl_0<TC<F,Fs...>,TC<Xs...>, TC<Ys...>,
    TC<TC<dFs_dXdYs...>,dFdXdYs...> >
{
  typedef
  typename
  Derivative2_impl_0<
  TC<Fs...>,
  TC<Xs...>,
  TC<Ys...>,
  TC<Diff2<F,TC<Xs...>,TC<Ys...>>,dFdXdYs...,dFs_dXdYs...>
  >::type type;

};


template<template<typename...> class TC,
         typename F,
         typename... Fs,
         typename... Xs,
         typename...Ys>
struct Derivative2_impl<TC<F,Fs...>,TC<Xs...>, TC<Ys...>>
{
  typedef
  typename
  Derivative2_impl_0<
  TC<Fs...>,
  TC<Xs...>,
  TC<Ys...>,
  TC<Diff2<F,TC<Xs...>, TC<Ys...>>>
  >::type
  type;

};



template<>
struct Derivative2_<double,double,double>
{
  typedef double field_value_type;
};

template<>
struct Derivative2_<double,std::vector<double>,double>
{
  typedef std::vector<double> field_value_type;
};

template<>
struct Derivative2_<double,double,std::vector<double>>
{
  typedef std::vector<double> field_value_type;
};


template<>
struct Derivative2_<std::vector<double>,double,double>
{
  typedef std::vector<double> field_value_type;
};
template<>
struct Derivative2_<std::vector<double>,std::vector<double>,double>
{
  typedef std::vector<std::vector<double>> field_value_type;
};

template<>
struct Derivative2_<std::vector<double>,double,std::vector<double>>
{
  typedef std::vector<std::vector<double>> field_value_type;
};

template<>
struct Derivative2_<std::vector<double>,std::vector<double>,std::vector<double>>
{
  typedef std::vector<std::vector<std::vector<double>>> field_value_type;
};


template<>
struct Derivative2_<double,M_Matrix<double>,double>
{
  typedef M_Matrix<double> field_value_type;
};

template<>
struct Derivative2_<double,double,M_Matrix<double>>
{
  typedef M_Matrix<double> field_value_type;
};

template<>
struct Derivative2_<double,M_Matrix<double>,M_Matrix<double>>
{
  typedef M_Matrix<double> field_value_type;
};

template<>
struct Derivative2_<M_Matrix<double>,double,double>
{
  typedef M_Matrix<double> field_value_type;
};






#endif // DERIVATIVES_H
