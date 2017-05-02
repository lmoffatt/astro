#ifndef CARTESIAN_H
#define CARTESIAN_H

#include "myTuples.h"
#include "Matrix.h"

template<typename...Ts>
class cartesian_;

template<typename...Ts>
class cartesian_impl;

template<typename...Ts>
using cartesian=typename cartesian_impl<Ts...>::type;

template<typename FieldA, typename FieldB>
struct cartesian_<FieldA,FieldB>
{
  typedef typename cartesian<typename FieldA::field_value_type,typename FieldB::field_value_type>::field_value_type field_value_type;
};

template<template<typename...> class TC,
         typename A,  typename ...AxBs>
struct cartesian_impl<A,TC<>, TC<AxBs...>>
{
    typedef TC<AxBs...> type;
};



template<template<typename...> class TC,
         typename A, typename B,
         typename... Bs, typename ...AxBs>
struct cartesian_impl<A,TC<B,Bs...>, TC<AxBs...>>
{
    typedef
  cartesian<
  A,
  TC<Bs...>,
  TC<AxBs...,cartesian<A,B>>
  > type;
};


template<template<typename...> class TC,
         typename A, typename B,
         typename... Bs>
struct cartesian_impl<A,TC<B,Bs...>>
{
  typedef
cartesian<
A,
TC<Bs...>,
TC<cartesian<A,B>>
> type;
};


template<template<typename...> class TC,
         typename... Bs,
         typename... AxBs_s>
struct cartesian_impl<TC<>,TC<Bs...>, TC<AxBs_s...>>
{
  typedef  TC<AxBs_s...> type;

};


template<template<typename...> class TC,
         typename A,
         typename... As,
         typename... Bs,
         typename... AxBs_s>
struct cartesian_impl<TC<A,As...>,TC<Bs...>, TC<AxBs_s...>>
{
  typedef
cartesian<
TC<As...>,
TC<Bs...>,
TC<AxBs_s...,cartesian<A,TC<Bs...>>>
> type;

};




template<template<typename...> class TC,
         typename A,
         typename... As, typename... Bs>
struct cartesian_impl<TC<A,As...>,TC<Bs...>>
{
  typedef
cartesian<
TC<As...>,
TC<Bs...>,
TC<cartesian<A,TC<Bs...>>>
> type;

};



template<>
struct cartesian_<double,double>
{
   typedef double field_value_type;
};

template<typename T>
struct cartesian_<double,T>
{
   typedef T field_value_type;
};

template<typename T>
struct cartesian_<T,double>
{
   typedef T field_value_type;
};

template<>
struct cartesian_<M_Matrix<double>,M_Matrix<double>>
{
   typedef M_Matrix<double> field_value_type;
};

template<>
struct cartesian_<std::vector<double>,std::vector<double>>
{
   typedef M_Matrix<double> field_value_type;
};









template <typename... Ts>
class commonTypes_imple;

template <class C, class R>
using commonTypes=commonTypes_imple<C,R>;


template <template<typename...>class C,typename... Rs, typename R>
class commonTypes_imple<C<>,R,C<Rs...>>
{
   typedef C<Rs...,R> type;
};

template <template<typename...>class C,typename... Ts,typename... Rs, typename R>
class commonTypes_imple<C<R,Ts...>,R,C<Rs...>>
{
   typedef C<Ts...,Rs...> type;
};

template <template<typename...>class C,typename T,typename... Ts,typename... Rs, typename R>
class commonTypes_imple<C<T,Ts...>,R,C<Rs...>>
{
   typedef typename commonTypes_imple<C<Ts...>,R,C<Rs...,T>>::type type;
};


template <template<typename...>class C,typename... Ts,typename R>
class commonTypes_imple<C<Ts...>,R>
{
   typedef typename commonTypes_imple<C<Ts...>,R,C<>>::type type;
};


template <template<typename...>class C,typename... Ts>
class commonTypes_imple<C<Ts...>,C<>>
{
   typedef C<Ts...> type;
};




template <template<typename...>class C,typename... Ts,typename R,typename... Rs>
class commonTypes_imple<C<Ts...>,C<R,Rs...>>
{
   typedef typename commonTypes_imple<
  typename commonTypes_imple<C<Ts...>,R>::type,
  C<Rs...>
  >::type type;
};





template<template<typename...> class C, typename...Ts>
decltype(auto)
commonTypes_impl(C<Ts...>,C<>,top)
{
  return C<Ts...>();
}




template<template<typename...> class C, typename...Ts, typename... Rs, typename R,  typename=decltype(std::get<R>(std::tuple<Ts...,R>()))>
decltype(auto)
commonTypes_impl(C<Ts...>,C<R,Rs...>,top)
{
  return commonTypes_impl(C<Ts...,R>(),C<Rs...>(), top());
}

template<template<typename...> class C, typename...Ts, typename... Rs, typename R>
decltype(auto)
commonTypes_impl(C<Ts...>,C<R,Rs...>,bottom)
{
  return commonTypes_impl(C<Ts...>(),C<Rs...>(),top());
}







#endif // CARTESIAN_H
