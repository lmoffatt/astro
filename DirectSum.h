#ifndef DIRECTSUM_H
#define DIRECTSUM_H
#include <type_traits>

template <class...>
class direct_sum_impl;

template <class...>
class direct_sums_impl;

template <class... C>
using direct_sums=typename direct_sums_impl<C...>::type;




template <class C, class R>
using direct_sum=typename direct_sum_impl<C,R>::type;


template <template<class...> class TC,
          class... I0,
          class... Cs>
class direct_sum_impl<TC<I0...>,TC<>,TC<>,TC<Cs...>>
{
   typedef    TC<Cs...>   type;

};


template <template<class...> class TC,
          class... I0,
          class J,
          class... Js,
          class... Cs>
class direct_sum_impl<TC<I0...>,TC<>,TC<J,Js...>,TC<Cs...>>
{
   typedef
   typename direct_sum_impl<TC<>,TC<I0...>,TC<Js...>,TC<Cs...,J>>::type
   type;

};


template <template<class...> class TC,
          class... I0,
          class I,
          class... Is,
          class J,
          class... Js,
          class... Cs>
class direct_sum_impl<TC<I0...>,TC<I,Is...>,TC<J,Js...>,TC<Cs...>>
{
   typedef
   std::conditional_t<std::is_same<I,J>::value,
   typename direct_sum_impl<TC<I0...>,TC<Is...>,TC<Js...>,TC<Cs...>>::type,
   typename direct_sum_impl<TC<I0...,I>,TC<Is...>,TC<J,Js...>,TC<Cs...>>::type>
   type;

};




template <template<class...> class TC,class... Is, class... Js>
class direct_sum_impl<TC<Is...>,TC<Js...>>
{
  typedef
typename
direct_sum_impl<TC<Is...>,TC<>,TC<Js...>,TC<Is...>>
::type
type;
};


template <class A>
struct direct_sums_impl<A>
{
  typedef
A
type;

};



template <class A, class B, class...C>
struct direct_sums_impl<A,B,C...>
{
  typedef
typename
direct_sums_impl<direct_sum<A,B>,C...>
::type
type;

};



#endif // DIRECTSUM_H
