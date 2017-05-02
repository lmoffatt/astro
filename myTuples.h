#ifndef MYTUPLES_H
#define MYTUPLES_H
#include <tuple>


template<class...>
using void_t=void;


template<typename...Ts>
class TypeCo
{

};


struct bottom{
  int x;
};

struct top:public bottom{};












template<typename T>
constexpr std::size_t sum_Is(T x)
{
  return x;
}


template<typename I,typename...Is>
constexpr std::size_t sum_Is(I x, Is... xs)
{
  return x+sum_Is(xs...);
}



// ------------- UTILITY---------------



template<int...> struct index_tuple{};

template<int I, typename IndexTuple, typename... Types>
struct make_indexes_impl;

template<int I, int... Indexes, typename T, typename ... Types>
struct make_indexes_impl<I, index_tuple<Indexes...>, T, Types...>
{
    typedef typename make_indexes_impl<I + 1, index_tuple<Indexes..., I>, Types...>::type type;
};

template<int I, int... Indexes>
struct make_indexes_impl<I, index_tuple<Indexes...> >
{
    typedef index_tuple<Indexes...> type;
};

template<typename ... Types>
struct make_indexes : make_indexes_impl<0, index_tuple<>, Types...>
{};

// ----------- FOR EACH -----------------
template<typename Func, typename Last>
void for_each_impl(Func&& f, Last&& last)
{
    f(last);
}

template<typename Func, typename First, typename ... Rest>
void for_each_impl(Func&& f, First&& first, Rest&&...rest)
{
    f(first);
    for_each_impl( std::forward<Func>(f), rest...);
}

template<typename Func, int ... Indexes, typename ... Args>
void for_each_helper( Func&& f, index_tuple<Indexes...>, std::tuple<Args...>&& tup)
{
    for_each_impl( std::forward<Func>(f), std::forward<Args>(std::get<Indexes>(tup))...);
}

template<typename Func, typename ... Args>
void for_each( std::tuple<Args...>& tup, Func&& f)
{
   for_each_helper(std::forward<Func>(f),
                   typename make_indexes<Args...>::type(),
                   std::forward<std::tuple<Args...>>(tup) );
}

template<typename Func, typename ... Args>
void for_each( std::tuple<Args...>&& tup, Func&& f)
{
   for_each_helper(std::forward<Func>(f),
                   typename make_indexes<Args...>::type(),
                   std::forward<std::tuple<Args...>>(tup) );
}


template<typename T, typename...Ts>
constexpr std::size_t index_of_impl(std::size_t I,TypeCo<T>,TypeCo<>)
{
   return I+1;
}



template<typename T, typename...Ts>
constexpr std::size_t index_of_impl(std::size_t I,TypeCo<T>,TypeCo<T,Ts...>)
{
   return I;
}

template<typename T,typename S, typename...Ts>
constexpr std::size_t index_of_impl(std::size_t I,TypeCo<T>,TypeCo<S,Ts...>)
{
   return index_of_impl(I+1,TypeCo<T>(),TypeCo<Ts...>());
}


template<typename T,template<typename...>class C,typename...Ts>
constexpr std::size_t index_of(C<Ts...>)
{
   return index_of_impl(0,TypeCo<T>(),TypeCo<Ts...>());
}


inline std::size_t mySize(double){return 1;}

template<typename M>
std::size_t mySize(const M& x){return x.size();}


template <typename... Ts>
std::tuple<Ts...>& self_add_impl(std::tuple<Ts...>& me, const std::tuple<Ts...> , std::index_sequence<>)
{
    return me;
}


template <std::size_t I,std::size_t...Is,typename... Ts>
std::tuple<Ts...>& self_add_impl(std::tuple<Ts...>& me, const std::tuple<Ts...> other, std::index_sequence<I,Is...>)
{
    std::get<I>(me)+=std::get<I>(other);
    return self_add_impl(me, other, std::index_sequence<Is...>());
}




template <typename... Ts>
std::tuple<Ts...>& operator+=(std::tuple<Ts...>& me, const std::tuple<Ts...> other)
{
    return self_add_impl(me, other, std::index_sequence_for<Ts...>());
}






#endif // MYTUPLES_H
