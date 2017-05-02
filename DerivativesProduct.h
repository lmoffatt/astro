#ifndef DERIVATIVESPRODUCT_H
#define DERIVATIVESPRODUCT_H
#include "Derivatives.h"
#include "ParametersT.h"

/// gradiente por Jacoviano
/// dlogL/dmi * dmi/dbetai
///

template<class Gr,class Ja, typename L,typename B>
typename Diff<L,B>::field_value_type
get_Diff_LB_(const Gr& G, const Ja& J, TypeCo<L>,TypeCo<B>,TypeCo<>, typename Diff<L,B>::field_value_type sum)
{
    return std::move(sum);
};


template<class Gr,class Ja, typename L,typename B, typename M,typename... Ms>
typename Diff<L,B>::field_value_type
get_Diff_LB_(const Gr& G, const Ja& J, TypeCo<L>,TypeCo<B>,TypeCo<M,Ms...>, typename Diff<L,B>::field_value_type sum)
{
    sum+=G.template get<Diff<L,M>>()*J.template get<Diff<M,B>>();
    return P_Mult(G,J, TypeCo<L>(),TypeCo<B>(),TypeCo<Ms...>(), std::move(sum));
};

template<class Gr,class Ja, typename L,typename B, typename M,typename... Ms>
typename Diff<L,B>::field_value_type
get_Diff_LB(const Gr& G, const Ja& J, TypeCo<L>,TypeCo<B>,TypeCo<M,Ms...>)
{
  typename Diff<L,B>::field_value_type sum=G.template get<Diff<L,M>>()*J.template get<Diff<M,B>>();
  return get_Diff_LB_(G,J, TypeCo<L>(),TypeCo<B>(),TypeCo<Ms...>(), std::move(sum));
};


template<class Gr,class Ja, class D,typename L,typename... Ms>
D P_Mult_(const Gr& G, const Ja& J, TypeCo<L>,TypeCo<>,TypeCo<Ms...>,D out)
{
      return std::move(out);

}


template<class Gr,class Ja, class D,typename L,typename B, typename...Bs,typename... Ms>
D P_Mult_(const Gr& G, const Ja& J, TypeCo<L>,TypeCo<B, Bs...>,TypeCo<Ms...>,
        D out)
{

   out.template get<Diff<L,B>>()=get_Diff_LB(G,J,TypeCo<L>(),TypeCo<B>(),TypeCo<Ms...>());
   return P_Mult_(G,J, TypeCo<L>(),TypeCo<Bs...>(),TypeCo<Ms...>(),std::move(out));

}



template<typename L,typename...Ms, typename B,typename...Bs>
Diff<L,ParametersT<Bs...>>
P_Mult(const Diff<L,ParametersT<Ms...>>& G, const Diff<ParametersT<Ms...>,ParametersT<B,Bs...>>& J)
{
   Diff<L,ParametersT<Bs...>> out;
   out.template get<Diff<L,B>>()=get_Diff_LB(G,J,TypeCo<L>(),TypeCo<B>(),TypeCo<Ms...>());
   return P_Mult_(G, J, TypeCo<L>(),TypeCo<Bs...>(),TypeCo<Ms...>(),std::move(out));

}


template<class X0,class X1, class Xout,typename I,typename J, typename K,typename...Ks>
  Xout get_Diff_IJ_(const X0& A, const X1& B, TypeCo<I>,TypeCo<J>,TypeCo<K,Ks...>)
{
  typename Diff<I,J>::field_value_type sum=A.template get<Diff<I,K>>()*B.template get<Diff<K,J>>();
  return get_Diff_IJ_(A,B, TypeCo<I>(),TypeCo<J>(),TypeCo<Ks...>(), std::move(sum));
};



template<class X0,class X1, class Xout,typename I,typename J, typename K,typename...Ks>
Xout get_Diff_IJ(const X0& A, const X1& B, TypeCo<I>,TypeCo<J>,TypeCo<K,Ks...>)
{
  typename Diff<I,J>::field_value_type sum=A.template get<Diff<I,K>>()*B.template get<Diff<K,J>>();
  return get_Diff_IJ_(A,B, TypeCo<I>(),TypeCo<J>(),TypeCo<Ks...>(), std::move(sum));
};









template<class X0,class X1, class Xout,typename I,typename... Is,typename J,typename...Js, typename...Ks>
Xout P_Mult__(const X0& A, const X1& B, TypeCo<I,Is...>,TypeCo<Js...>,TypeCo<Ks...>,
        Xout out)
{
  out.template get<Diff<I,J>>()=get_Diff_IJ();
}


template<typename... Is,typename...Js, typename...Ks>
Diff<ParametersT<Is...>,ParametersT<Js...>>
P_Mult(const Diff<ParametersT<Is...>,ParametersT<Ks...>>& A, const Diff<ParametersT<Ks...>,ParametersT<Js...>>& B)
{
   Diff<ParametersT<Is...>,ParametersT<Js...>> out;
   return P_Mult__(A, B, TypeCo<Is...>(),TypeCo<Js...>(),TypeCo<Ks...>(),std::move(out));

}



template<typename L,typename...Ms, typename B,typename...Bs>
Diff<L,ParametersT<Bs...>>
P_QuadraticForm(const Diff<ParametersT<Ms...>,ParametersT<B,Bs...>>& J, const Diff2<L,ParametersT<Ms...>,ParametersT<Ms...>>& H )
{
   Diff2<L,ParametersT<Bs...>,ParametersT<Bs...>> out;
   return P_Mult_(J, H,TypeCo<L>(),TypeCo<Bs...>(),TypeCo<Ms...>(),std::move(out));

}













#endif // DERIVATIVESPRODUCT_H
