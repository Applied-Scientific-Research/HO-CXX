#ifndef __forall_hpp
#define __forall_hpp

#include <functional>
#include <algorithm>

namespace HighOrderFEM
{

template <class Functor>
void forall( const int I, const Functor& f )
{
   for (int i = 0; i < I; ++i)
      f(i);
}

template <class Functor>
void forall( const int I, const int J, const Functor& f )
{
   for (int j = 0; j < J; ++j)
      for (int i = 0; i < I; ++i)
         f(i,j);
}

/*
template < typename Int, class Functor, class BinaryOp >
   typename std::enable_if< std::is_integral<Int>::value and std::is_class<Functor>::value,
                            typename BinaryOp::result_type >:: type
forall_reduce( const Int I, const Functor& f, const BinaryOp& op = std::plus<double>() )
{
   typedef typename BinaryOp::result_type value_type;

   value_type v;

   for (int i = 0; i < I; ++i)
      v = op( f(i), v );

   return v;
}
*/

namespace Reduction
{

template <typename T>
struct maxval
{
   typedef T result_type;
   constexpr T identity() const { return std::numeric_limits<T>::min(); }
   T operator()( const T& a, const T& b ) const { return std::max(a,b); }
};

template <typename T>
struct minval
{
   typedef T result_type;
   constexpr T identity() const { return std::numeric_limits<T>::max(); }
   T operator()( const T& a, const T& b ) const { return std::min(a,b); }
};

template <typename T>
struct sum
{
   typedef T result_type;
   constexpr T identity() const { return T(0); }
   T operator()( const T& a, const T& b ) const { return a+b; }
};

template <typename T>
struct product
{
   typedef T result_type;
   constexpr T identity() const { return T(1); }
   T operator()( const T& a, const T& b ) const { return a*b; }
};

}

template < typename Int, class Functor, class BinaryOp >
auto forall_reduce( const Int I, const Functor& f, const BinaryOp& op ) -> decltype( f(0) )
{
   static_assert( std::is_same< decltype(f(0)), decltype(op(0,0)) >::value, "Type mismatch between function return and binary operator");

   typedef decltype( f(0) ) value_type;

   value_type v( op.identity() );

   for (int i = 0; i < I; ++i)
      v = op( f(i), v );

   return v;
}

template < typename Int, class Functor >
auto forall_reduce( const Int I, Functor& f ) -> decltype( f(0) )
{
   typedef decltype( f(0) ) value_type;

   value_type v(0);
   std::plus<value_type> op;

   for (int i = 0; i < I; ++i)
      v = op( f(i), v );

   return v;
}

template < class Functor, class BinaryOp >
   typename BinaryOp::result_type
forall_reduce( const int I, const int J, const Functor& f, const BinaryOp& op = std::plus<double>() )
{
   typedef typename BinaryOp::result_type value_type;

   value_type v;

   for (int j = 0; j < J; ++j)
      for (int i = 0; i < I; ++i)
         v = op( f(i,j), v );

   return v;
}

} // end namespace

#endif // #forall_hpp
