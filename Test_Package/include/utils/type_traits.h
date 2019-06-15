#ifndef __utils_types_h
#define __utils_types_h

#include <tr1/type_traits>
//#include <type_traits>

namespace utils
{
/*
   typedef std::tr1::integral_constant<bool, true>      __true_type;
   typedef std::tr1::integral_constant<bool, false>     __false_type;

   template <typename T1, typename T2>
   struct is_same
   {
      typedef typename std::tr1::is_same<T1,T2>::type type;
      //enum { value = typename std::tr1::is_same<T1,T2>::value };
      enum { value = std::tr1::is_same<T1,T2>::value };
   }; */

   template <bool, typename T = void>
   struct enable_if {};

   template <typename T>
   struct enable_if<true,T> { typedef T type; };

   typedef std::tr1::integral_constant<bool, true>	true_type;
   typedef std::tr1::integral_constant<bool, false>	false_type;

   using std::tr1::is_same;
   using std::tr1::is_integral;
   using std::tr1::is_floating_point;
   //using std::tr1::enable_if;
}

#endif
