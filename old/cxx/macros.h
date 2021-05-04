#ifndef __macros_h
#define __macros_h

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>

#define K_VALUES (1)(3)(5)
#define L_VALUES (1)(3)(4)(5)

#define __make_unique_from_K_and_L(_K,_L) ( (_L) + (_K)*100 )

#define __for_each_K_and_L__(r, _K, _L) \
    case __make_unique_from_K_and_L((_K),(_L)): { \
       _op(_K,_L) \
       }

#define __for_each_K_and_L(r, KL) \
           __for_each_K_and_L__(r, BOOST_PP_SEQ_ELEM(0, KL), BOOST_PP_SEQ_ELEM(1, KL))

#define __case_for_each_K_and_L(Kv,Lv) \
      switch( __make_unique_from_K_and_L(Kv,Lv) ) { \
        BOOST_PP_SEQ_FOR_EACH_PRODUCT( __for_each_K_and_L, (K_VALUES)(L_VALUES) ) \
        default: \
          fprintf(stderr,"Unsupported K/L run-time values %d/%d\n", Kv,Lv); \
          exit(-1); \
      }

#endif
