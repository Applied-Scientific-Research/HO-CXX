#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <type_traits>
#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>

#include "wtimer.hpp"
#include "memory.hpp"
#include "array.hpp"
#include "basis.hpp"
#include "geometry.hpp"

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include "macros.h"


namespace HighOrderFEM
{

namespace // private
{
   inline int make_hash( const int K, const int L ) 
   {
      return __make_unique_from_K_and_L(K,L);
   }

}

GeometryFactoryType::
   GeometryFactoryType() : ptrs()
   {
      printf("Inside GeometryFactoryType%x\n", this);
   }

GeometryFactoryType::
   ~GeometryFactoryType()
   {
      printf("Inside ~GeometryFactoryType: %x\n", this);
      for ( auto obj : ptrs )
      {
         auto sptr = obj.second;
         printf("obj: %d %x\n", obj.first, sptr.get(), sptr->getK(), sptr->getL());
      }
   }

   std::shared_ptr< BaseGeometryType >
GeometryFactoryType::
   allocate( const int K, const int L )
   {
      printf("Inside GeometryFactoryType::allocate( %d %d )\n", K, L);

#define _op(_K,_L) { \
       auto sptr = std::make_shared< GeometryType<_K,_L> >(); \
       this->ptrs[hash] = std::dynamic_pointer_cast< BaseGeometryType >(sptr); \
       break; \
   }

      const auto hash = make_hash(K,L);
      auto out = this->ptrs.find( hash );
      if ( out == this->ptrs.end() )
      {
         __case_for_each_K_and_L(K,L)

         printf("Created new BaseGeometryType %d %d %x\n", K,L,this->ptrs[hash].get());

         return this->ptrs[hash];
      }
      else
      {
         printf("Found existing BaseGeometryType %d %d %x\n", K,L,out->second.get());
         return out->second;
      }

#undef _op
   }

   std::shared_ptr< BaseGeometryType >
GeometryFactoryType::
   get( const int K, const int L )
   {
      printf("Inside GeometryFactoryType::get( %d %d )\n", K, L);

      const auto hash = make_hash(K,L);
      auto out = this->ptrs.find( hash );
      if ( out == this->ptrs.end() )
      {
         fprintf(stderr,"Requested BaseGeometryType object (%d,%d) does not exist\n", K,L);
         exit(1);
      }

      return out->second;
   }

GeometryFactoryType* GeometryFactoryType::m_instance = NULL;

   GeometryFactoryType*
GeometryFactoryType::
   Instance()
   {
      if ( m_instance == NULL)
      {
         m_instance = new GeometryFactoryType;
      }

      return m_instance;
   }

} // namespace HO
