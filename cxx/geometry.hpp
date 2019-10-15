#ifndef __geometry_hpp
#define __geometry_hpp

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

namespace HighOrderFEM
{

template <int _K, int _L, typename T = double>
struct GeometryType
{
   typedef T value_type;

   enum : int { K = _K,
                      L = _L };

   typedef StaticArrayType< value_type[K] > K_ArrayType;
   typedef StaticArrayType< value_type[K][K] > KxK_ArrayType;
   typedef StaticArrayType< value_type[2][2][K][K] > KxKx2x2_ArrayType;
   typedef StaticArrayType< value_type[4][K] > Kx4_ArrayType;

   DynamicArrayType< KxK_ArrayType > Vol_Jac;
   DynamicArrayType< KxKx2x2_ArrayType > Vol_Dx_iDxsi_j;
   DynamicArrayType< Kx4_ArrayType > Face_Jac;
   DynamicArrayType< Kx4_ArrayType > Face_Acoef;
   DynamicArrayType< Kx4_ArrayType > Face_Bcoef;
   DynamicArrayType< Kx4_ArrayType > Face_Norm;

   explicit GeometryType ( const size_t Nel,
                       double* Vol_Jac_in,
                       double* Vol_Dx_iDxsi_j_in,
                       double* Face_Jac_in,
                       double* Face_Acoef_in,
                       double* Face_Bcoef_in,
                       double* Face_Norm_in )
      : 
         Vol_Jac( (KxK_ArrayType*)Vol_Jac_in, Nel ),
         Vol_Dx_iDxsi_j( (KxKx2x2_ArrayType*)Vol_Dx_iDxsi_j_in, Nel ),
         Face_Jac( (Kx4_ArrayType*)Face_Jac_in, Nel ),
         Face_Acoef( (Kx4_ArrayType*)Face_Acoef_in, Nel ),
         Face_Bcoef( (Kx4_ArrayType*)Face_Bcoef_in, Nel ),
         Face_Norm( (Kx4_ArrayType*)Face_Norm_in, Nel )
   {
      std::cout << "explicit GeometryType::GeometryType(...)" << std::endl;
   }

   GeometryType (void)
   {
      std::cout << "GeometryType::GeometryType(void)" << std::endl;
   }

   ~GeometryType ()
   {
      std::cout << "GeometryType::~GeometryType()" << std::endl;
   }

   template <int K, int L, typename value_type>
   struct element_metrics
   {
      const size_t el;
      const GeometryType<K,L,value_type>& parent;

      element_metrics(const GeometryType<K,L,value_type>& g, const size_t el) : parent(g), el(el)
      {}

      const value_type& xxsi( const int i, const int j ) const { return parent.Vol_Dx_iDxsi_j[this->el](i,j,1-1,1-1); };
      const value_type& yxsi( const int i, const int j ) const { return parent.Vol_Dx_iDxsi_j[this->el](i,j,2-1,1-1); };
      const value_type& xeta( const int i, const int j ) const { return parent.Vol_Dx_iDxsi_j[this->el](i,j,1-1,2-1); };
      const value_type& yeta( const int i, const int j ) const { return parent.Vol_Dx_iDxsi_j[this->el](i,j,2-1,2-1); };
      const value_type& jac ( const int i, const int j ) const { return parent.Vol_Jac[this->el](i,j); };

      const value_type faceA ( const int i, const int f, const int e) const { return parent.Face_Acoef[e](i,f) / parent.Face_Jac[e](i,f); };
      const value_type faceB ( const int i, const int f, const int e) const { return parent.Face_Bcoef[e](i,f) / parent.Face_Jac[e](i,f); };
      const value_type faceN ( const int i, const int f, const int e) const { return parent.Face_Norm[e](i,f); };

      const value_type faceA ( const int i, const int f) const { return this->faceA(i,f,this->el); }
      const value_type faceB ( const int i, const int f) const { return this->faceB(i,f,this->el); }
      const value_type faceN ( const int i, const int f) const { return this->faceN(i,f,this->el); }
   };

   const auto getElementMetrics(const size_t el) const
   {
      return element_metrics<K,L,value_type>( *this, el );
   }
};

} // end namespace HighOrderFEM

#endif
