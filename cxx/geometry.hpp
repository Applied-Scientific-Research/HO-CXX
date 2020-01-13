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
#include <array>
#include <memory>

#include "wtimer.hpp"
#include "memory.hpp"
#include "array.hpp"
#include "basis.hpp"

namespace HighOrderFEM
{

enum class BndryType : int { Default = 0, Dirichlet = 1, Neumann };

constexpr std::array<int,4> nghbrFace { 1, 0, 3, 2 }; // reflect across the shared face.
constexpr std::array<int,4> tensor2fem{ 3, 1, 0, 2 }; // switch from indexed (WESN) to FEM CCW notation with S=0

struct BaseGeometryType
{
   typedef double value_type;
   typedef size_t index_type;

   virtual int getK(void) const = 0;
   virtual int getL(void) const = 0;

   virtual index_type getNumElements(void) const = 0;
   virtual index_type getNumBoundaryElements(void) const = 0;

   virtual void import ( const size_t Nel, const size_t NelB,
                       double* Vol_Jac_in,
                       double* Vol_Dx_iDxsi_j_in,
                       double* Face_Jac_in,
                       double* Face_Acoef_in,
                       double* Face_Bcoef_in,
                       double* Face_Norm_in,
                       int*    elemID_in,
                       int*    bndryElementID_in,
                       int*    BC_Switch_in) = 0;
};

template <int _K, int _L>
struct GeometryType : BaseGeometryType
{
   enum : int { K = _K,
                L = _L };

   index_type Nel, NelB;

   int getK(void) const { return K; }
   int getL(void) const { return L; }

   index_type getNumElements(void) const { return this->Nel; }
   index_type getNumBoundaryElements(void) const { return this->NelB; }

   typedef StaticArrayType< value_type[K] > K_ArrayType;
   typedef StaticArrayType< value_type[K][K] > KxK_ArrayType;
   typedef StaticArrayType< value_type[2][2][K][K] > KxKx2x2_ArrayType;
   typedef StaticArrayType< value_type[4][K] > Kx4_ArrayType;
   typedef StaticArrayType< int[4] > NghbrListType;

   DynamicArrayType< BndryType > bndryType;

   DynamicArrayType< NghbrListType > elemNghborID; ///< List of element neighbors for each element in FEM order.
   DynamicArrayType< int > bndryElementID;      ///< Volume element ID's for each boundary element.

   DynamicArrayType< KxK_ArrayType > Vol_Jac;
   DynamicArrayType< KxKx2x2_ArrayType > Vol_Dx_iDxsi_j;
   DynamicArrayType< Kx4_ArrayType > Face_Jac;
   DynamicArrayType< Kx4_ArrayType > Face_Acoef;
   DynamicArrayType< Kx4_ArrayType > Face_Bcoef;
   DynamicArrayType< Kx4_ArrayType > Face_Norm;

   GeometryType (void)
   {
      printf("GeometryType<%d,%d>::GeometryType(void)\n", K, L);
   }

   ~GeometryType ()
   {
      printf("GeometryType<%d,%d>::~GeometryType()\n", K, L);
   }

   template <int K, int L>
   struct element_metrics
   {
      typedef double value_type;

      const size_t el;
      const GeometryType<K,L>& parent;

      element_metrics(const GeometryType<K,L>& g, const size_t el) : parent(g), el(el)
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
      return element_metrics<K,L>( *this, el );
   }

   inline size_t getNeighborID(const size_t el, const int f, const bool fem_order = true) const
   {
      if (fem_order)
         return this->elemNghborID(el)[ tensor2fem[f] ];
      else
         return this->elemNghborID(el)[ f ];
   }

   inline BndryType getBoundaryType(const int bel) const
   {
      return this->bndryType(bel);
   }

   //! Convert a neighbor element index to a boundary index.
   inline int getBoundaryID( const size_t nghbr_id ) const
   {
      return nghbr_id - this->Nel;
   }

   //! Test if an element index represents a boundary.
   inline bool isBoundaryElement( const size_t elem_id ) const
   {
      return elem_id < this->Nel;
   }

   void import ( const size_t Nel, const size_t NelB,
                       double* Vol_Jac_in,
                       double* Vol_Dx_iDxsi_j_in,
                       double* Face_Jac_in,
                       double* Face_Acoef_in,
                       double* Face_Bcoef_in,
                       double* Face_Norm_in,
                       int*    elemID_in,
                       int*    bndryElementID_in,
                       int*    BC_Switch_in)
   {
      printf("Inside GeometryType<%d,%d>::import(...)\n", K, L);

      this->Nel  = Nel;
      this->NelB = NelB;

      this->Vol_Jac       .setData( (KxK_ArrayType*)Vol_Jac_in, Nel );
      this->Vol_Dx_iDxsi_j.setData( (KxKx2x2_ArrayType*)Vol_Dx_iDxsi_j_in, Nel );
      this->Face_Jac      .setData( (Kx4_ArrayType*)Face_Jac_in, Nel );
      this->Face_Acoef    .setData( (Kx4_ArrayType*)Face_Acoef_in, Nel );
      this->Face_Bcoef    .setData( (Kx4_ArrayType*)Face_Bcoef_in, Nel );
      this->Face_Norm     .setData( (Kx4_ArrayType*)Face_Norm_in, Nel );

      this->elemNghborID  .resize(Nel);
      this->bndryElementID.resize(NelB);
      this->bndryType     .resize(NelB);

      for (int el(0); el < Nel; ++el)
         for (int f(0); f < 4; ++f)
         {
            const int nghbor = elemID_in[ f + 4*el ];
            if ( nghbor < 0 )
            {
               // Boundary element id.
               const int bel = (-nghbor) - 1;
               elemNghborID(el)[f] = bel + Nel; // append the boundary id's *after* the elements so both lists can be 0-based (more easily).
               bndryType(bel) = ( BC_Switch_in[bel] == 1 ) ? BndryType::Dirichlet :
                                ( BC_Switch_in[bel] == 2 ) ? BndryType::Neumann   :
                                                             BndryType::Default;
            }
            else
            {
               // Normal element-to-element connectivity in FEM format.
               elemNghborID(el)[f] = nghbor - 1;
            }
         }

      for (int bel(0); bel < NelB; ++bel)
         bndryElementID(bel) = bndryElementID_in[bel] - 1;
   }
};

struct GeometryFactoryType
{
   std::shared_ptr< BaseGeometryType > allocate( const int K, const int L );

   std::shared_ptr< BaseGeometryType > get( const int K, const int L );

   static GeometryFactoryType* Instance(void);

private:
   std::map< int, std::shared_ptr< BaseGeometryType > > ptrs;

   GeometryFactoryType();
   ~GeometryFactoryType();
   GeometryFactoryType& operator=(const GeometryFactoryType& f);

   static GeometryFactoryType* m_instance;
};

extern GeometryFactoryType GeometryFactory;

} // end namespace HighOrderFEM

#endif
