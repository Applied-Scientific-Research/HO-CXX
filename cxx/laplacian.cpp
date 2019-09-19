#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <type_traits>
#include <algorithm>

#include "wtimer.hpp"
#include "memory.hpp"
#include "array.hpp"

#include "aplles_interface.h"

using namespace HighOrderFEM;

template <int _Knod>
struct LaplacianDataType
{
   enum { Knod = _Knod };

   typedef StaticArrayType< double[Knod][Knod] > ElementNodalType;

   int Nel; /// Number of elements in mesh.
   DynamicArrayType< ElementNodalType > Vol_Jac; /// Volumetric Jacobian of the element.

   StaticArrayType< double[Knod] > wgt;

   int NelB; /// Number of boundary
   DynamicArrayType< int > BoundaryPointElementIDs;

   typedef StaticArrayType< double[Knod*Knod][Knod] > BoundarySourceType;
   DynamicArrayType< BoundarySourceType > BoundarySource;

   typedef StaticArrayType< double[Knod] > BoundaryValuesType;
   DynamicArrayType< BoundaryValuesType > BoundaryValues;

   /// APLLES handles
   APLLES_MatrixHandle_t MatrixHandle;
   APLLES_SolverHandle_t SolverHandle;
   APLLES_PreconHandle_t PreconHandle;
};

typedef LaplacianDataType<3> LaplacianType;
LaplacianType LaplacianData;

template <typename T>
T sqr( const T& x ) { return x*x; }

template <class A, class B>
   typename std::enable_if< A::rank == B::rank and A::rank == 1 and A::I0 == B::I0,
                                   typename getBaseValueType<A>::type >::type
dot_product ( const A& a, const B& b )
{
   typedef typename getBaseValueType<A>::type value_type;

   value_type prod(0);
   const int len = A::I0;
   for (int i = 0; i < len; ++i)
      prod += a(i) * b(i);

   return prod;
}

extern "C"
{

void setLaplacian( const int Knod, const int Nel,
                     double wgt[], double Vol_Jac[],
                     const int NelB, int BoundaryPointElementIDs1[], double BoundarySource[], double BoundaryValues[],
                     void *A_handle, void *S_handle, void *P_handle )
{
   printf("Knod: %d\n", Knod);
   printf("Nel: %d\n", Nel);
   printf("NelB: %d\n", NelB);

   if (Knod != 3) {
      fprintf(stderr,"Knod != 3 %d\n", Knod);
      exit(-1);
   }

   LaplacianData.Nel = Nel;

   for (int i(0); i < Knod; i++) {
      LaplacianData.wgt(i) = wgt[i];
      printf("wgt[%d]= %f\n", i, LaplacianData.wgt[i]);
   }

   LaplacianData.Vol_Jac.setData( (LaplacianType::ElementNodalType *) Vol_Jac, Nel );

   LaplacianData.NelB = NelB;

   LaplacianData.BoundarySource.setData( (LaplacianType::BoundarySourceType *) BoundarySource, NelB );
   LaplacianData.BoundaryValues.setData( (LaplacianType::BoundaryValuesType  *) BoundaryValues, NelB );

   // Make a copy of this since it changes in the driver code.
   //LaplacianData.BoundaryValues.resize( NelB );
   //memcpy( LaplacianData.BoundaryValues.getPointer(), BoundaryValues, sizeof(LaplacianType::BoundaryValuesType)*NelB );

   LaplacianData.BoundaryPointElementIDs.resize( NelB );
   for (int i(0); i < NelB; ++i)
      LaplacianData.BoundaryPointElementIDs(i) = BoundaryPointElementIDs1[i] - 1;

   LaplacianData.MatrixHandle = (APLLES_MatrixHandle_t) A_handle;
   LaplacianData.SolverHandle = (APLLES_SolverHandle_t) S_handle;
   LaplacianData.PreconHandle = (APLLES_PreconHandle_t) P_handle;

   return;
}

void getLaplacian( double VorticityIn[], double psiIn[], double* BoundarySourceIn, double *BoundaryValuesIn)
{
   const int Knod = LaplacianData.Knod;
   const int Nel  = LaplacianData.Nel;
   const int NelB = LaplacianData.NelB;

   const auto& Vol_Jac = LaplacianData.Vol_Jac;
   const auto& wgt = LaplacianData.wgt;

   typedef StaticArrayType< double[Knod][Knod] > NodalType;

   const DynamicArrayType< NodalType > Vorticity( (NodalType *) VorticityIn, Nel );
   const DynamicArrayType< NodalType > psi( (NodalType *) psiIn, Nel );

   DynamicArrayType< NodalType > x(Nel), b(Nel);

   auto t_start = getTimeStamp();

   /// Evaluate the RHS vorticity for each element's nodes.

   double resid = 0, vol = 0;

   #pragma omp parallel for default(shared) //reduction(+: resid, vol)
   for (int el = 0; el < Nel; ++el)
      for (int j = 0; j < Knod; ++j)
         for (int i = 0; i < Knod; ++i)
         {
            auto& vort_el = Vorticity(el);
            auto& jac_el = Vol_Jac(el);

            //resid += wgt(i) * wgt(j) * vort_el(i,j) * jac_el(i,j);
            //vol   += wgt(i) * wgt(j) *                jac_el(i,j);

            b[el](i,j) = -jac_el(i,j) * vort_el(i,j);
            x[el](i,j) = 0.0;
         }

   //LaplacianData.BoundarySource.setData( (LaplacianType::BoundarySourceType *) BoundarySourceIn, NelB );
   //LaplacianData.BoundaryValues.setData( (LaplacianType::BoundaryValuesType  *) BoundaryValuesIn, NelB );

   /// Factor in the boundary face/element.

   //double sum_bnd(0), max_bnd(0);
   for (int bnd_el_id = 0; bnd_el_id < NelB; ++bnd_el_id)
   {
      /// Volumetric element ID
      const auto el_id = LaplacianData.BoundaryPointElementIDs(bnd_el_id);

      /// References to the element's objects.
            auto& b_el = b[el_id];
      const auto& bndry_src = LaplacianData.BoundarySource[bnd_el_id];
      const auto& bndry_val = LaplacianData.BoundaryValues[bnd_el_id];

      for (int j = 0; j < Knod; ++j)
         for (int i = 0; i < Knod; ++i)
         {
            const auto ij = i + j * Knod;

            //const auto slice = bndry_src.slice<0>( ij );

            //double dotp(0);
            //for (int l = 0; l < Knod; ++l)
            //   dotp += slice(l) * bndry_val(l);
            //   //dotp += bndry_src(l,ij) * bndry_val(l);

            double dotp = dot_product( bndry_src.slice<0>( ij ), bndry_val );

            b_el(i,j) -= dotp;

            //sum_bnd += sqr( dotp );
            //max_bnd = std::max( std::fabs( dotp ), max_bnd );
         }
   }

   double sum_b(0), max_b(0);
   //for (int el = 0; el < Nel; ++el)
   //   for (int j = 0; j < Knod; ++j)
   //      for (int i = 0; i < Knod; ++i)
   //      {
   //         sum_b += sqr( b[el](i,j) );
   //         max_b = std::max( std::fabs( b[el](i,j) ), max_b );
   //      }

   //printf("Resid, vol, eps: %e %e %e %e %e %e %e\n", resid, vol, (1.0+resid)/vol, sqrt(sum_bnd), max_bnd, sqrt(sum_b), max_b);

   auto t_middle = getTimeStamp();

   int ierr = APLLES_Solve( LaplacianData.MatrixHandle,
                            x.getRawPointer(), b.getRawPointer(),
                            LaplacianData.SolverHandle, LaplacianData.PreconHandle );

   auto t_end = getTimeStamp();

   double sum_x(0), max_x(0), err(0), ref(0);
   for (int el = 0; el < Nel; ++el)
      for (int j = 0; j < Knod; ++j)
         for (int i = 0; i < Knod; ++i)
         {
            sum_x += sqr( x[el](i,j) );
            max_x = std::max( std::fabs( x[el](i,j) ), max_x );

            auto diff = x[el](i,j) - psi[el](i,j);
            err += sqr( diff );
            ref += sqr( psi[el](i,j) );
         }

   printf("cxx solved: %d %e %e %f %f %e %e %e\n", ierr, sqrt(sum_x), max_x, getElapsedTime( t_start, t_middle ), getElapsedTime( t_middle, t_end ), sqrt(err), sqrt(ref), sqrt(err/ref));

   return;
}

}
