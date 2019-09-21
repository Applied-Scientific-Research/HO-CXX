#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <type_traits>
#include <algorithm>
#include <cmath>

#include "wtimer.hpp"
#include "memory.hpp"
#include "array.hpp"

#include "aplles_interface.h"
#include "splib/csr_matrix.h"
#include "splib/multiply.h"
#include "splib/vector.h"
#include "splib/spblas.h"

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

// Test AMGCL

#ifndef SOLVER_BACKEND_BUILTIN
#  define SOLVER_BACKEND_BUILTIN
#endif

#include <tuple>
#include <vector>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/backend/interface.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

template <typename T>
struct PointerWrapper
{
   T *m_data;
   size_t len;

   typedef T value_type;
   typedef T type;

   PointerWrapper( T* ptr, size_t len ) : m_data(ptr), len(len) {}

   const T* end(void) const { return m_data + len; }
   const T* begin(void) const { return m_data; }

   const T operator[](const size_t& i) const { return m_data[i]; }
};

template <class CSRMatrix, class NodalVector>
void eval_amgcl ( const CSRMatrix& A, NodalVector& x_ref, const NodalVector& b, const int Knod )
{
   typedef amgcl::backend::builtin<double> Backend;

   const size_t nrows = A.num_rows, nnz = A.num_values;

   Backend::params backend_prm;
   //backend_prm.block_size = 4;

   auto rowptr = &A.rowptr[0];
   auto colidx = &A.colidx[0];
   auto values = &A.values[0];

   amgcl::backend::crs< double, int > crs( nrows, nrows, PointerWrapper<int   >(rowptr,nrows+1),
                                                         PointerWrapper<int   >(colidx,nnz),
                                                         PointerWrapper<double>(values,nnz) );

   {
      typedef amgcl::amg<
           Backend,
           amgcl::coarsening::smoothed_aggregation,
           //amgcl::relaxation::spai0
           //amgcl::relaxation::gauss_seidel
           amgcl::relaxation::damped_jacobi
           > AMG;

      AMG::params amg_prm;
      amg_prm.coarse_enough = 500;
    //amg_prm.coarsening.aggr.eps_strong = 0;
    //amg_prm.coarsening.aggr.block_size = 9;
      amg_prm.npre = 1; amg_prm.npost = 2;
      amg_prm.relax.damping = 0.53333333;

      AMG P( crs, amg_prm, backend_prm );
      std::cout << "AMG: " << std::endl;
      std::cout << P << std::endl;

      typedef amgcl::solver::fgmres<AMG::backend_type> SolverType;

      SolverType::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;
      solver_prms.maxiter = 500;

      SolverType solver( nrows, solver_prms );

      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;

      std::vector< double > f( nrows );
      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];
      for (int i = 0; i < nrows; ++i)
         f[i] = bptr[i];

      int niters = 0;
      double resid = 0;
      std::tie( niters, resid ) = solver( P, f, x );

      std::vector< double > r( nrows );
      amgcl::backend::residual( f, crs, x, r);
      auto norm_r = std::sqrt( amgcl::backend::inner_product( r, r ) );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                << "Norm(r):    " << norm_r << std::endl
                << std::endl;
   }

/*
   {
      boost::property_tree::ptree prm;

      typedef amgcl::make_solver<
                     amgcl::runtime::preconditioner<Backend>,
                     amgcl::runtime::solver::wrapper<Backend>
                     > Solver;

      prm.put("solver.type", "fgmres");
      prm.put("solver.tol", 0);
      prm.put("solver.abstol", 1e-10);
      prm.put("solver.maxiter", 100);
      prm.put("solver.M", 16);
      prm.put("precond.class", "amg");
      //prm.put("precond.coarsening.type", "smoothed_aggregation");
      prm.put("precond.coarsening.type", "aggregation");
      prm.put("precond.coarse_enough", "500");
      prm.put("precond.relax.type", "spai0");

      Solver solver( crs, prm, backend_prm );
      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;
   }
*/

/*
   {
      typedef amgcl::relaxation::ilu0< Backend > ilu0Type;
      typedef amgcl::relaxation::as_preconditioner<
         Backend,
         amgcl::relaxation::ilu0
         > ILU;

      ILU P( crs );
      std::cout << "ILU: " << std::endl;
      std::cout << P << std::endl;

      typedef amgcl::solver::fgmres< ILU::backend_type > SolverType;

      SolverType::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;

      SolverType solver( nrows, solver_prms );

      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;

      std::vector< double > f( nrows );
      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];
      for (int i = 0; i < nrows; ++i)
         f[i] = bptr[i];

      int niters = 0;
      double resid = 0;
      std::tie( niters, resid ) = solver( P, f, x );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                << std::endl;
   }
*/
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

   {
      typedef splib::csr_matrix<double> csr_matrix;
      csr_matrix *csr = reinterpret_cast< csr_matrix* >( A_handle );
      printf("CSR matrix: %d %d %d\n", csr->num_rows, csr->num_columns, csr->num_values);
   }

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

   {
      typedef splib::csr_matrix<double> csr_matrix;
      csr_matrix *csr = reinterpret_cast< csr_matrix* >( LaplacianData.MatrixHandle );
      printf("CSR matrix: %d %d %d\n", csr->num_rows, csr->num_columns, csr->num_values);

      auto xptr = x.getRawPointer();
      auto bptr = b.getRawPointer();

      splib::Vector<double> _x( xptr, csr->num_rows );
      splib::Vector<double> _b( bptr, csr->num_rows );
      splib::Vector<double> _r( csr->num_rows );

      splib::gemv( 1.0, *csr, _x, -1.0, _b, _r );
      auto norm = splib::spblas::norm( _r );
      auto norm_b = splib::spblas::norm( _b );
      printf("norm: %e %e\n", norm, norm_b);

      eval_amgcl ( *csr, x, b, Knod );
   }

   return;
}

}
