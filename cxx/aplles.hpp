#ifndef __aplles_hpp
#define __aplles_hpp

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <sstream>
#include <algorithm>

#include "aplles_interface.h"
#include "linsolver.hpp"

namespace HighOrderFEM
{
namespace LinearSolver
{

namespace aplles_details
{

#include <solver/base_solver.h>
#include <splib/csr_matrix.h>
#include <splib/vector.h>
#include <splib/multiply.h>
#include <splib/spblas.h>

   int getNumIterations(APLLES_SolverHandle_t &S_handle)
   {
      solver::BaseSolver *baseSolver = reinterpret_cast<solver::BaseSolver *>( S_handle );
      return baseSolver->num_iters;
   }

   template <typename T>
   double getResidual( APLLES_MatrixHandle_t &A_handle, T* x_ptr, T* b_ptr )
   {
      using namespace splib;

      typedef csr_matrix<T> MatrixType;
      typedef Vector<T> VectorType;

      BaseMatrix *baseA = reinterpret_cast<BaseMatrix *>(A_handle);

      const auto num_rows = baseA->num_rows;

      //VectorType r( num_rows );
      std::vector< T > r( num_rows );

      const MatrixType &A = *( reinterpret_cast<MatrixType *>(baseA) );
      const VectorType x( x_ptr, num_rows );
      const VectorType b( b_ptr, num_rows );

      gemv( 1.0, A, x, 1.0, b, r );

      double normb = spblas::norm2( b );
      double normr = spblas::norm2( r );

      double relnorm = ( normb > 0.0 ) ? normr / normb : normr;

      printf("Final actual residual = %e %e %e\n", normr, normb, relnorm );

      return normr;
   }

} // aplles_details

struct ApllesSolver : BaseLinearSolver
{
   APLLES_MatrixHandle_t A;
   APLLES_SolverHandle_t S;
   APLLES_PreconHandle_t P;

   ApllesSolver(void)
      : BaseLinearSolver( LinearSolverTag::APLLES ),
        A(NULL), S(NULL), P(NULL)
   {}

   ~ApllesSolver()
   {
      std::cout << "Destroying ApllesSolver: " << this << "\n";
      if (A) APLLES_Destroy_Matrix( A );
      if (S) APLLES_Destroy_Solver( S );
      if (P) APLLES_Destroy_Precon( P );
   }

   int build( const index_type nrows, const std::vector<index_type>& rowptr,
                                      const std::vector<index_type>& colidx,
                                      const std::vector<value_type>& values )
   {
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      auto t_start = getTimeStamp();

      APLLES_Initialize();

      this->num_rows = nrows;

      APLLES_MatrixHandle_t A_view;
      APLLES_Setup_Matrix_CSR( nrows, const_cast<index_type*>( rowptr.data() ),
                                      const_cast<index_type*>( colidx.data() ),
                                      const_cast<value_type*>( values.data() ),
                                      A_view );

      std::string format = "csr";
      APLLES_Matrix_Copy_CSR_To_Other( A_view, const_cast<char*>( format.c_str() ), this->A );

      // Default is FlexGMRES with BoomerAMG.
      const char* solver = "fgmres";
      const char* precon = "amg";

      APLLES_Setup_Solver( this->A, const_cast<char*>(solver), this->S );
      APLLES_Setup_Precon( this->A, const_cast<char*>(precon), this->P );

      auto t_end = getTimeStamp();

      printf("Aplles setup time (sec): %f\n", getElapsedTime( t_start, t_end ));
   }

   int solve ( const std::vector<value_type>& b_in, std::vector<value_type>& x_out )
   {
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      auto t_start = getTimeStamp();

      value_type *b_ptr = const_cast<value_type*>( b_in.data() );
      value_type *x_ptr = x_out.data();

      int ret = APLLES_Solve( this->A, x_ptr, b_ptr, this->S, this->P );

      auto t_end = getTimeStamp();

      //if (this->verbosity)
         printf("Aplles solve time: %f %d\n", getElapsedTime( t_start, t_end ), aplles_details::getNumIterations(this->S) );

      //if (this->verbosity)
      {
         //aplles_details::getResidual( this->A, x_ptr, b_ptr );
      }

      return 1;
   }
};

} // LinearSolver
} // HO

#endif
