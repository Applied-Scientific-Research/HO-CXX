#ifndef __eigen_hpp
#define __eigen_hpp

#include <type_traits>
#include <memory>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "my_IncompleteLUT.h"

namespace details
{

template <typename MatrixType>
void print_info( const MatrixType& A, const std::string name = "Undefined" )
{
   std::cout << "Eigen Matrix: " << name << std::endl;
   std::cout << "  rows()         " << A.rows()         << std::endl;         // Number of rows
   std::cout << "  cols()         " << A.cols()         << std::endl;         // Number of columns 
   std::cout << "  nonZeros()     " << A.nonZeros()     << std::endl;     // Number of non zero values   
   std::cout << "  outerSize()    " << A.outerSize()    << std::endl;    // Number of columns (resp. rows) for a column major (resp. row major )
   std::cout << "  innerSize()    " << A.innerSize()    << std::endl;    // Number of rows (resp. columns) for a row major (resp. column major)
   std::cout << "  isVector()     " << A.isVector()     << std::endl;     // Check if A is a sparse vector or a sparse matrix
   std::cout << "  isCompressed() " << A.isCompressed() << std::endl; // Check if A is in compressed form
}

template <typename IndexType, typename ValueType>
std::shared_ptr< Eigen::MappedSparseMatrix<ValueType, Eigen::RowMajor, IndexType> >
   build_eigen_matrix ( const size_t nrows, IndexType* rowptr, IndexType* colidx, ValueType* values )
{
   typedef ValueType value_type;
   typedef IndexType index_type;

   typedef Eigen::MappedSparseMatrix<value_type, Eigen::RowMajor, index_type> matrix_type;

   typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> dynamic_vector_type;

   auto ncols = nrows;
   auto nnz   = rowptr[nrows];

   /// Copy matrix from builtin backend.
   auto A = std::make_shared< matrix_type > ( nrows, ncols, nnz, rowptr, colidx, values );

   print_info ( *A );

   return A;
}

struct EigenSolver
{
   typedef double value_type;
   typedef int    index_type;

   typedef Eigen::SparseMatrix<value_type>              matrix_type;

   typedef Eigen::BiCGSTAB< matrix_type, Eigen::my_IncompleteLUT<value_type> > solver_type;

   solver_type solver;
   matrix_type A;

   double abstol, reltol;
   int maxiters;
   double fill_factor, droptol;

   template <typename IndexType, typename ValueType>
   EigenSolver( const size_t nrows, IndexType* rowptr, IndexType* colidx, ValueType* values )
      : abstol(1e-10), reltol(0), maxiters(200), fill_factor(10), droptol(1e-5),
        solver()
   {
      this->solver.preconditioner().setFillfactor( 10 );
      this->solver.preconditioner().setDroptol( 1e-5 );

      this->build( nrows, rowptr, colidx, values );
   }

   ~EigenSolver()
   {
      std::cout << "~EigenSolver" << std::endl;
   }

   template <typename IndexType, typename ValueType>
   int build( const size_t nrows, IndexType* rowptr, IndexType* colidx, ValueType* values )
   {
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      auto A_map = build_eigen_matrix( nrows, rowptr, colidx, values );

      this->A = *A_map;

      auto t_start_compute = getTimeStamp();
      this->solver.compute( this->A ); // Compute the ILUT factorization and other solver/precond initializations.
      auto t_stop_compute = getTimeStamp();

      std::cout << "Build info: " << (solver.preconditioner().info() == Eigen::Success) << std::endl;
      std::cout << "Build time: " << getElapsedTime( t_start_compute, t_stop_compute ) << std::endl;

      auto LU_nnz = solver.preconditioner().m_lu.nonZeros();
      std::cout << "fill factor: " << double(LU_nnz) / this->A.nonZeros() << std::endl;

      return (solver.preconditioner().info() == Eigen::Success);
   }

   template < typename ValueType >
   int solve ( const size_t nrows, const ValueType* fp, const ValueType* up = NULL )
   {
      using namespace Eigen;
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      std::cout << "Inside EigenSolver::solve" << std::endl;

      Map<const VectorXd> f_map(fp, nrows);
      VectorXd f = f_map;

      auto normf = f.norm();
      std::cout << "norm(f): " << normf << std::endl;

      this->reltol = this->abstol / normf; // to get abs tolerance in Eigen which only looks at reltol.

      std::cout << "reltol: " << this->reltol << ", abstol: " << this->abstol << std::endl;
      std::cout << "maxiters: " << this->maxiters << std::endl;

      solver.setTolerance( this->reltol );
      solver.setMaxIterations( this->maxiters );

      auto t_start = getTimeStamp();
      VectorXd u = solver.solve(f);
      auto t_stop = getTimeStamp();

      bool failed = solver.info() != Eigen::Success;

      if (failed)
         std::cout << "Solver failed to converge!" << std::endl;

      std::cout << "eigen solver time: " << getElapsedTime( t_start, t_stop ) << std::endl;
      std::cout << "             tolr: " << solver.tolerance() << std::endl;
      std::cout << "             erro: " << solver.error() << std::endl;
      std::cout << "             iter: " << solver.iterations() << std::endl;

      if (true)
      std::cout << " true enorm: " << ( this->A * u - f ).norm() << std::endl;
      std::cout << "\n";

      return not(failed);
   }
};

std::shared_ptr< EigenSolver > eigen_solver;

} // namespace details

template <typename IndexType, typename ValueType>
void build_eigen ( const size_t nrows, IndexType* rowptr, IndexType* colidx, ValueType* values )
{
   details::eigen_solver = std::make_shared< details::EigenSolver > ( nrows, rowptr, colidx, values );

   return;
}

template <typename ValueType>
int test_eigen ( const size_t nrows, const ValueType* fp, const ValueType* up )
{
   return details::eigen_solver->solve( nrows, fp, up );
}

#endif
