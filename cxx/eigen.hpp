#ifndef __eigen_hpp
#define __eigen_hpp

#include <type_traits>
#include <memory>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "my_IncompleteLUT.h"

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

   std::cout << A->rows() << std::endl;         // Number of rows
   std::cout << A->cols() << std::endl;         // Number of columns 
   std::cout << A->nonZeros() << std::endl;     // Number of non zero values   
   std::cout << A->outerSize() << std::endl;    // Number of columns (resp. rows) for a column major (resp. row major )
   std::cout << A->innerSize() << std::endl;    // Number of rows (resp. columns) for a row major (resp. column major)
   std::cout << A->isVector() << std::endl;     // Check if A is a sparse vector or a sparse matrix
   std::cout << A->isCompressed() << std::endl; // Check if A is in compressed form

   return A;
}

template <typename IndexType, typename ValueType>
void test_eigen ( const size_t nrows, IndexType* rowptr, IndexType* colidx, ValueType* values, const ValueType* fp, const ValueType* xp )
{
   using namespace Eigen;
   using namespace HighOrderFEM;

   typedef ValueType value_type;
   typedef IndexType index_type;

   typedef MappedSparseMatrix<value_type, RowMajor, index_type> mapped_matrix_type;
   typedef SparseMatrix<value_type> matrix_type;
   typedef Matrix<value_type, Dynamic, 1> vector_type;

   auto A_map = build_eigen_matrix( nrows, rowptr, colidx, values );
   Map<const VectorXd> f_map(fp,nrows);
   VectorXd f = f_map;

   std::cout << A_map->rows() << std::endl;
   std::cout << f.rows() << std::endl;
   std::cout << f.cols() << std::endl;

   matrix_type A = *A_map;
   std::cout << A.rows() << std::endl;

   BiCGSTAB< matrix_type, my_IncompleteLUT<ValueType> > solver;
   solver.preconditioner().setFillfactor( 10 );
   solver.preconditioner().setDroptol( 1e-5 );

// BiCGSTAB< matrix_type, DiagonalPreconditioner<ValueType> > solver;

   auto normb = f.norm();
   std::cout << "normb: " << normb << std::endl;
   const value_type abstol = 1e-10;
   const value_type tol = abstol / normb; // to get abs tolerance in Eigen which only looks at reltol.

   solver.setTolerance( tol );
   solver.setMaxIterations( 200 );

   auto t_start_compute = getTimeStamp();
   solver.compute( A ); // Compute the ILUT factorization and other solver/precond initializations.
   auto t_stop_compute = getTimeStamp();

   std::cout << "Build info: " << (solver.preconditioner().info() == Success) << std::endl;
   std::cout << "Build time: " << getElapsedTime( t_start_compute, t_stop_compute ) << std::endl;

   auto LU_nnz = solver.preconditioner().m_lu.nonZeros();
   std::cout << "fill factor: " << double(LU_nnz) / A.nonZeros() << std::endl;

   auto t_start = getTimeStamp();
   VectorXd x = solver.solve(f);
   auto t_stop = getTimeStamp();
   std::cout << "solver info: " << solver.info() << std::endl;
   std::cout << "solver time: " << getElapsedTime( t_start, t_stop ) << std::endl;
   std::cout << "solver tolr: " << solver.tolerance() << std::endl;
   std::cout << "solver erro: " << solver.error() << std::endl;
   std::cout << "solver iter: " << solver.iterations() << std::endl;
   //std::cout << "type: " << typeid(x).name() << std::endl;
   //std::cout << x.rows() << std::endl;
   std::cout << "resid norm: " << ( A * x - f ).norm() << std::endl;

   //{
   //   Map<const VectorXd> _x( xp, 20 );
   //   Map<const VectorXd> _f( fp, 20 );
   //   std::cout << _x << std::endl;
   //   std::cout << _f << std::endl;
   //}

   //for (int i = 0; i < 20; i++) printf("f[%d] =  %e %e\n", i, f(i), x(i));
   //for (int i = 0; i < 20; i++) printf("f[%d] =  %e %e\n", i, fp[i], xp[i]);

   return;
}

#endif
