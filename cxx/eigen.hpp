#ifndef __eigen_hpp
#define __eigen_hpp

#include <type_traits>
#include <memory>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "linsolver.hpp"

namespace HighOrderFEM
{
namespace LinearSolver
{

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

   //print_info ( *A );

   return A;
}

} // details

struct EigenSolver : BaseLinearSolver
{
   typedef Eigen::SparseMatrix<value_type> matrix_type;

   typedef Eigen::IncompleteLUT<value_type> precon_type;
   typedef Eigen::BiCGSTAB< matrix_type, precon_type > solver_type;

   solver_type solver;
   matrix_type A;

   double fill_factor, droptol;

   EigenSolver()
      : fill_factor(10), droptol(1e-5),
        solver(),
        BaseLinearSolver( LinearSolverTag::EIGEN )
   {
      std::cout << "EigenSolver: " << this << std::endl;
   }

   ~EigenSolver()
   {
      std::cout << "~EigenSolver: " << this << std::endl;
   }

   void setFillFactor( const double fill_factor ) { this->fill_factor = fill_factor; }
   void setDropTol( const double droptol ) { this->droptol = droptol; }

   int build( const index_type nrows, const std::vector<index_type>& rowptr,
                                      const std::vector<index_type>& colidx,
                                      const std::vector<value_type>& values )
   {
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      this->num_rows = nrows;

      this->solver.preconditioner().setFillfactor( this->fill_factor );
      this->solver.preconditioner().setDroptol( this->droptol );

      auto A_map = details::build_eigen_matrix( nrows, const_cast<index_type*>(rowptr.data()),
                                                       const_cast<index_type*>(colidx.data()),
                                                       const_cast<value_type*>(values.data()) );

      this->A = *A_map;

      auto t_start_compute = getTimeStamp();
      this->solver.compute( this->A ); // Compute the ILUT factorization and other solver/precond initializations.
      auto t_stop_compute = getTimeStamp();

      std::cout << "Eigen Build info: " << (solver.preconditioner().info() == Eigen::Success) << std::endl;
      std::cout << "Eigen Build time: " << getElapsedTime( t_start_compute, t_stop_compute ) << std::endl;

      //auto LU_nnz = solver.preconditioner().m_lu.nonZeros();
      //std::cout << "Eigen fill factor: " << double(LU_nnz) / this->A.nonZeros() << std::endl;

      return (solver.preconditioner().info() == Eigen::Success);
   }

   int solve ( const std::vector<value_type>& fv, std::vector<value_type>& out )
   {
      using namespace Eigen;
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      //std::cout << "Inside EigenSolver::solve" << std::endl;

      Map<const VectorXd> f_map(fv.data(), this->num_rows);
      VectorXd f = f_map;

      // Eigen only has a relative tolerance test.
      auto tol = this->reltol;
      if ( this->abstol > 0.0 )
      {
         auto normf = f.norm();

         if ( this->reltol > 0.0 )
         {
            tol = std::min( this->reltol, this->abstol / normf );
            fprintf(stderr,"Eigen3 only supports relative tolerances. Both reltol and abstol were specified. Selecting more strigent: %e %e %e\n", this->reltol, this->abstol, tol);
         }
         else
            tol = this->abstol / normf;
      }

      std::cout << "reltol: " << this->reltol << ", abstol: " << this->abstol << ", tol: " << tol << std::endl;
      std::cout << "maxiters: " << this->maxiters << std::endl;

      solver.setTolerance( tol );
      solver.setMaxIterations( this->maxiters );

      auto t_start = getTimeStamp();
      VectorXd u = solver.solve(f);
      auto t_stop = getTimeStamp();

      bool failed = solver.info() != Eigen::Success;

      //Map<VectorXd> u_map(up, this->num_rows);
      //VectorXd u_out = u_map;
      //for (int i = 0; i < num_rows; ++i)
      //   u_out[i] = u[i];
      for (int i = 0; i < num_rows; ++i)
         out[i] = u[i];

      if (failed)
         std::cout << "Solver failed to converge!" << std::endl;

      if (this->verbosity > 1)
      {
         std::cout << "eigen solver time: " << getElapsedTime( t_start, t_stop ) << std::endl;
         std::cout << "              tol: " << solver.tolerance() << std::endl;
         std::cout << "            error: " << solver.error() << std::endl;
         std::cout << "           niters: " << solver.iterations() << std::endl;
         std::cout << "         ||Ax-b||: " << ( this->A * u - f ).norm() << std::endl;
         std::cout << "\n";
      }

      return not(failed);
   }
};

} // LinearSolver
} // HighOrder

#endif
