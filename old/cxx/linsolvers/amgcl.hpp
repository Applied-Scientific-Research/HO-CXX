#ifndef __amgcl_hpp
#define __amgcl_hpp

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
#include <amgcl/solver/gmres.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/block_matrix.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// AMGCL interfaces.

namespace HighOrderFEM
{

namespace LinearSolver
{

namespace details
{

template <typename T>
struct ArrayPointerWrapper
{
   T *m_data;
   size_t len;

   typedef T value_type;
   typedef T type;

   ArrayPointerWrapper( T* ptr, size_t len ) : m_data(ptr), len(len) {}

   const T* end  (void) const { return m_data + len; }
   const T* begin(void) const { return m_data; }

   const T& operator[](const size_t& i) const { return m_data[i]; }
         T& operator[](const size_t& i)       { return m_data[i]; }
};

struct IdentityPrecon
{
   template <typename Vec1, typename Vec2>
   void apply( const Vec1& f, Vec2& x ) const
   {
      const int n = f.size();

      for (int i = 0; i < n; ++i)
         x[i] = f[i];
   }
};

template <class Solver, class Matrix>
struct SolverAsPrecon
{
   typedef Solver solver_type;
   typedef Matrix matrix_type;

   solver_type &S;
   matrix_type &A;

   SolverAsPrecon ( solver_type& solver, matrix_type& matrix ) : S(solver), A(matrix)
   {
      std::cout << "SolverAsPrecon" << std::endl
                << this->S << std::endl
                << std::endl;
   }

   template <typename Vec1, typename Vec2>
   void apply( const Vec1& f, Vec2& x ) const
   {
      const int n = f.size();

      int niters = 0;
      double resid = 0;
      std::tie( niters, resid ) = this->S( this->A, IdentityPrecon(), f, x );
      //std::cout << "SolverAsPrecon::apply" << std::endl
      //          << "Iters " << niters << std::endl
      //          << "Resid " << resid << std::endl;
   }
};

template <class Precon>
struct MixedPrecisionPrecon
{
   typedef Precon precon_type;
   typedef amgcl::backend::numa_vector<float> vector_type;

   const precon_type &P;
   vector_type f, x;

   MixedPrecisionPrecon ( const precon_type& precon ) : P(precon)
   {
      size_t nrows = amgcl::backend::rows(P.system_matrix());
      this->f.resize(nrows);
      this->x.resize(nrows);

      //std::cout << "MixedPrecisionPrecon(const precon&)" << std::endl
      //          << this->P << std::endl
      //          << nrows << std::endl
      //          << std::endl;
   }

   MixedPrecisionPrecon ( MixedPrecisionPrecon&& other ) : P(other.P), f( std::move(other.f) ), x( std::move(other.x) )
   {
      //std::cout << "MixedPrecisionPrecon(&&)" << std::endl
      //          << this->P << std::endl
      //          << this->x.size() << std::endl
      //          << std::endl;
   }

   ~MixedPrecisionPrecon()
   {
      //std::cout << "~MixedPrecisionPrecon" << std::endl;
   }

   friend std::ostream& operator<<(std::ostream &os, const MixedPrecisionPrecon &p)
   {
      os << "MixedPrecisionPrecon of ... " << std::endl;
      os << p.P << std::endl;

      return os;
   }

   template <typename Vec1, typename Vec2>
   void apply( const Vec1& f, Vec2& x ) const
   {
      const size_t n = f.size();

      auto &fout = const_cast<vector_type&>(this->f);
      auto &xout = const_cast<vector_type&>(this->x);

      int niters = 0;
      double resid = 0;

      #pragma omp parallel for
      for (int i = 0; i < n; ++i) {
         fout[i] = f[i];
         xout[i] = x[i];
      }

      this->P.apply( fout, xout );

      //this->P.apply( f, x );
      //std::cout << "MixedPrecisionPrecon::apply" << std::endl
      //          << "typeid(f): " << typeid(f).name() << std::endl
      //          << "typeid(x): " << typeid(x).name() << std::endl;

      #pragma omp parallel for
      for (int i = 0; i < n; ++i)
         x[i] = xout[i];
   }
};

template < typename U, typename V >
typename U::value_type norm2( const U& u, const V& v )
{
   const int n = u.size();
   typename U::value_type s(0);
   for (int i = 0; i < n; ++i)
      s += u[i] * v[i];

   return std::sqrt(s);
}

template <typename Solver, typename Vector1, typename Vector2, typename Precon>
void exec_solver ( Solver& S, const Vector1& f, Vector2& x, Precon* P)
{
   std::cout << "Solver: " << std::endl;
   std::cout << S << std::endl;

   if (P) {
      std::cout << "Precon: " << std::endl;
      std::cout << *P << std::endl;
   }

   auto t_start = getTimeStamp();

   int niters = 0;
   double resid = 0;

//   if (P)
//      std::tie( niters, resid ) = S( *P, f, x );
//   else
      std::tie( niters, resid ) = S( f, x );

   auto t_stop = getTimeStamp();

   //std::vector< double > r( nrows );
   const int nrows = x.size();
   Vector2 r( nrows );
   amgcl::backend::residual( f, S.system_matrix(), x, r);
   //auto norm_r = std::sqrt( amgcl::backend::inner_product( r, r ) );
   auto norm_r = norm2( r, r );

   std::cout << "Iterations: " << niters << std::endl
             << "Error:      " << resid << std::endl
             << "Norm(r):    " << norm_r << std::endl
             << "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
             << std::endl;
}

template <typename Solver, typename Vector1, typename Vector2>
void exec_solver ( Solver& S, const Vector1& f, Vector2& x)
{
   int *P = NULL;
   exec_solver(S, f, x, P);
}

boost::property_tree::ptree load_amgcl_params (void)
{
   boost::property_tree::ptree prm;

   std::string config_filename = "config.amgcl";

   char *env = getenv("AMGCL_CONFIG_FILE");
   if ( env )
      config_filename = std::string(env);

   std::cout << "config " << config_filename << std::endl;

   std::ifstream config_file( config_filename );
   if ( config_file.fail() )
   {
      fprintf(stderr,"Failed to open config.amgcl ... using default options\n");

      prm.put("solver.type", "fgmres");
      prm.put("solver.tol", 0);
      prm.put("solver.abstol", 1e-10);
      prm.put("solver.maxiter", 100);
      prm.put("solver.M", 16);
      prm.put("precond.class", "amg");
      prm.put("precond.coarsening.type", "smoothed_aggregation");
    //prm.put("precond.coarsening.type", "aggregation");
      prm.put("precond.coarse_enough", 500);
    //prm.put("precond.relax.type", "spai0");
      prm.put("precond.relax.type", "damped_jacobi");
      prm.put("precond.relax.damping", 0.53333333);
      prm.put("precond.npre", 1);
      prm.put("precond.npost", 2);
   }
   else
   {
      for( std::string line; std::getline( config_file, line ); )
      {
         //std::cout << "line " << line << std::endl;

         // Remove whitespaces.
         auto end_pos = std::remove( line.begin(), line.end(), ' ' );
         line.erase( end_pos, line.end() );

         // Ignore everything after comment symbol (#)
         auto comment_pos = line.find( '#' );
         if ( comment_pos != std::string::npos )
            line.erase( line.begin() + comment_pos, line.end() );
         //auto comment_pos = std::remove( line.begin(), line.end(), '#' );
         //line.erase( comment_pos, line.end() );

         // Find the (key, value) seperator.
         auto sep = line.find( ',' );
         if ( sep == std::string::npos )
            continue;

         std::string key = std::string( line.begin(), line.begin()+sep );
         std::string val = std::string( line.begin()+sep+1, line.end() );

         //std::cout << "key " << key << std::endl;
         //std::cout << "val " << val << std::endl;
         prm.put( key, val );
      }
   }

   boost::property_tree::write_json( std::cout, prm );

   return prm;
}

template <typename ColType, typename ValType>
struct RowSorter
{
   size_t n;
   std::vector<ColType> m_col;
   std::vector<ValType> m_val;

   template <typename T>
   bool operator() ( const T& i, const T& j ) const
   {
      return this->m_col[i] < this->m_col[j];
   }

   RowSorter( ColType *col, ValType *val, const size_t n ) : m_col(col,col+n), m_val(val,val+n), n(n)
   {
      std::vector<int> index(n);

      for (int i(0); i < n; ++i)
         index[i] = i;

      std::sort( index.begin(), index.end(), *this );

      for (int i = 0; i < n; ++i)
      {
         col[i] = this->m_col[ index[i] ];
         val[i] = this->m_val[ index[i] ];
      }
   }
};

template <typename R, typename C, typename V>
void sort_crs_rows( const int nrows, const R *rowptr, C *colidx, V *values )
{
   for (int i = 0; i < nrows; ++i)
   {
      const int j = rowptr[i];
      const int nnz = rowptr[i+1] - j;

      //printf("%lu, %lu\n", i, nnz);
      //for (int k = j; k < (j+nnz); ++k)
      //{
      //   std::cout << colidx[k];
      //   if ((k-j) % 9 == 8 or k == (j+nnz-1))
      //      std::cout << std::endl;
      //   else
      //      std::cout << ",";
      //}

      RowSorter<C,V> sorter( colidx + j, values + j, nnz );

      //printf("sorted\n");
      //for (int k = j; k < (j+nnz); ++k)
      //{
      //   std::cout << colidx[k];
      //   if ((k-j) % 9 == 8 or k == (j+nnz-1))
      //      std::cout << std::endl;
      //   else
      //      std::cout << ",";
      //}
   }
}

} // namespace details

/* Hold the preconditioner and iterative solver for the AMGCL library.
 */
struct AMGCLSolver : BaseLinearSolver
{
   typedef amgcl::backend::builtin<value_type> BackendType;		///< Primary back-end
   typedef amgcl::backend::builtin<float > mixed_BackendType;	///< Secondary back-end for the lower precision operations.

   typedef amgcl::runtime::preconditioner< mixed_BackendType > PreconditionerType;
   typedef amgcl::runtime::solver::wrapper< BackendType >      IterativeSolverType;
   typedef details::MixedPrecisionPrecon< PreconditionerType > MixedPreconditionerType;
   
   typedef amgcl::make_solver< PreconditionerType, IterativeSolverType > MakeSolverType;

   typedef typename MakeSolverType::build_matrix build_matrix;

   std::shared_ptr< MakeSolverType > shared_solver;
   std::shared_ptr< MixedPreconditionerType > shared_mixed_precon;

   build_matrix system_matrix; //!< A native copy of the system matrix since the internal
                               //!< storage will be lower precision.

   decltype( details::load_amgcl_params() ) prm;

   value_type resid;
   int num_iters;

   /// Constructor accepts any valid input matrix format.
   AMGCLSolver (void) : BaseLinearSolver( LinearSolver::LinearSolverTag::AMGCL ),
                        resid(0), num_iters(0)
   {
      std::cout << "AMGCLSolver: " << this << std::endl;

      this->prm = details::load_amgcl_params();
   }

   /// Constructor accepts any valid input matrix format.
   int build( const index_type num_rows, const std::vector<index_type>& rowptr,
                                         const std::vector<index_type>& colidx,
                                         const std::vector<value_type>& values )
   {
      std::cout << "AMGCLSolver::build: " << typeid(build_matrix).name() << std::endl;

      this->num_rows = num_rows;

      auto A = std::tie( num_rows, rowptr,
                                   colidx,
                                   values );

      this->system_matrix = A;

      this->prm.put("solver.tol", this->reltol);
      this->prm.put("solver.abstol", this->abstol);
      this->prm.put("solver.maxiter", this->maxiters);

      auto t_build_start = getTimeStamp();

      this->shared_solver = std::make_shared< MakeSolverType >( this->system_matrix, this->prm );

      auto &P = this->shared_solver->precond();
      this->shared_mixed_precon = std::make_shared< MixedPreconditionerType > ( P );

      if ( this->verbosity )
      {
         auto t_build_stop = getTimeStamp();

         std::cout << "AMGCL Build time: " << getElapsedTime( t_build_start, t_build_stop ) << std::endl;

         std::cout << "AMGCL Solver: " << std::endl << shared_solver->solver() << std::endl;
         std::cout << "AMGCL Precon: " << std::endl << *shared_mixed_precon << std::endl;
      }
   }

   ~AMGCLSolver()
   {
      std::cout << "~AMGCLSolver: " << this << std::endl;
   }

   template <typename VectorX, typename VectorF>
   int apply( VectorX& x, const VectorF& f )
   {
      auto t_start = getTimeStamp();

      auto& _A = this->shared_solver->system_matrix();
      auto& A = this->system_matrix;
      auto& P = *(this->shared_mixed_precon.get());

      auto& solver = this->shared_solver->solver();

    //std::cout << "typeid(solver) " << typeid(solver).name() << std::endl;
    //std::cout << "typeid(_A) " << typeid(_A).name() << std::endl;
    //std::cout << "typeid(A) " << typeid(A).name() << std::endl;
    //std::cout << "typeid(P) " << typeid(P).name() << std::endl;
    //std::cout << "typeid(x) " << typeid(x).name() << std::endl;
    //std::cout << "typeid(f) " << typeid(f).name() << std::endl;

      auto ret = solver.operator()( A, P, f, x );

      auto t_stop = getTimeStamp();

      auto niters = std::get<0>(ret);
      auto resid  = std::get<1>(ret);

      this->num_iters = niters;
      this->resid = resid;

      if (this->verbosity)
      {
         size_t nrows = amgcl::backend::rows(A);
         std::vector<value_type> r(nrows);
         amgcl::backend::residual( f, A, x, r );
         auto norm_r = details::norm2( r, r );

         std::cout << "AMGCL Iterations: " << niters << std::endl
                   << "      Error:      " << resid << std::endl
                   << "      Norm(r):    " << norm_r << std::endl
                   << "      Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                   << std::endl;
      }

      //int iret = ( this->maxiters > 0 ) ? ( niters < this->maxiters ) : 1;
      //return iret;
      //return ( this->resid < this->reltol ) ?
      return ( this->num_iters < this->maxiters ) ?
                   LinearSolver::SolverStatusFlags::Success :
                   LinearSolver::SolverStatusFlags::Failure;
   }

   double getResidual(void) { return this->resid; }
   int getNumIterations(void) { return this->num_iters; }

   //template <typename X, typename F>
   //int operator()( X& x, const F& f )
   //{
   //   return this->apply( x, f );
   //}

   int solve( const std::vector<value_type>& f, std::vector<value_type>& u )
   {
      return this->apply( u, f );
   }
};

} // namespace LinearSolver
} // namespace HighOrder

#endif
