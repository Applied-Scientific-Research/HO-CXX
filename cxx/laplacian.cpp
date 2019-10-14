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
#include <amgcl/solver/gmres.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/block_matrix.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

      std::cout << "MixedPrecisionPrecon(const precon&)" << std::endl
                << this->P << std::endl
                << nrows << std::endl
                << std::endl;
   }

   MixedPrecisionPrecon ( MixedPrecisionPrecon&& other ) : P(other.P), f( std::move(other.f) ), x( std::move(other.x) )
   {
      std::cout << "MixedPrecisionPrecon(&&)" << std::endl
                << this->P << std::endl
                << this->x.size() << std::endl
                << std::endl;
   }

   ~MixedPrecisionPrecon()
   {
      std::cout << "~MixedPrecisionPrecon" << std::endl;
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

struct AMGCL_SolverType
{
   typedef amgcl::backend::builtin<double> BackendType;
   typedef amgcl::backend::builtin<float > mixed_BackendType;

   typedef amgcl::runtime::preconditioner< mixed_BackendType > PreconditionerType;
   typedef amgcl::runtime::solver::wrapper< BackendType >      IterativeSolverType;
   typedef MixedPrecisionPrecon< PreconditionerType > MixedPreconditionerType;
   
   typedef amgcl::make_solver< PreconditionerType, IterativeSolverType > MakeSolverType;

   typedef typename MakeSolverType::build_matrix build_matrix;

   std::shared_ptr< MakeSolverType > shared_solver;
   std::shared_ptr< MixedPreconditionerType > shared_mixed_precon;

   template <typename Matrix>
   AMGCL_SolverType ( const Matrix& A )
   {
      std::cout << "typeid(build_matrix): " << typeid(build_matrix).name() << std::endl;

      auto prm = load_amgcl_params();

      auto t_build_start = getTimeStamp();

      shared_solver = std::make_shared< MakeSolverType >( A, prm );

      auto &P = shared_solver->precond();
    //MixedPreconditionerType m_precon(P);
      shared_mixed_precon = std::make_shared< MixedPreconditionerType > ( P );

      auto t_build_stop = getTimeStamp();

      std::cout << "Build time: " << getElapsedTime( t_build_start, t_build_stop ) << std::endl;

      std::cout << "Solver: " << std::endl << shared_solver->solver() << std::endl;
      std::cout << "Precon: " << std::endl << *shared_mixed_precon << std::endl;
   }

   ~AMGCL_SolverType()
   {
      std::cout << "~AMGCL_SolverType" << std::endl;
   }

   template <typename VectorX, typename VectorF>
   void apply( VectorX& x, const VectorF& f )
   {
      auto t_start = getTimeStamp();

      int niters = 0;
      double resid = 0;

      auto& A = this->shared_solver->system_matrix();

      std::tie< niters, resid > = this->shared_solver->solver()( A, this->shared_mixed_precon.get(), f, x );

      auto t_stop = getTimeStamp();

      size_t nrows = amgcl::backend::rows(A);
      decltype(x) r(nrows);
      amgcl::backend::residual( f, A, x, r);
      auto norm_r = norm2( r, r );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                << "Norm(r):    " << norm_r << std::endl
                << "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                << std::endl;
   }

   template <typename X, typename F>
   void operator()( X& x, const F& f )
   {
      this->apply( x, f );
   }
};

std::shared_ptr< AMGCL_SolverType > amgcl_solver;

template <class CSRMatrix, class NodalVector>
void eval_amgcl ( const CSRMatrix& A_in, NodalVector& x_ref, const NodalVector& b, const int Knod )
{
   typedef amgcl::backend::builtin<double> Backend;

   const size_t nrows = A_in.num_rows, nnz = A_in.num_values;

   Backend::params backend_prm;
   //backend_prm.block_size = 4;

   auto rowptr = &A_in.rowptr[0];
   auto colidx = &A_in.colidx[0];
   auto values = &A_in.values[0];

#if 0
   std::vector<int> A_rowptr( rowptr, rowptr + nrows+1 );
   std::vector<int> A_colidx( colidx, colidx + nnz );
   std::vector<double> A_values( values, values + nnz );

   sort_crs_rows( nrows, A_rowptr.data(), A_colidx.data(), A_values.data() );
#else
   //amgcl::backend::crs< double, int > A( nrows, nrows, ArrayPointerWrapper<int   >(rowptr,nrows+1),
   //                                                      ArrayPointerWrapper<int   >(colidx,nnz),
   //                                                      ArrayPointerWrapper<double>(values,nnz) );

   ArrayPointerWrapper<int   > A_rowptr(rowptr,nrows+1);
   ArrayPointerWrapper<int   > A_colidx(colidx,nnz);
   ArrayPointerWrapper<double> A_values(values,nnz);
#endif

   auto A = std::tie( nrows, A_rowptr,
                             A_colidx,
                             A_values );

   //typedef amgcl::backend::builtin<float> fBackend;

   //std::vector<float> values_f( nnz );
   //for (int i = 0; i < nnz; ++i)
   //   values_f[i] = values[i];

   //amgcl::backend::crs< float, int > A_f( nrows, nrows, ArrayPointerWrapper<int  >(rowptr,nrows+1),
   //                                                       ArrayPointerWrapper<int  >(colidx,nnz),
   //                                                       values_f );

/*
   {
      typedef amgcl::amg<
           Backend,
           //fBackend,
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

      //AMG P( A_f, amg_prm, backend_prm );
      AMG P( A, amg_prm, backend_prm );
      std::cout << "AMG: " << std::endl;
      std::cout << P << std::endl;

      //typedef amgcl::solver::fgmres< AMG::backend_type > SolverType;
      typedef amgcl::solver::fgmres< Backend > SolverType;

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
      amgcl::backend::residual( f, A, x, r);
      auto norm_r = std::sqrt( amgcl::backend::inner_product( r, r ) );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                << "Norm(r):    " << norm_r << std::endl
                << std::endl;
   }
*/

   auto prm = load_amgcl_params();

   if (0)
   {
      typedef amgcl::backend::builtin<float> fBackend;

      typedef amgcl::make_solver<
                     amgcl::runtime::preconditioner<Backend>,
                   //amgcl::runtime::preconditioner<fBackend>,
                     amgcl::runtime::solver::wrapper<Backend>
                     > Solver;

      auto t_build_start = getTimeStamp();

      Solver S( A, prm, backend_prm );

      auto t_build_stop = getTimeStamp();

      std::cout << "Build time: " << getElapsedTime( t_build_start, t_build_stop ) << std::endl;

      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];
      //ArrayPointerWrapper< double > f( bptr, nrows );
      std::vector< double > f( bptr, bptr + nrows );

      //for (int j = 0; j < 3; ++j)
      //{
      //   for (int i = 0; i < nrows; ++i)
      //      x[i] = 0.0;

         exec_solver( S, f, x );
      //}
   }

   if (1)
   {
      typedef amgcl::backend::builtin<float> fBackend;

      typedef amgcl::make_solver<
                     amgcl::runtime::preconditioner<fBackend>,
                     amgcl::runtime::solver::wrapper<Backend>
                     > Solver;

      auto t_build_start = getTimeStamp();

      Solver S( A, prm, backend_prm );

      auto &P = S.precond();
      MixedPrecisionPrecon< decltype(P) > P_f(P);

      auto t_build_stop = getTimeStamp();

      std::cout << "Build time: " << getElapsedTime( t_build_start, t_build_stop ) << std::endl;

      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];
      //ArrayPointerWrapper< double > f( bptr, nrows );
      std::vector< double > f( bptr, bptr + nrows );

      //exec_solver( S, f, x, &P_f );
      std::cout << "Solver: " << std::endl << S.solver() << std::endl;
      std::cout << "Precon: " << std::endl << P_f << std::endl;

      auto t_start = getTimeStamp();

      int niters = 0;
      double resid = 0;

      std::tie( niters, resid ) = S.solver()( A, P_f, f, x );

      auto t_stop = getTimeStamp();

      decltype(x) r(nrows);
      amgcl::backend::residual( f, A, x, r);
      auto norm_r = norm2( r, r );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                << "Norm(r):    " << norm_r << std::endl
                << "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                << std::endl;
   }

   if (0)
   {
      typedef amgcl::backend::builtin<float> fBackend;

      typedef amgcl::amg<
           fBackend,
           amgcl::coarsening::smoothed_aggregation,
           amgcl::relaxation::damped_jacobi
           > Precon;

      typedef amgcl::solver::fgmres< Backend > Solver;

      Precon::params precon_prms;
      precon_prms.coarse_enough = 500;
    //precon_prms.coarsening.aggr.eps_strong = 0;
    //precon_prms.coarsening.aggr.block_size = 9;
      precon_prms.npre = 1; precon_prms.npost = 2;
      precon_prms.relax.damping = 0.53333333;

      Solver::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;
      solver_prms.maxiter = 500;

      auto t_build_start = getTimeStamp();

      auto A = std::tie( nrows, A_rowptr,
                                A_colidx,
                                A_values );

      std::vector<float> A_values_f( nnz );
      for (int i = 0; i < nnz; ++i)
         A_values_f[i] = A_values[i];

      auto A_f = std::tie( nrows, A_rowptr, A_colidx, A_values_f );

      Precon P( A_f, precon_prms, backend_prm );
      MixedPrecisionPrecon<Precon> P_f(P);

      Solver S( nrows, solver_prms, backend_prm );

      auto t_build_stop = getTimeStamp();

      std::cout << "Build time: " << getElapsedTime( t_build_start, t_build_stop ) << std::endl;
      std::cout << "Solver: " << std::endl << S << std::endl;
      std::cout << "Precon: " << std::endl << P << std::endl;

      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];
      std::vector< double > f( bptr, bptr + nrows );

      auto t_start = getTimeStamp();

      int niters = 0;
      double resid = 0;

      std::tie( niters, resid ) = S( A, P_f, f, x );

      auto t_stop = getTimeStamp();

      decltype(x) r(nrows);
      amgcl::backend::residual( f, A, x, r);
      auto norm_r = norm2( r, r );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                << "Norm(r):    " << norm_r << std::endl
                << "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                << std::endl;
   }

/*
   {
      typedef amgcl::relaxation::as_preconditioner<
         Backend,
         //amgcl::relaxation::ilu0
         amgcl::solver::gmres
         > Precon;

      Precon P( A );
      std::cout << "Precon: " << std::endl;
      std::cout << P << std::endl;

      typedef amgcl::solver::fgmres< Backend > Solver;

      Solver::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;

      Solver solver( nrows, solver_prms );

      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;

      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];

      std::vector< double > f( bptr, bptr + nrows );

      exec_solver( solver, f, x );
   }
*/

   if (0)
   {
      IdentityPrecon P;

      typedef amgcl::solver::fgmres< Backend > Solver;

      Solver::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;
      solver_prms.maxiter = 500;

      Solver solver( nrows, solver_prms );

      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;

      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];

      std::vector< double > f( bptr, bptr + nrows );

      int niters = 0;
      double resid = 0;

      std::tie( niters, resid ) = solver( A, P, f, x );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                //<< "Norm(r):    " << norm_r << std::endl
                //<< "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                << std::endl;
   }

   if (0)
   {
      IdentityPrecon P;

      const int BS = 9;

      typedef amgcl::static_matrix<double, BS, BS> value_type;
      typedef amgcl::static_matrix<double, BS,  1> rhs_type;

      typedef amgcl::backend::builtin< value_type > BBackend;

      typedef amgcl::solver::fgmres< BBackend > Solver;
      //typedef amgcl::make_solver<
      //               amgcl::runtime::preconditioner<BBackend>,
      //               amgcl::runtime::solver::wrapper<BBackend>
      //            > Solver;

      auto Ab = amgcl::adapter::block_matrix< value_type >(std::tie( nrows, A_rowptr, A_colidx, A_values ));

      Solver::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;
      solver_prms.maxiter = 500;

      Solver solver( nrows/BS, solver_prms );
      //Solver solver( Ab, prm );

      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;

      std::vector< double > x( nrows, 0 );

      rhs_type const *fptr = reinterpret_cast<rhs_type const *>( &b[0] );
      rhs_type       *xptr = reinterpret_cast<rhs_type       *>( &x[0] );

      amgcl::backend::numa_vector<rhs_type> F(fptr, fptr + nrows/BS);
      amgcl::backend::numa_vector<rhs_type> X(xptr, xptr + nrows/BS);

      int niters = 0;
      double resid = 0;

      std::tie( niters, resid ) = solver( Ab, P, F, X );
      //std::tie( niters, resid ) = solver( F, X );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                //<< "Norm(r):    " << norm_r << std::endl
                //<< "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                << std::endl;
   }

   if (0)
   {
      typedef amgcl::solver::gmres< Backend > PreconSolver;
      //typedef amgcl::solver::bicgstab< Backend > PreconSolver;
      typedef SolverAsPrecon< PreconSolver, decltype(A) > Precon;

      PreconSolver::params precon_solver_prms;
      precon_solver_prms.M = 16;
      precon_solver_prms.tol = 0.1;
      precon_solver_prms.maxiter = 32;

      PreconSolver precon_solver( nrows, precon_solver_prms );

      Precon P( precon_solver, A );

      typedef amgcl::solver::fgmres< Backend > Solver;

      Solver::params solver_prms;
      solver_prms.M = 16;
      solver_prms.tol = 0; // to get only abstol
      solver_prms.abstol = 1e-10;
      solver_prms.maxiter = 500;

      Solver solver( nrows, solver_prms );

      std::cout << "Solver: " << std::endl;
      std::cout << solver << std::endl;

      std::vector< double > x( nrows, 0 );

      double *bptr = (double *) &b[0];

      std::vector< double > f( bptr, bptr + nrows );

      int niters = 0;
      double resid = 0;

      std::tie( niters, resid ) = solver( A, P, f, x );

      std::cout << "Iterations: " << niters << std::endl
                << "Error:      " << resid << std::endl
                //<< "Norm(r):    " << norm_r << std::endl
                //<< "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                << std::endl;
   }
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

   {
      BasisFunctions<3,4> basis;
   }

   return;
}

void assembleLaplacian( const int _Knod, const int Nel, const int NelB,
                        double Vol_Jac_in[], double Vol_Dx_iDxsi_j_in[],
                        double Face_Acoef_in[], double Face_Bcoef_in[], double Face_Jac_in[], double Face_Norm_in[],
                        double lapCenter_in[], double lapNghbor_in[], double bndrySrc_in[],
                        int elemID_in[], int BC_Switch_in[]
                      )
{
   printf("Inside assembleLaplacian\n");
   printf("Knod: %d\n", _Knod);
   printf("Nel: %d\n", Nel);
   printf("NelB: %d\n", NelB);

   if (_Knod != 3) {
      fprintf(stderr,"Knod != 3 %d\n", _Knod);
      exit(-1);
   }

   const int K = 3, L = 4;
   const int Ksq = K*K;

   typedef StaticArrayType< double[K] > K_ArrayType;
   typedef StaticArrayType< double[K][K] > KxK_ArrayType;
   typedef StaticArrayType< double[2][2][K][K] > KxKx2x2_ArrayType;
   typedef StaticArrayType< double[4][K] > Kx4_ArrayType;
   typedef StaticArrayType< double[Ksq][Ksq] > KsqxKsq_ArrayType;
   typedef StaticArrayType< double[4][Ksq][Ksq] > KsqxKsqx4_ArrayType;
   typedef StaticArrayType< double[Ksq][K] > KsqxK_ArrayType;
   typedef StaticArrayType< double[Ksq][K] > KsqxK_ArrayType;

   DynamicArrayType< KxK_ArrayType > Vol_Jac( (KxK_ArrayType*)Vol_Jac_in, Nel ) ;
   DynamicArrayType< KxKx2x2_ArrayType > Vol_Dx_iDxsi_j( (KxKx2x2_ArrayType*)Vol_Dx_iDxsi_j_in, Nel ) ;
   DynamicArrayType< Kx4_ArrayType > Face_Jac( (Kx4_ArrayType*)Face_Jac_in, Nel ) ;
   DynamicArrayType< Kx4_ArrayType > Face_Acoef( (Kx4_ArrayType*)Face_Acoef_in, Nel ) ;
   DynamicArrayType< Kx4_ArrayType > Face_Bcoef( (Kx4_ArrayType*)Face_Bcoef_in, Nel ) ;
   DynamicArrayType< Kx4_ArrayType > Face_Norm( (Kx4_ArrayType*)Face_Norm_in, Nel ) ;
   DynamicArrayType< KsqxKsq_ArrayType > lapCenter_ref( (KsqxKsq_ArrayType*)lapCenter_in, Nel ) ;
   DynamicArrayType< KsqxKsqx4_ArrayType > lapNghbor_ref( (KsqxKsqx4_ArrayType*)lapNghbor_in, Nel ) ;

   DynamicArrayType< KsqxK_ArrayType > bndrySrc_ref( (KsqxK_ArrayType*)bndrySrc_in, NelB );

   DynamicArrayType< StaticArrayType< int[4] > > elemNghborID(Nel);

   enum class BndryTypes : int { Default = 0, Dirichlet = 1, Neumann };

   DynamicArrayType< StaticArrayType< double[Ksq][K] > > bndrySrc( NelB );
   DynamicArrayType< BndryTypes > bndryType(NelB);

   for (int bel(0); bel < NelB; ++bel)
      bndrySrc(bel).set(0.0);

   for (int el(0); el < Nel; ++el)
      for (int f(0); f < 4; ++f)
      {
         const int nghbor = elemID_in[ f + 4*el ];
         if ( nghbor < 0 )
         {
            // Boundary element id.
            const int bel = (-nghbor) - 1;
            elemNghborID(el)[f] = bel + Nel; // append the boundary id's *after* the elements so both lists can be 0-based (more easily).
            bndryType(bel) = ( BC_Switch_in[bel] == 1 ) ? BndryTypes::Dirichlet :
                             ( BC_Switch_in[bel] == 2 ) ? BndryTypes::Neumann   :
                                                          BndryTypes::Default;
         }
         else
         {
            // Normal element-to-element connectivity in FEM format.
            elemNghborID(el)[f] = nghbor - 1;
         }
      }

   // Build up some common basis function terms first.
   KxK_ArrayType E00, E01, E10, E11, E, C, D, F00, F01, F10, F11;

   const double gLprime[] = { -0.25*K*(K+1),
                               0.25*K*(K+1) };

   BasisFunctions<K,L> Basis;

   forall( K, K, [&]( const int i, const int j )
      {
         E00(i,j) = 0.5 * (Basis.SolnBndryLgrangeBasis(i,0) * Basis.NodesGradRadau(j,0));
         E01(i,j) = 0.5 * (Basis.SolnBndryLgrangeBasis(i,0) * Basis.NodesGradRadau(j,1));
         E10(i,j) = 0.5 * (Basis.SolnBndryLgrangeBasis(i,1) * Basis.NodesGradRadau(j,0));
         E11(i,j) = 0.5 * (Basis.SolnBndryLgrangeBasis(i,1) * Basis.NodesGradRadau(j,1));
         E  (i,j) = 2.0 * (E00(i,j) + E11(i,j));
         C  (i,j) = Basis.SolnNodesGradLgrangeBasis(i,j) - (E00(i,j) + E11(i,j));
         D  (i,j) = C(i,j) - (E00(i,j) + E11(i,j));
         F00(i,j) = 0.5 * (Basis.SolnBndryGradLgrangeBasis(i,0) - gLprime[0] * Basis.SolnBndryLgrangeBasis(i,0)) * Basis.NodesGradRadau(j,0);
         F01(i,j) = 0.5 * (Basis.SolnBndryGradLgrangeBasis(i,0) - gLprime[0] * Basis.SolnBndryLgrangeBasis(i,0)) * Basis.NodesGradRadau(j,1);
         F10(i,j) = 0.5 * (Basis.SolnBndryGradLgrangeBasis(i,1) - gLprime[1] * Basis.SolnBndryLgrangeBasis(i,1)) * Basis.NodesGradRadau(j,0);
         F11(i,j) = 0.5 * (Basis.SolnBndryGradLgrangeBasis(i,1) - gLprime[1] * Basis.SolnBndryLgrangeBasis(i,1)) * Basis.NodesGradRadau(j,1);
      });

   // Now build up the matrix coefficients on a per-element basis. This will look like a block-matrix.
   auto blk_idx = [&]( const size_t blk_row, const size_t blk_col ) -> size_t { return blk_row * size_t(Nel) + blk_col; };
   auto blk_ij  = [&]( const size_t blk_idx ) -> std::pair<size_t,size_t> {
         size_t blk_row = blk_idx / size_t(Nel);
         size_t blk_col = ( blk_idx - blk_row * size_t(Nel) );
         return std::make_pair( blk_row, blk_col );
      };

   const int nghbrFace[] = { 1, 0, 3, 2 }; // reflect across the shared face.
   const int tensor2fem[] = { 3, 1, 0, 2 }; // switch from indexed (WESN) to FEM CCW notation with S=0

   std::map<size_t, KsqxKsq_ArrayType > coefs;

   auto t_start = getTimeStamp();

   for (int el(0); el < Nel; ++el)
   {
      // 2nd-order volumetric and boundary metrics
      auto xxsi = [&]( const int i, const int j ) { return Vol_Dx_iDxsi_j[el](i,j,1-1,1-1); };
      auto yxsi = [&]( const int i, const int j ) { return Vol_Dx_iDxsi_j[el](i,j,2-1,1-1); };
      auto xeta = [&]( const int i, const int j ) { return Vol_Dx_iDxsi_j[el](i,j,1-1,2-1); };
      auto yeta = [&]( const int i, const int j ) { return Vol_Dx_iDxsi_j[el](i,j,2-1,2-1); };
      auto jac  = [&]( const int i, const int j ) { return Vol_Jac[el](i,j); };

      KxK_ArrayType Ax, Ay, Bx, By;

      forall( K, K, [&] ( const int i, const int j ) {
            Ax(i,j) = ( xeta(i,j)*xeta(i,j) + yeta(i,j)*yeta(i,j) ) / jac(i,j);
          //Ay(i,j) = ( xxsi(i,j)*xxsi(i,j) + yxsi(i,j)*yxsi(i,j) ) / jac(i,j);
            Ay(j,i) = ( xxsi(i,j)*xxsi(i,j) + yxsi(i,j)*yxsi(i,j) ) / jac(i,j);
            Bx(i,j) = ( xxsi(i,j)*xeta(i,j) + yxsi(i,j)*yeta(i,j) ) / jac(i,j);
            By(j,i) = Bx(i,j); // tranposed only
         });

      typedef std::function<double(const int, const int, const int)> face_func_type;
      face_func_type faceA = [&]( const int i, const int fid, const int elid) { return Face_Acoef[elid](i,fid) / Face_Jac[elid](i,fid); };
      face_func_type faceB = [&]( const int i, const int fid, const int elid) { return Face_Bcoef[elid](i,fid) / Face_Jac[elid](i,fid); };
      face_func_type faceN = [&]( const int i, const int fid, const int elid) { return Face_Norm[elid](i,fid); };

      //if (el < 10) {
      //   auto print_array = [&]( const KxK_ArrayType& A, const char* S ) {
      //         std::cout << S << std::endl;
      //         forall( K, K, [&]( const int i, const int j ) { printf("%d %d %e\n", i, j, A(i,j)); });
      //      };

      //   std::cout << "el metrics " << el << std::endl;
      //   print_array( Ax, "Ax");
      //   print_array( Ay, "Ay");
      //   print_array( Bx, "Bx");
      //   print_array( By, "By");

      //   auto print_face = [&]( const face_func_type& A, const char *S ) {
      //         std::cout << S << std::endl;
      //         forall( K, [&](const int i) {
      //               printf("%e %e %e %e\n", A(i,0,el), A(i,1,el), A(i,2,el), A(i,3,el));
      //            });
      //      };

      //   print_face( faceA, "faceA" );
      //   print_face( faceB, "faceB" );
      //   print_face( faceN, "faceN" );
      //}

      // Form diagonal terms (self element dependencies)
      auto& diag_coef = coefs[ blk_idx(el,el) ];

      diag_coef.set(0.0);

      forall( K, K, K, [&]( const int i, const int j, const int k ) {
            const int ij = i + j * K;
            const int kj = k + j * K;
            const int ik = i + k * K;

            auto sumx = [&]( const int l ) -> double { return Ax(l,j) * D(l,i) * C(k,l); };
            auto sumy = [&]( const int l ) -> double { return Ay(l,i) * D(l,j) * C(k,l); };

            diag_coef(kj,ij) += forall_reduce( K, sumx );
            diag_coef(ik,ij) += forall_reduce( K, sumy );

            forall( K, [&]( const int l ) {
                  const int lk = l + k * K;
                  diag_coef(lk,ij) -= (  Bx(l,j) * D(l,i) * C(k,j)
                                       + By(k,i) * D(k,j) * C(l,i) );
               });
         });

      // Add in neighbor terms (and boundary terms) in WESN order.
      for (int f = 0; f < 4; ++f)
      {
         const int LorR = f % 2; // 0 = left, 1 = right
         const int XorY = f / 2; // 0 is x, 1 is y;

         const auto nghbr_elem_id = elemNghborID(el)[ tensor2fem[f] ];

       //if (el < 10) printf("el: %d nghbr %d %d\n", el, nghbr_elem_id, f);

         if ( nghbr_elem_id < Nel )
         {
            auto& nghbr_coef = coefs[ blk_idx(el,nghbr_elem_id) ];
            nghbr_coef.set(0.0);

            /// Element neighbor on the f^th boundary.
            const int nghbr_face_id = nghbrFace[f];
            //printf("%d %d %d %d\n", f, nghbr_face_id, nghbr_elem_id, el);

            // Form the common geometric terms between the elements.
            typedef decltype(E00) Etype;
            typedef decltype(Ax)  Atype;
            auto kernel = [&] ( const Etype& Ebb, const Etype& Eab, const Etype& Fbb, const Etype& Fab, const Atype& A, const Atype& B)
               {
                  forall( K, K, K, [&] ( const int i, const int j, const int k )
                  {
                     auto Aj = faceA(j, nghbr_face_id, nghbr_elem_id);
                     auto Bj = 0.5 * ( faceB(j, f, el) + faceB(j, nghbr_face_id, nghbr_elem_id) );
                     //if (el == 0 and i == 0 and k == 0) printf("Aj %e %e\n", Aj, Bj);

                     const int ij = i + j * K + XorY * (K-1) * (i-j);
                     const int kj = k + j * K + XorY * (K-1) * (k-j);

                     {
                        auto sumx = [&] ( const int l ) -> double { return A(l,j) * D(l,i) * Eab(k,l); };

                        nghbr_coef(kj,ij) += (  forall_reduce( K, sumx )
                                              + Aj * Fab(k,i)
                                              + gLprime[LorR] * faceA(j,f,el) * Eab(k,i) );
                        diag_coef(kj,ij) += (  faceA(j,f,el) * Fbb(k,i)
                                             - gLprime[LorR] * Aj * Ebb(k,i) );
                     }

                     {
                        auto tmp = Bj * Basis.SolnNodesGradLgrangeBasis(k,j);
                        auto tmpy = tmp + B(i,k) * D(k,j);

                        forall( K, [&] ( const int l ) {
                              const int lk = l + k * K + XorY * (K-1) * (l-k);
                              nghbr_coef(lk,ij) -= tmpy * Eab(l,i);
                              diag_coef(lk,ij)  -= tmp  * Ebb(l,i);
                           });
                     }
                  });
               };

            if ( XorY == 0 )
               if ( LorR == 0 ) kernel( E00, E10, F00, F10, Ax, Bx );
               else             kernel( E11, E01, F11, F01, Ax, Bx );
            else
               if ( LorR == 0 ) kernel( E00, E10, F00, F10, Ay, By );
               else             kernel( E11, E01, F11, F01, Ay, By );
         }
         else
         {
            /// No element neighbor on the f^th surface. It's a boundary surface.
            /// Apply Neumann or Dirichlet conditions.

            const int bel = nghbr_elem_id - Nel;
          //if (el < 10) printf("bel: %d %d %d %d %d\n", f, nghbr_elem_id, el, bel, bndryType(bel));

            auto& bndry_src = bndrySrc[bel];

            if ( bndryType(bel) == BndryTypes::Dirichlet )
            {
               forall( K, K, [&] ( const int i, const int j ) {
                     const int ij = i + j * K + XorY * (K-1) * (i-j); // this transposes the storage for the x or y direction.
                     bndry_src(j,ij) += ( 2.0 * gLprime[LorR] * faceA(j,f,el) ) * Basis.NodesGradRadau(i,LorR);

                     forall( K, [&] ( const int k )
                        {
                           const int kj = k + j * K + XorY * (K-1) * (k-j);
                           const auto& A = ( XorY == 0 ) ? Ax : Ay;
                           const auto& B = ( XorY == 0 ) ? Bx : By;
                           const auto& Eaa = ( LorR == 0 ) ? E00 : E11;
                           const auto& Faa = ( LorR == 0 ) ? F00 : F11;

                           bndry_src(j,ij) += ( A(k,j) * D(k,i) * Basis.NodesGradRadau(k,LorR) );
                           bndry_src(k,ij) -= ( Basis.NodesGradRadau(i,LorR)
                                                  * (  B(i,k) * D(k,j)
                                                     + faceB(j,f,el) * Basis.SolnNodesGradLgrangeBasis(k,j) ) );

                           auto sumx = [&]( const int l ) { return A(l,j) * D(l,i) * Eaa(k,l); };

                           diag_coef(kj,ij) += ( 2.0 * faceA(j,f,el) * ( Faa(k,i) - gLprime[LorR] * Eaa(k,i) )
                                                  - forall_reduce( K, sumx ) );

                           forall( K, [&]( const int l ) {
                                 const int lk = l + k * K + XorY * (K-1) * (l-k);
                                 diag_coef(lk,ij) += ( B(i,k) * D(k,j) * Eaa(l,i) );
                              });
                        });
                  });
            }
            else if ( bndryType(bel) == BndryTypes::Neumann )
            {
               fprintf(stderr,"BC == Neumann not implemented\n");
               exit(1);
            }
            else
            {
               fprintf(stderr,"BC == default invalid\n");
               exit(1);
            }
         }
      }
   }

   auto t_end = getTimeStamp();
   std::cout << "Laplacian assembly 1 time: " << getElapsedTime( t_start, t_end ) << std::endl;

   if (false)
   {
      FILE *fp = fopen("cxx_assembly.out","w");

      for (int el(0); el < Nel; ++el)
      {
         fprintf(fp, "el: %d\n", el);

         {
            // insert the diagonal block.
            auto it = coefs.find( blk_idx(el,el) );
            if ( it == coefs.end() )
            {
               fprintf(stderr,"Diagonal block at %d not found\n", el);
               exit(1);
            }

            const auto& diag_coef = it->second;

            fprintf(fp,"diag: %d\n", el);
            forall( Ksq, [&](const int i) {
                  fprintf(fp,"%d,", i);
                  forall( Ksq, [&](const int j){ fprintf(fp,"%20.13e%c", diag_coef(i,j), (j == Ksq-1) ? '\n' : ','); });
               });
         }

         // Now add any neighbor terms.
         for (int f(0); f < 4; ++f)
         {
            const auto nghbr_elem_id = elemNghborID(el)[ tensor2fem[f] ];

            const auto idx = blk_idx(el,nghbr_elem_id);

            auto it = coefs.find( idx );
            if ( it == coefs.end() ) continue;

            const auto& nghbr_coef = it->second;
            const auto ij = blk_ij(it->first);

            fprintf(fp,"neigh: %d %d\n", f, nghbr_elem_id);

            forall( Ksq, [&](const int i) {
                  fprintf(fp,"%d,", i);
                  forall( Ksq, [&](const int j){ fprintf(fp,"%20.13e%c", nghbr_coef(i,j), (j == Ksq-1) ? '\n' : ','); });
               });
         }
      }

      fprintf(fp, "bndrySrc: \n");
      for (int bel = 0; bel < NelB; bel++)
      {
         auto& bndry_src = bndrySrc[bel];

         fprintf(fp, "bnd: %d\n", bel);
         forall( Ksq, [&](const int j) {
               fprintf(fp,"%d, ", j);
               forall( K, [&](const int i) { fprintf(fp,"%20.13e%c", bndry_src(i,j), (i == K-1) ? '\n' : ','); });
            });
      }

      fclose(fp);
   }

   // Convert block-matrix to scalar matrix for now.

   const size_t nrows = Nel * Ksq;
   const size_t maxnnz = Nel * 5 * (Ksq * Ksq);

   std::vector<int   > rowptr( nrows+1, 0 );
   std::vector<double> values; values.reserve( maxnnz );
   std::vector<int   > colidx; colidx.reserve( maxnnz );

   for (int el(0); el < Nel; ++el)
   {
      const int row0 = el * Ksq;

      // insert the diagonal block.
      auto it = coefs.find( blk_idx(el,el) );
      if ( it == coefs.end() )
      {
         fprintf(stderr,"Diagonal block at %d not found\n", el);
         exit(1);
      }

      const auto& diag_coef = it->second;

      for (int b_row(0); b_row < Ksq; b_row++)
      {
         rowptr[row0 + b_row] = values.size();

         // Add all central terms directly.
         for (int b_col(0); b_col < Ksq; b_col++)
         {
            values.push_back( diag_coef(b_col,b_row) ); // This seems transposed. Was the original storage backwards?
            colidx.push_back( row0 + b_col );
         }

         // Now add any neighbor terms.
         for (int f(0); f < 4; ++f)
         {
            const auto nghbr_elem_id = elemNghborID(el)[ tensor2fem[f] ];
            const auto nghbr_col0 = nghbr_elem_id * Ksq;

            const auto idx = blk_idx(el,nghbr_elem_id);

            auto nghbr_it = coefs.find( idx );
            if ( nghbr_it == coefs.end() ) continue;

            const auto& nghbr_coef = nghbr_it->second;

            for (int b_col(0); b_col < Ksq; b_col++)
            {
               values.push_back( nghbr_coef(b_col,b_row) ); // This seems transposed. Was the original storage backwards?
               colidx.push_back( nghbr_col0 + b_col );
            }
         }
      }
   }

   rowptr[nrows] = values.size();
   const size_t nnz = values.size();

   std::cout << "nrows: " << nrows << std::endl;
   std::cout << "nnz:   " << nnz << std::endl;

   // Construct the AMGCL solver/preconditioner.
   {
      auto A = std::tie( nrows, rowptr,
                                colidx,
                                values );

      amgcl_solver = std::make_shared<AMGCL_SolverType>( A );
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
      printf("norm: %e %e %e\n", norm, norm_b, norm / norm_b);

      eval_amgcl ( *csr, x, b, Knod );
   }

   return;
}

}
