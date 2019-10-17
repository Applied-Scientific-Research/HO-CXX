#ifndef __laplacian_hpp
#define __laplacian_hpp

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

namespace HighOrderFEM
{

//template <typename T>
//T sqr( const T& x ) { return x*x; }

// AMGCL interfaces.

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

/* Hold the preconditioner and iterative solver for the AMGCL library.
 */
struct AMGCL_SolverType
{
   typedef double value_type;

   typedef amgcl::backend::builtin<value_type> BackendType;		///< Primary back-end
   typedef amgcl::backend::builtin<float > mixed_BackendType;	///< Secondary back-end for the lower precision operations.

   typedef amgcl::runtime::preconditioner< mixed_BackendType > PreconditionerType;
   typedef amgcl::runtime::solver::wrapper< BackendType >      IterativeSolverType;
   typedef MixedPrecisionPrecon< PreconditionerType > MixedPreconditionerType;
   
   typedef amgcl::make_solver< PreconditionerType, IterativeSolverType > MakeSolverType;

   typedef typename MakeSolverType::build_matrix build_matrix;

   std::shared_ptr< MakeSolverType > shared_solver;
   std::shared_ptr< MixedPreconditionerType > shared_mixed_precon;

   build_matrix system_matrix; //!< A native copy of the system matrix since the internal
                               //!< storage will be lower precision.

   //!< Maximum # of iterations allowed by the iterative solver.
   int maxiters;

   /*! Target tolerance normalized by the RHS norm:
    * @f[
    * \frac{|| Ax - f ||}{||f||} < relTol
    * @f]
    */
   double reltol;

   /*! Target absolute tolerance.
    * @f[
    * || Ax - f || < absTol
    * @f]
    */
   double abstol;

   bool verbose;

   /// Constructor accepts any valid input matrix format.
   template <typename Matrix>
   AMGCL_SolverType ( const Matrix& A ) : system_matrix(A)
   {
      //std::cout << "typeid(build_matrix): " << typeid(build_matrix).name() << std::endl;

      auto prm = load_amgcl_params();

      this->maxiters = prm.get("solver.maxiter", 0 );
      this->reltol   = prm.get("solver.tol", 0.0 );
      this->abstol   = prm.get("solver.abstol", 0.0 );
      this->verbose  = prm.get("verbosity", 0 );

      auto t_build_start = getTimeStamp();

      shared_solver = std::make_shared< MakeSolverType >( system_matrix, prm );

      auto &P = shared_solver->precond();
      shared_mixed_precon = std::make_shared< MixedPreconditionerType > ( P );

      if ( this->verbose )
      {
         auto t_build_stop = getTimeStamp();

         std::cout << "Build time: " << getElapsedTime( t_build_start, t_build_stop ) << std::endl;

         std::cout << "Solver: " << std::endl << shared_solver->solver() << std::endl;
         std::cout << "Precon: " << std::endl << *shared_mixed_precon << std::endl;
      }
   }

   ~AMGCL_SolverType()
   {
      std::cout << "~AMGCL_SolverType" << std::endl;
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

      if (verbose)
      {
         size_t nrows = amgcl::backend::rows(A);
         std::vector<value_type> r(nrows);
         amgcl::backend::residual( f, A, x, r );
         auto norm_r = norm2( r, r );

         std::cout << "Iterations: " << niters << std::endl
                   << "Error:      " << resid << std::endl
                   << "Norm(r):    " << norm_r << std::endl
                   << "Solve time: " << getElapsedTime( t_start, t_stop ) << std::endl
                   << std::endl;
      }

      int iret = ( this->maxiters > 0 ) ? ( niters < this->maxiters ) : 1;
      return iret;
   }

   template <typename X, typename F>
   int operator()( X& x, const F& f )
   {
      return this->apply( x, f );
   }
};

} // namespace details

template <int _K, int _L, typename T = double>
struct LaplacianType
{
   typedef T value_type;
   enum : int { K = _K,
                L = _L };

   const GeometryType<K,L,value_type>& geometry;

   BasisFunctions<K,L,value_type> Basis;

   DynamicArrayType< StaticArrayType< value_type[K*K][K] > > bndrySrc;
   DynamicArrayType< StaticArrayType< value_type[K] > > bndryVal;

   std::shared_ptr< details::AMGCL_SolverType > amgcl_solver;

   LaplacianType( const GeometryType<K,L,value_type>& geo )
      : geometry(geo)
   {
      printf("Inside Laplacian::Laplacian\n");
      this->assemble();
   }

   void assemble(void)
   {
      using namespace details;

      const size_t Nel  = geometry.Nel;
      const size_t NelB = geometry.NelB;

      printf("Inside Laplacian::assemble\n");
      printf("K: %d\n", K);
      printf("Nel: %d\n", Nel);
      printf("NelB: %d\n", NelB);

      const int Ksq = K*K;

      typedef StaticArrayType< value_type[K] > K_ArrayType;
      typedef StaticArrayType< value_type[K][K] > KxK_ArrayType;
      typedef StaticArrayType< value_type[Ksq][Ksq] > KsqxKsq_ArrayType;

      this->bndrySrc.resize(NelB);

      for (int bel(0); bel < NelB; ++bel)
         bndrySrc(bel).set(0.0);

      // Build up some common basis function terms first.
      KxK_ArrayType E00, E01, E10, E11, E, C, D, F00, F01, F10, F11;

      const value_type gLprime[] = { -0.25*K*(K+1),
                                      0.25*K*(K+1) };

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

      std::map<size_t, KsqxKsq_ArrayType > coefs;

      auto t_start = getTimeStamp();

      for (int el(0); el < Nel; ++el)
      {
         const auto& metrics = geometry.getElementMetrics(el);

         // Shorthand to avoid object name.
         auto xxsi  = [&]( const int i, const int j ) { return metrics.xxsi(i,j); };
         auto yxsi  = [&]( const int i, const int j ) { return metrics.yxsi(i,j); };
         auto xeta  = [&]( const int i, const int j ) { return metrics.xeta(i,j); };
         auto yeta  = [&]( const int i, const int j ) { return metrics.yeta(i,j); };
         auto jac   = [&]( const int i, const int j ) { return metrics.jac (i,j); };
         auto faceA = [&]( const int i, const int f ) { return metrics.faceA(i,f); };
         auto faceB = [&]( const int i, const int f ) { return metrics.faceB(i,f); };
         auto faceN = [&]( const int i, const int f ) { return metrics.faceN(i,f); };

         // 2nd-order volumetric and boundary metric terms

         KxK_ArrayType Ax, Ay, Bx, By;

         forall( K, K, [&] ( const int i, const int j ) {
               Ax(i,j) = ( xeta(i,j)*xeta(i,j) + yeta(i,j)*yeta(i,j) ) / jac(i,j);
               Ay(j,i) = ( xxsi(i,j)*xxsi(i,j) + yxsi(i,j)*yxsi(i,j) ) / jac(i,j); // tranposed indexes!!!
               Bx(i,j) = ( xxsi(i,j)*xeta(i,j) + yxsi(i,j)*yeta(i,j) ) / jac(i,j);
               By(j,i) = Bx(i,j); // tranposed only
            });

         // Form diagonal terms (self element dependencies)
         auto& diag_coef = coefs[ blk_idx(el,el) ];

         diag_coef.set(0.0);

         forall( K, K, K, [&]( const int i, const int j, const int k )
            {
               const int ij = i + j * K;
               const int kj = k + j * K;
               const int ik = i + k * K;

               auto sumx = [&]( const int l ) -> value_type { return Ax(l,j) * D(l,i) * C(k,l); };
               auto sumy = [&]( const int l ) -> value_type { return Ay(l,i) * D(l,j) * C(k,l); };

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

            const auto nghbr_elem_id = geometry.getNeighborID(el,f);

            const auto& nghbr_metrics = geometry.getElementMetrics(nghbr_elem_id);

            if ( geometry.isBoundaryElement( nghbr_elem_id ) )
            {
               auto& nghbr_coef = coefs[ blk_idx(el,nghbr_elem_id) ];
               nghbr_coef.set(0.0);

               /// Element neighbor on the f^th boundary.
               const int nghbr_face_id = nghbrFace[f];

               // Form the common geometric terms between the elements.
               typedef decltype(E00) Etype;
               typedef decltype(Ax)  Atype;
               auto kernel = [&] ( const Etype& Ebb, const Etype& Eab, const Etype& Fbb, const Etype& Fab, const Atype& A, const Atype& B)
                  {
                     forall( K, K, K, [&] ( const int i, const int j, const int k )
                     {
                        auto Aj = nghbr_metrics.faceA(j, nghbr_face_id);
                        auto Bj = 0.5 * ( faceB(j, f) + nghbr_metrics.faceB(j, nghbr_face_id) );

                        const int ij = i + j * K + XorY * (K-1) * (i-j);
                        const int kj = k + j * K + XorY * (K-1) * (k-j);

                        {
                           auto sumx = [&] ( const int l ) -> value_type { return A(l,j) * D(l,i) * Eab(k,l); };

                           nghbr_coef(kj,ij) += (  forall_reduce( K, sumx )
                                                 + Aj * Fab(k,i)
                                                 + gLprime[LorR] * faceA(j,f) * Eab(k,i) );
                           diag_coef(kj,ij) += (  faceA(j,f) * Fbb(k,i)
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

               const auto bel = geometry.getBoundaryID( nghbr_elem_id );
               const BndryType bndry_type = geometry.getBoundaryType(bel);

               auto& bndry_src = bndrySrc[bel];

               if ( bndry_type == BndryType::Dirichlet )
               {
                  forall( K, K, [&] ( const int i, const int j ) {
                        const int ij = i + j * K + XorY * (K-1) * (i-j); // this transposes the storage for the x or y direction.
                        bndry_src(j,ij) += ( 2.0 * gLprime[LorR] * faceA(j,f) ) * Basis.NodesGradRadau(i,LorR);

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
                                                        + faceB(j,f) * Basis.SolnNodesGradLgrangeBasis(k,j) ) );

                              auto sumx = [&]( const int l ) { return A(l,j) * D(l,i) * Eaa(k,l); };

                              diag_coef(kj,ij) += ( 2.0 * faceA(j,f) * ( Faa(k,i) - gLprime[LorR] * Eaa(k,i) )
                                                     - forall_reduce( K, sumx ) );

                              forall( K, [&]( const int l ) {
                                    const int lk = l + k * K + XorY * (K-1) * (l-k);
                                    diag_coef(lk,ij) += ( B(i,k) * D(k,j) * Eaa(l,i) );
                                 });
                           });
                     });
               }
               else if ( bndry_type == BndryType::Neumann )
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

      // Convert block-matrix to scalar matrix for now.

      const size_t nrows = Nel * Ksq;
      const size_t maxnnz = Nel * 5 * (Ksq * Ksq);

      std::vector<int> rowptr( nrows+1, 0 );
      std::vector<value_type> values;
      std::vector<int> colidx;

      values.reserve( maxnnz );
      colidx.reserve( maxnnz );

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
               const auto nghbr_elem_id = geometry.getNeighborID(el,f);
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
      auto A = std::tie( nrows, rowptr,
                                colidx,
                                values );

      this->amgcl_solver = std::make_shared<AMGCL_SolverType>( A );

      return;
   }

   void setBoundaryValues( const size_t inNelB, value_type inBoundaryValuesPtr[] )
   {
      const size_t NelB = geometry.NelB;
      if ( inNelB != NelB )
      {
         std::cerr << "NelB does not match GeometryType value: " << inNelB << ", " << NelB << std::endl;
         exit(-1);
      }

      this->bndryVal.resize( NelB );

      typedef decltype( *bndryVal.getPointer() ) bndryValType;
      const decltype( bndryVal ) inBoundaryValues( inBoundaryValuesPtr, NelB );

      for (int bel(0); bel < NelB; ++bel)
         bndryVal(bel).copy( inBoundaryValues[bel] );
   }

 //void solve( value_type VorticityIn[] )
   template < typename VorticityType, typename PsiType >
   void solve( const VorticityType& Vorticity, PsiType& psi )
   {
      const size_t Nel  = geometry.Nel;
      const size_t NelB = geometry.NelB;

      printf("Inside Laplacian::solve\n");
      printf("K: %d\n", K);
      printf("Nel: %d\n", Nel);
      printf("NelB: %d\n", NelB);

      typedef StaticArrayType< value_type[K][K] > NodeArrayType;

    //const DynamicArrayType< NodeArrayType > Vorticity( (NodeArrayType *) VorticityIn, Nel );
    //const DynamicArrayType< NodeArrayType > psi( (NodeArrayType *) psiIn, Nel );

      const size_t nrows = Nel * (K*K);
      std::vector<value_type> xflat(nrows), bflat(nrows);

      DynamicArrayType< NodeArrayType > x( (NodeArrayType*) xflat.data(), Nel),
                                        b( (NodeArrayType*) bflat.data(), Nel);

      auto t_start = getTimeStamp();

      /// Evaluate the RHS vorticity for each element's nodes.

      #pragma omp parallel for
      for (int el = 0; el < Nel; ++el)
      {
         const auto& vort = Vorticity(el);
         const auto& metrics = geometry.getElementMetrics(el);

         forall( K, K, [&]( const int i, const int j )
            {
               b[el](i,j) = -metrics.jac(i,j) * vort(i,j);
             //x[el](i,j) = 0.0;
            });
      }

      /// Factor in the boundary face/element.

      for (int bel = 0; bel < NelB; ++bel)
      {
         /// Volumetric element ID
         const auto el_id = geometry.bndryElementID(bel);

         /// References to the element's objects.
               auto& b_el = b[el_id];
         const auto& bndry_src = this->bndrySrc[bel];
         const auto& bndry_val = this->bndryVal[bel];

         forall( K, K, [&] ( const int i, const int j )
            {
               const auto ij = i + j * K;

               b_el(i,j) -= dot_product( bndry_src.template slice<0>( ij ), bndry_val );
            });
      }

      auto t_middle = getTimeStamp();

      int ierr = (*this->amgcl_solver)( xflat, bflat );

      std::copy( xflat.data(), xflat.data() + nrows, psi.getRawPointer() );

      auto t_end = getTimeStamp();

      printf("cxx solved: %d %f %f\n", ierr, getElapsedTime( t_start, t_middle ), getElapsedTime( t_middle, t_end ));

      return;
   }

};

}

#endif
