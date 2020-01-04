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
#include <vector>

#include "wtimer.hpp"
#include "memory.hpp"
#include "array.hpp"
#include "basis.hpp"
#include "geometry.hpp"

#include "linsolvers/linsolver.hpp"

namespace HighOrderFEM
{

namespace details
{

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

   std::shared_ptr< LinearSolver::BaseLinearSolver > linear_solver;

#ifdef ENABLE_LINEAR_SOLVER_BENCHMARK
   // For an internal benchmark ...
   std::vector< std::shared_ptr< LinearSolver::BaseLinearSolver > > all_linear_solvers;
#endif

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

/*#ifdef ENABLE_AMGCL

      // Construct the AMGCL solver/preconditioner.
      //auto A = std::tie( nrows, rowptr,
      //                          colidx,
      //                          values );

      //this->amgcl_solver = std::make_shared< amgcl_solver_type >( A );
      //this->amgcl_solver = std::make_shared< amgcl_solver_type >( );
      //this->amgcl_solver->build( nrows, rowptr, colidx, values );
      this->amgcl_solver.build( nrows, rowptr, colidx, values );
#endif*/

      static bool write_matrix = true;
      if (write_matrix) {
         FILE *fp = fopen("a.out","wb");
         const int n = nrows;
         fwrite( &n, sizeof(int), 1, fp );
         fwrite( rowptr.data(), sizeof(int), (nrows+1), fp );
         fwrite( colidx.data(), sizeof(int), nnz, fp );
         fwrite( values.data(), sizeof(double), nnz, fp );
         fclose(fp);
         write_matrix = false;
         std::cout << "Wrote CSR matrix to binary file\n";
      }

/*#ifdef ENABLE_EIGEN
      {
         this->eigen_solver.build( nrows, rowptr, colidx, values );
      }
#endif
#ifdef ENABLE_HYPRE
      {
         //this->hypre_solver.setVerbosity(1);
         this->hypre_solver.build( nrows, rowptr, colidx, values );
         this->linear_solver = &this->hypre_solver;
      }
#endif
#ifdef ENABLE_APLLES
      {
         this->aplles_solver.build( nrows, rowptr, colidx, values );
      }
#endif*/

      this->linear_solver = LinearSolver::CreateLinearSolver();
      printf("Instantiated linear solver %s\n", linear_solver->name().c_str());

      this->linear_solver->build( nrows, rowptr, colidx, values );
      printf("Built linear solver %s\n", linear_solver->name().c_str());

#ifdef ENABLE_LINEAR_SOLVER_BENCHMARK
      this->all_linear_solvers = LinearSolver::CreateAllLinearSolvers();

      for (auto solver : this->all_linear_solvers)
      {
         printf("Instantiated linear solver %s for benchmark\n", solver->name().c_str());

         solver->build( nrows, rowptr, colidx, values );
         printf("Built linear solver %s for benchmark\n", solver->name().c_str());
      }
#endif

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
   template < typename VectorType >
   VectorType solve( const VectorType& Vorticity )
   {
      const size_t Nel  = geometry.Nel;
      const size_t NelB = geometry.NelB;

      printf("Inside Laplacian::solve\n");
      printf("K: %d\n", K);
      printf("Nel: %d\n", Nel);
      printf("NelB: %d\n", NelB);

      typedef StaticArrayType< value_type[K][K] > NodeArrayType;

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

      //int ierr = (*this->amgcl_solver)( xflat, bflat );
      //int ierr = this->amgcl_solver->solve( bflat, xflat );
      //int ierr = this->amgcl_solver.solve( bflat, xflat );
      int ierr = this->linear_solver->solve( bflat, xflat );

      VectorType psi(Nel);

      std::copy( xflat.data(), xflat.data() + nrows, psi.getRawPointer() );

      auto t_end = getTimeStamp();

      printf("linear solver %s time: %d %f %f (sec)\n", linear_solver->name().c_str(), ierr, getElapsedTime( t_start, t_middle ), getElapsedTime( t_middle, t_end ));

#ifdef ENABLE_LINEAR_SOLVER_BENCHMARK
      {
         static bool finished = false;
         if (not(finished))
         {
            for (auto solver : this->all_linear_solvers)
            {
               auto t_begin = getTimeStamp();

               std::vector<value_type> _x(nrows);

               const int nsolves = 10;

               int ret = LinearSolver::SolverStatusFlags::Success;
               for (int iter = 0; iter < nsolves; iter++)
               {
                  std::fill(_x.begin(), _x.end(), 0.0);
                  ret = solver->solve( bflat, _x );
                  if ( ret != LinearSolver::SolverStatusFlags::Success )
                     break;
               }

               auto t_end = getTimeStamp();

               if ( ret != LinearSolver::SolverStatusFlags::Success )
                  printf("solver %s failed to solve the problem %d\n", solver->name().c_str(), ret);
               else
                  printf("solver %s finished in %f (sec)\n", solver->name().c_str(),
                                                                   getElapsedTime( t_begin, t_end ) / nsolves);
            }
            finished = true;

            this->all_linear_solvers.clear();
         }
      }
#endif

      {
         static bool write_rhs = true;
         if (write_rhs) {
            FILE *fp = fopen("b.out","wb");
            const int n = nrows;
            fwrite( &n, sizeof(int), 1, fp );
            fwrite( bflat.data(), sizeof(double), nrows, fp );
            fclose(fp);
            write_rhs = false;
         }
      }

/*#ifdef ENABLE_AMGCL
      if (linear_solver != &amgcl_solver)
      {
         auto t_begin = getTimeStamp();

         std::vector<value_type> _x(nrows);
         int ierr = this->amgcl_solver.solve( bflat, _x );

         auto t_end = getTimeStamp();
         printf("amgcl solver time: %d %f (sec)\n", ierr, getElapsedTime( t_begin, t_end ));

         double err2 = 0., ref2 = 0.;
         for (int i = 0; i < nrows; ++i)
         {
            double diff = xflat[i] - _x[i];
            err2 += diff*diff;
            ref2 += xflat[i]*xflat[i];
         }
         printf("amgcl regtest: %e %e %e\n", sqrt(err2), sqrt(ref2), sqrt(err2/ref2));
      }
#endif

#ifdef ENABLE_EIGEN
      if (linear_solver != &eigen_solver)
      {
         auto t_begin = getTimeStamp();

         std::vector<value_type> _x(nrows);
         int ierr = this->eigen_solver.solve( bflat, _x );

         auto t_end = getTimeStamp();
         printf("eigen solver time: %d %f (sec)\n", ierr, getElapsedTime( t_begin, t_end ));

         double err2 = 0., ref2 = 0.;
         for (int i = 0; i < nrows; ++i)
         {
            double diff = xflat[i] - _x[i];
            err2 += diff*diff;
            ref2 += xflat[i]*xflat[i];
         }
         printf("eigen regtest: %e %e %e\n", sqrt(err2), sqrt(ref2), sqrt(err2/ref2));
      }
#endif

#ifdef ENABLE_HYPRE
      if (linear_solver != &hypre_solver)
      {
         auto t_begin = getTimeStamp();

         std::vector<value_type> _x(nrows);
         int ierr = this->hypre_solver.solve( bflat, _x );

         auto t_end = getTimeStamp();
         printf("hypre solver time: %d %f (sec)\n", ierr, getElapsedTime( t_begin, t_end ));

         double err2 = 0., ref2 = 0.;
         for (int i = 0; i < nrows; ++i)
         {
            double diff = xflat[i] - _x[i];
            err2 += diff*diff;
            ref2 += xflat[i]*xflat[i];
         }
         printf("hypre regtest: %e %e %e\n", sqrt(err2), sqrt(ref2), sqrt(err2/ref2));
      }
#endif

#ifdef ENABLE_APLLES
      if (linear_solver != &aplles_solver)
      {
         auto t_begin = getTimeStamp();

         std::vector<value_type> _x(nrows);
         int ierr = this->aplles_solver.solve( bflat, _x );

         auto t_end = getTimeStamp();
         printf("aplles solver time: %d %f (sec)\n", ierr, getElapsedTime( t_begin, t_end ));

         double err2 = 0., ref2 = 0.;
         for (int i = 0; i < nrows; ++i)
         {
            double diff = xflat[i] - _x[i];
            err2 += diff*diff;
            ref2 += xflat[i]*xflat[i];
         }
         printf("aplles regtest: %e %e %e\n", sqrt(err2), sqrt(ref2), sqrt(err2/ref2));
      }
#endif */

      return psi;
   }

};

}

#endif
