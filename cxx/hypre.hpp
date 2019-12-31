#ifndef __hypre_hpp
#define __hypre_hpp

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <sstream>
#include <algorithm>

#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "seq_mv.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_parcsr_ls.h"

#include "linsolver.hpp"

namespace HighOrderFEM
{
namespace LinearSolver
{

namespace hypre_details
{

decltype(hypre_MPI_COMM_WORLD) mpi_comm = hypre_MPI_COMM_WORLD;
int mpi_nprocs = 1, mpi_id = 0;
int ilower = -1, iupper = -1;
int global_num_rows = -1, local_num_rows = -1;

struct HypreSolverWrapper
{
   enum Tags : int { None, GMRES, FlexGMRES, BiCGSTAB, BoomerAMG, DiagScaling };

   const char* getName( Tags tag )
   {
      switch( tag )
      {
         case(GMRES): return "GMRES";
         case(FlexGMRES): return "FlexGMRES";
         case(BiCGSTAB): return "BiCGSTAB";
         case(BoomerAMG): return "BoomerAMG";
         case(DiagScaling): return "DiagScaling";
      }
   }
   const char* getName(void) { return this->getName( this->tag ); }

   Tags getTag( const int id )
   {
      switch ( id )
      {
         case(GMRES): return GMRES;
         case(FlexGMRES): return FlexGMRES;
         case(BiCGSTAB): return BiCGSTAB;
         case(BoomerAMG): return BoomerAMG;
         case(DiagScaling): return DiagScaling;
         default: {
            fprintf(stderr,"Unknown Solver id %d\n", id);
            exit(1);
         }
      }
   }

   typedef HYPRE_Int (*DestructorFcn)( HYPRE_Solver solver );
   typedef HYPRE_Int (*GetNumIterationsFcn)( HYPRE_Solver solver, int* );

   Tags tag;
   HYPRE_Solver solver_ptr;

   HYPRE_PtrToSolverFcn solve_func_ptr;
   HYPRE_PtrToSolverFcn setup_func_ptr;
   DestructorFcn destr_func_ptr;
   GetNumIterationsFcn niter_func_ptr;

   //HypreSolverWrapper(void)
   //   : tag(None), solver_ptr(NULL),
   //     solve_func_ptr(NULL), setup_func_ptr(NULL), destr_func_ptr(NULL), niter_func_ptr(NULL)
   //{}

   HypreSolverWrapper( Tags tag, HYPRE_Solver solver_ptr,
                           HYPRE_PtrToSolverFcn solve_func_ptr,
                           HYPRE_PtrToSolverFcn setup_func_ptr,
                           DestructorFcn destr_func_ptr,
                           GetNumIterationsFcn niter_func_ptr )
      : tag(tag), solver_ptr(solver_ptr),
        solve_func_ptr(solve_func_ptr), setup_func_ptr(setup_func_ptr),
        destr_func_ptr(destr_func_ptr), niter_func_ptr(niter_func_ptr)
   {}

   ~HypreSolverWrapper()
   {
      std::cout << "Destroying HypreSolverWrapper: " << this << " with tag: " << int(this->tag) << " name: " << this->getName() << "\n";
      if ( this->solver_ptr )
      {
         this->destr_func_ptr( this->solver_ptr );
         this->solver_ptr = NULL;
      }
   }

   int getNumIterations(void)
   {
      if (this->niter_func_ptr)
      {
         int n;
         this->niter_func_ptr( this->solver_ptr, &n );
         return n;
      }
      else
         return -1;
   }
};

std::shared_ptr<HypreSolverWrapper>
build_BoomerAMG( std::vector<std::string>& args )
{
   HYPRE_Solver solver;

   // Defaults
   int coarsen_type = 10; // hmis
   int coarse_threshold = 500; // this is small enough.
   int npre = 1, npost = 1;
   int relax_type[] = { 0, 13, 14, 9 }; // all, down, up, coarse solvers
   double relax_wt = 1.0;
   int interp_type = 6; // extended classical
   double theta = 0.25;
   double trunc_factor = 0.; // truncation factor for coarsening.
   double agg_trunc_factor = 0.; // aggressive truncation factor for interpolation.

   int argc = args.size();
   for (int arg_index = 0; arg_index < argc; arg_index++ )
   {
      std::cout << args[arg_index] << std::endl;

      if ( args[arg_index].compare("-coarsen_type") == 0 )
      {
         coarsen_type = atoi( args[ ++arg_index ].c_str() );
      }
      else if ( args[arg_index].compare("-cljp") == 0 )
         coarsen_type      = 0;
      else if ( args[arg_index].compare("-cljp1") == 0 )
         coarsen_type      = 7;
      else if ( args[arg_index].compare("-pmis") == 0 )
         coarsen_type      = 8;
      else if ( args[arg_index].compare("-pmis1") == 0 )
         coarsen_type      = 9;
      else if ( args[arg_index].compare("-hmis") == 0 )
         coarsen_type      = 10;
      else if ( args[arg_index].compare("-ruge") == 0 )
         coarsen_type      = 1;
      else if ( args[arg_index].compare("-ruge1p") == 0 )
         coarsen_type      = 11;
      else if ( args[arg_index].compare("-ruge2b") == 0 )
         coarsen_type      = 2;
      else if ( args[arg_index].compare("-ruge3") == 0 )
         coarsen_type      = 3;
      else if ( args[arg_index].compare("-ruge3c") == 0 )
         coarsen_type      = 4;
      else if ( args[arg_index].compare("-rugerlx") == 0 )
         coarsen_type      = 5;
      else if ( args[arg_index].compare("-falgout") == 0 )
         coarsen_type      = 6;
      else if ( args[arg_index].compare("-coarse_enough") == 0 )
      {
         coarse_threshold = atoi(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-npre") == 0 )
      {
         npre = atoi(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-npost") == 0 )
      {
         npost = atoi(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-relax_type") == 0 )
      {
         relax_type[0] = atoi(args[ ++arg_index ].c_str());
         relax_type[1] = relax_type[0];
         relax_type[2] = relax_type[0];
      }
      else if ( args[arg_index].compare("-relax_wt") == 0 )
      {
         relax_wt = atof(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-interp_type") == 0 )
      {
         interp_type = atoi(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-theta") == 0 )
      {
         theta = atof(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-trunc_factor") == 0 )
      {
         trunc_factor = atof(args[ ++arg_index ].c_str());
      }
      else if ( args[arg_index].compare("-agg_trunc_factor") == 0 )
      {
         agg_trunc_factor = atof(args[ ++arg_index ].c_str());
      }
   }

   /* Now set up the AMG preconditioner and specify any parameters */
   HYPRE_BoomerAMGCreate(&solver);
   HYPRE_BoomerAMGSetPrintLevel(solver, 1); /* print amg solution info */
   HYPRE_BoomerAMGSetCoarsenType(solver, coarsen_type);
   HYPRE_BoomerAMGSetTol(solver, 0.0); /* conv. tolerance zero */
   HYPRE_BoomerAMGSetMaxIter(solver, 1); /* do only one iteration! */
   HYPRE_BoomerAMGSetMaxCoarseSize(solver, coarse_threshold);
   HYPRE_BoomerAMGSetInterpType(solver, interp_type);
   HYPRE_BoomerAMGSetStrongThreshold( solver, theta ); // strong threshold

   HYPRE_BoomerAMGSetTruncFactor(solver, trunc_factor);
   HYPRE_BoomerAMGSetAggTruncFactor(solver, agg_trunc_factor);

   HYPRE_BoomerAMGSetNumSweeps(solver, 1);
   HYPRE_BoomerAMGSetCycleNumSweeps(solver, 1, 3); // # coarse solves
   HYPRE_BoomerAMGSetCycleNumSweeps(solver, npre, 1); // # pre-smoother solves
   HYPRE_BoomerAMGSetCycleNumSweeps(solver, npost, 2); // # post-smoother solves

   for (int i = 0; i < 4; ++i)
      HYPRE_BoomerAMGSetCycleRelaxType(solver, relax_type[i], i);
   HYPRE_BoomerAMGSetRelaxWt(solver, relax_wt);

   {
      int ct;
      HYPRE_BoomerAMGGetCoarsenType( solver, &ct );
      printf("ct: %d %d\n", ct, coarsen_type);

      int rt[4];
      for (int i = 0; i < 4; ++i)
         HYPRE_BoomerAMGGetCycleRelaxType(solver, &rt[i], i);
      printf("rt: %d %d %d %d\n", rt[0], rt[1], rt[2], rt[3]);
   }

   return std::make_shared<HypreSolverWrapper>( HypreSolverWrapper::BoomerAMG,
                              solver,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                              HYPRE_BoomerAMGDestroy,
                              HYPRE_BoomerAMGGetNumIterations );
}

HypreSolverWrapper build_DiagScaling( std::vector<std::string>& args )
{
   return HypreSolverWrapper( HypreSolverWrapper::DiagScaling,
                              NULL,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
                              NULL, NULL );
}

HypreSolverWrapper build_GMRES( std::vector<std::string>& args, HypreSolverWrapper *precon, BaseLinearSolver *base)
{
   HYPRE_Solver solver;

   HYPRE_ParCSRGMRESCreate( hypre_MPI_COMM_WORLD, &solver );

   HYPRE_GMRESSetKDim(solver, base->restart_k);
   HYPRE_GMRESSetMaxIter(solver, base->maxiters);
   HYPRE_GMRESSetTol(solver, base->reltol);
   HYPRE_GMRESSetAbsoluteTol(solver, base->abstol);
   HYPRE_GMRESSetLogging(solver, int(base->verbosity > 0) );
   HYPRE_GMRESSetPrintLevel(solver, base->verbosity );

   if (precon != NULL)
      HYPRE_GMRESSetPrecond( solver,
                             (HYPRE_PtrToSolverFcn) precon->solve_func_ptr,
                             (HYPRE_PtrToSolverFcn) precon->setup_func_ptr,
                             (HYPRE_Solver) precon->solver_ptr );

   return HypreSolverWrapper( HypreSolverWrapper::GMRES, solver,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRGMRESSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRGMRESSetup,
                              HYPRE_ParCSRGMRESDestroy,
                              HYPRE_GMRESGetNumIterations );
}

std::shared_ptr<HypreSolverWrapper>
build_FlexGMRES( std::vector<std::string>& args, std::shared_ptr<HypreSolverWrapper> precon, BaseLinearSolver *base)
{
   HYPRE_Solver solver;

   HYPRE_ParCSRFlexGMRESCreate( hypre_MPI_COMM_WORLD, &solver );

   printf("created HYPRE_ParCSRFlexGMRES\n");

   HYPRE_FlexGMRESSetKDim(solver, base->restart_k);
   HYPRE_FlexGMRESSetMaxIter(solver, base->maxiters);
   HYPRE_FlexGMRESSetTol(solver, base->reltol);
   HYPRE_FlexGMRESSetAbsoluteTol(solver, base->abstol);
   HYPRE_FlexGMRESSetLogging(solver, int(base->verbosity > 0) );
   HYPRE_FlexGMRESSetPrintLevel(solver, base->verbosity );

   if (precon != NULL)
      HYPRE_FlexGMRESSetPrecond( solver,
                             (HYPRE_PtrToSolverFcn) precon->solve_func_ptr,
                             (HYPRE_PtrToSolverFcn) precon->setup_func_ptr,
                             (HYPRE_Solver) precon->solver_ptr );

   return std::make_shared<HypreSolverWrapper>( HypreSolverWrapper::FlexGMRES, solver,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRFlexGMRESSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRFlexGMRESSetup,
                              HYPRE_ParCSRFlexGMRESDestroy,
                              HYPRE_FlexGMRESGetNumIterations );
}

HypreSolverWrapper build_BiCGSTAB( std::vector<std::string>& args, HypreSolverWrapper *precon, BaseLinearSolver *base)
{
   HYPRE_Solver solver;

   HYPRE_ParCSRBiCGSTABCreate( hypre_MPI_COMM_WORLD, &solver);

   HYPRE_BiCGSTABSetMaxIter(solver, base->maxiters);
   HYPRE_BiCGSTABSetTol(solver, base->reltol);
   HYPRE_BiCGSTABSetAbsoluteTol(solver, base->abstol);
   HYPRE_BiCGSTABSetLogging(solver, int(base->verbosity > 0) );
   HYPRE_BiCGSTABSetPrintLevel(solver, base->verbosity );

   if (precon != NULL)
      HYPRE_BiCGSTABSetPrecond( solver,
                             (HYPRE_PtrToSolverFcn) precon->solve_func_ptr,
                             (HYPRE_PtrToSolverFcn) precon->setup_func_ptr,
                             (HYPRE_Solver) precon->solver_ptr );

   return HypreSolverWrapper( HypreSolverWrapper::BiCGSTAB, solver,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRBiCGSTABSolve,
                              (HYPRE_PtrToSolverFcn) HYPRE_ParCSRBiCGSTABSetup,
                              HYPRE_ParCSRBiCGSTABDestroy,
                              HYPRE_BiCGSTABGetNumIterations );
}

} // hypre_details

struct HypreSolver : BaseLinearSolver
{
   static_assert( std::is_same< value_type, double >::value, "Only double support by Hypre" );
   static_assert( std::is_same< index_type, int >::value, "Only int support by Hypre" );

   std::shared_ptr<hypre_details::HypreSolverWrapper> solver;
   std::shared_ptr<hypre_details::HypreSolverWrapper> precon;

   HYPRE_IJMatrix ij_A;
   HYPRE_IJVector ij_x;
   HYPRE_IJVector ij_b;
   HYPRE_IJVector ij_r;

   HYPRE_ParCSRMatrix  parcsr_A;
   HYPRE_ParVector     parvec_b;
   HYPRE_ParVector     parvec_x;
   HYPRE_ParVector     parvec_r;

   HypreSolver(void)
      : BaseLinearSolver( LinearSolverTag::HYPRE ),
        ij_A(NULL), ij_x(NULL), ij_b(NULL), ij_r(NULL)
   {}

   ~HypreSolver()
   {
      std::cout << "Destroying HypreSolver: " << this << "\n";
      if (ij_A) HYPRE_IJMatrixDestroy( ij_A );
      if (ij_x) HYPRE_IJVectorDestroy( ij_x );
      if (ij_b) HYPRE_IJVectorDestroy( ij_b );
      if (ij_r) HYPRE_IJVectorDestroy( ij_r );
   }

   HYPRE_IJMatrix build_CSRMatrix( const index_type num_rows,
                                   const std::vector<index_type>& rowptr,
                                   const std::vector<index_type>& colidx,
                                   const std::vector<value_type>& values )
   {
      using namespace hypre_details;

      int num_nonzeros = rowptr[ num_rows ];
      printf("%d %d\n", num_rows, num_nonzeros);

      ilower = 0;
      iupper = num_rows - 1;
      local_num_rows = iupper - ilower + 1;
      global_num_rows = num_rows;

      HYPRE_IJMatrix A;

      HYPRE_IJMatrixCreate( hypre_MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A );

      /* Choose a parallel csr format storage (see the User's Manual) */
      HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

      /* Initialize before setting coefficients */
      HYPRE_IJMatrixInitialize(A);

      for (int i = ilower; i <= iupper; ++i)
      {
         int k0 = rowptr[i];
         int nnz = rowptr[i+1] - k0;
         double *data_i = const_cast<double*>( values.data() ) + k0;
         int *cols_i = const_cast<int*>( colidx.data() ) + k0;

         HYPRE_IJMatrixSetValues( A, 1, &nnz, &i, cols_i, data_i );
      }

      /* Assemble after setting the coefficients */
      HYPRE_IJMatrixAssemble(A);

      return A;
   }

   double residual_2norm (void)
   {
      hypre_ParVectorCopy( parvec_b, parvec_r );
      hypre_ParCSRMatrixMatvec( -1.0, parcsr_A, parvec_x, 1.0, parvec_r );

      return vector_2norm( parvec_r );
   }

   double vector_2norm( HYPRE_ParVector &vec )
   {
      double norm;
      HYPRE_ParVectorInnerProd( vec, vec, &norm );
      return sqrt( norm );
   }

   int build( const index_type nrows, const std::vector<index_type>& rowptr,
                                      const std::vector<index_type>& colidx,
                                      const std::vector<value_type>& values )
   {
      using namespace hypre_details;
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      auto t_start = getTimeStamp();

      this->num_rows = nrows;

#ifndef HYPRE_SEQUENTIAL
#error 'Assuming HYPRE is built w/o MPI support and only running on 1 MPI rank'
#endif

      /* Initialize MPI */
      hypre_MPI_Init( NULL, NULL );

      hypre_MPI_Comm_size( hypre_MPI_COMM_WORLD, &mpi_nprocs );
      hypre_MPI_Comm_rank( hypre_MPI_COMM_WORLD, &mpi_id );

      /* Initialize Hypre */
      HYPRE_Init( 0, NULL );

      this->ij_A = build_CSRMatrix( nrows, rowptr, colidx, values );

      // Create solution and residual vectors.
      std::vector<double> x_values( local_num_rows );
      for (int i = 0; i < local_num_rows; i++)
          x_values[i] = 0.;

      HYPRE_IJVectorCreate( hypre_MPI_COMM_WORLD, ilower, iupper, &ij_x );
      HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_x);
      HYPRE_IJVectorSetValues( ij_x, local_num_rows, NULL, x_values.data() );
      HYPRE_IJVectorAssemble( ij_x );

      HYPRE_IJVectorCreate( hypre_MPI_COMM_WORLD, ilower, iupper, &ij_r );
      HYPRE_IJVectorSetObjectType(ij_r, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_r);
      HYPRE_IJVectorSetValues( ij_r, local_num_rows, NULL, x_values.data() );
      HYPRE_IJVectorAssemble( ij_r );

      HYPRE_IJVectorCreate( hypre_MPI_COMM_WORLD, ilower, iupper, &ij_b );
      HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_b);
      HYPRE_IJVectorSetValues( ij_b, local_num_rows, NULL, x_values.data() );
      HYPRE_IJVectorAssemble( ij_b );

#ifdef _OPENMP
      HYPRE_IJMatrixSetOMPFlag( this->ij_A, 1 );
#endif

      /* Get the parcsr matrix/vector objects to use */
      HYPRE_IJMatrixGetObject( ij_A, (void**) &parcsr_A);
      HYPRE_IJVectorGetObject( ij_x, (void**) &parvec_x );
      HYPRE_IJVectorGetObject( ij_b, (void**) &parvec_b );
      HYPRE_IJVectorGetObject( ij_r, (void**) &parvec_r );

      // Default is FlexGMRES with BoomerAMG.
      std::string options = "-pmis -interp_type 7 -theta 0.6 -relax_type 0 -relax_wt 0.533333 -npost 2 -trunc_factor 0.25";
      std::vector< std::string > args;

      {
         std::stringstream ss( options );
         std::string token;
         while ( std::getline(ss, token, ' ') )
            args.push_back(token);
      }

      this->precon = build_BoomerAMG( args );
      this->solver = build_FlexGMRES( args, precon, static_cast< BaseLinearSolver* >(this) );

      printf("instantiated objects ... calling setup ... \n");

      this->solver->setup_func_ptr( solver->solver_ptr,
                                  (HYPRE_Matrix) parcsr_A,
                                  (HYPRE_Vector) parvec_b,
                                  (HYPRE_Vector) parvec_x );

      auto t_end = getTimeStamp();

      printf("Hypre setup time (sec): %f\n", getElapsedTime( t_start, t_end ));
   }

   int solve ( const std::vector<value_type>& b, std::vector<value_type>& x_out )
   {
      using namespace hypre_details;
      using HighOrderFEM::getTimeStamp;
      using HighOrderFEM::getElapsedTime;

      auto t_start = getTimeStamp();

      HYPRE_IJVectorSetValues( ij_b, num_rows, NULL, const_cast<value_type*>( b.data() ) );

      double normb = vector_2norm( parvec_b );
      if (mpi_id == 0)
         printf("normb = %e\n", normb);

      if (not(normb > 0.0)) normb = 1.0;

      HYPRE_ParVectorSetConstantValues( parvec_x, 0.0 );
      double init_normr = residual_2norm();
      if (mpi_id == 0)
         printf("Initial residual = %e %e\n", init_normr, init_normr / normb);

      this->solver->solve_func_ptr( this->solver->solver_ptr,
                                   (HYPRE_Matrix) parcsr_A,
                                   (HYPRE_Vector) parvec_b,
                                   (HYPRE_Vector) parvec_x );

      double normr = residual_2norm();
      if (mpi_id == 0)
         printf("Final actual residual = %e %e\n", normr, normr / normb );

      HYPRE_ParVectorGetValues( parvec_x, num_rows, NULL, x_out.data() );

      auto t_end = getTimeStamp();

      printf("Hypre solve time: %f %d\n", getElapsedTime( t_start, t_end ), solver->getNumIterations() );

      return 1;
   }
};

} // LinearSolver
} // HO

#endif
