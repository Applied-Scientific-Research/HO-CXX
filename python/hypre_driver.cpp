/*--------------------------------------------------------------------------
 * Test driver for unstructured matrix interface (IJ_matrix interface).
 * Do `driver -help' for usage info.
 * This driver started from the driver for parcsr_linear_solvers, and it
 * works by first building a parcsr matrix as before and then "copying"
 * that matrix row-by-row into the IJMatrix interface. AJC 7/99.
 *--------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "seq_mv.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_parcsr_ls.h"

#include <assert.h>
#include <time.h>

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <string>
#include <iostream>

#include <ctime>
#include <ratio>
#include <chrono>

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef _OPENMP
//#warning 'Using OpenMP wall-clock timer'
typedef double TimerType;

TimerType getTimestamp(void)
{
   return omp_get_wtime();
}
double getElapsedTime( const TimerType& a, const TimerType& b )
{
   return 1000.*(b - a);
}
#else
typedef std::chrono::high_resolution_clock::time_point TimerType;

TimerType getTimestamp(void)
{
   return std::chrono::high_resolution_clock::now();
}

double getElapsedTime( const TimerType& a, const TimerType& b )
{
   using namespace std::chrono;

   std::chrono::duration<double> t_span = duration_cast< duration<double> >(b - a);

   return 1000.*t_span.count();
}
#endif

decltype(hypre_MPI_COMM_WORLD) mpi_comm = hypre_MPI_COMM_WORLD;
int mpi_nprocs = 1, mpi_id = 0;
int ilower = -1, iupper = -1;
int global_num_rows = -1, local_num_rows = -1;

struct SolverParams
{
   int max_iters;
   double rtol, atol;
   int verbosity;
   double normb;
   int k;
   int two_norm;

   SolverParams() : max_iters(100), rtol(1e-8), atol(0), verbosity(0), normb(1), k(16), two_norm(1) {}
};

SolverParams solver_params;

bool compareStrings( const std::string a, const std::string b )
{
   if ( a.length() == b.length() )
   {
      for (size_t i = 0; i < a.length(); ++i)
         if ( tolower(a[i]) != tolower(b[i]))
            return false;

      return true;
   }
   else
      return false;
}

struct SolverType
{
   enum Tags : int { None, GMRES, FlexGMRES, BiCGSTAB, BoomerAMG, DiagScaling };
   static std::vector< std::string > Names;

   static std::string getName( Tags tag )
   {
      return Names[(int)tag];
   }

   static Tags getTag( const int id )
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

   static Tags getTag( const std::string& name )
   {
      for (size_t i = 1; i < Names.size(); ++i)
         if ( compareStrings( Names[i], name ) )
            return getTag(i);

      return getTag(-1);
   }

   typedef HYPRE_Int (*DestructorFcn)( HYPRE_Solver solver );
   typedef HYPRE_Int (*GetNumIterationsFcn)( HYPRE_Solver solver, int* );

   Tags tag;
   HYPRE_Solver object;

   HYPRE_PtrToSolverFcn solve;
   HYPRE_PtrToSolverFcn setup;
   DestructorFcn destructor;
   SolverType *precon;
   GetNumIterationsFcn getNumIterationsFcn;

   int getNumIterations(void)
   {
      if (this->getNumIterationsFcn)
      {
         int n;
         this->getNumIterationsFcn( this->object, &n );
         return n;
      }
      else
         return -1;
   }

   void setNumIterationsFcn( GetNumIterationsFcn fp )
   {
      this->getNumIterationsFcn = fp;
   }

   SolverType(void)
      : tag(None), object(NULL), solve(NULL), setup(NULL), destructor(NULL), precon(NULL), getNumIterationsFcn(NULL)
   {}

   SolverType(Tags tag, HYPRE_Solver obj,
                              HYPRE_PtrToSolverFcn solve,
                              HYPRE_PtrToSolverFcn setup,
                              DestructorFcn destructor,
              SolverType* precon = NULL)
      : tag(tag), object(obj), solve(solve), setup(setup), destructor(destructor), precon(precon), getNumIterationsFcn(NULL)
   {}

   ~SolverType()
   {
      std::cout << "Destroying " << this << " with tag: " << int(this->tag) << " name: " << this->name() << "\n";
      if ( this->object ) {
         this->destructor( this->object );
         this->object = NULL;
      }
   }

   std::string name(void) const { return getName( this->tag ); }
};

std::vector< std::string > SolverType::Names = { "none", "GMRES", "FlexGMRES", "BiCGSTAB", "BoomerAMG", "DiagScaling" };

SolverType create_BoomerAMG( int argc, char* argv[], const bool as_precon = true)
{
   HYPRE_Solver solver;

   int coarsen_type = 10; // hmis
   int coarse_threshold = 500; // this is small enough.
   int npre = 1, npost = 1;
   int relax_type[] = { 0, 13, 14, 9 }; // all, down, up, coarse solvers
   double relax_wt = 1.0;
   int interp_type = 6; // extended classical
   double theta = 0.25;
   double trunc_factor = 0.; // truncation factor for coarsening.
   double agg_trunc_factor = 0.; // aggressive truncation factor for interpolation.

   for (int arg_index = 0; arg_index < argc; arg_index++ )
   {
      if ( strcmp(argv[arg_index], "-coarsen_type") == 0 )
      {
         coarsen_type = atoi( argv[ ++arg_index ] );
      }
      else if ( strcmp(argv[arg_index], "-cljp") == 0 )
         coarsen_type      = 0;
      else if ( strcmp(argv[arg_index], "-cljp1") == 0 )
         coarsen_type      = 7;
      else if ( strcmp(argv[arg_index], "-pmis") == 0 )
         coarsen_type      = 8;
      else if ( strcmp(argv[arg_index], "-pmis1") == 0 )
         coarsen_type      = 9;
      else if ( strcmp(argv[arg_index], "-hmis") == 0 )
         coarsen_type      = 10;
      else if ( strcmp(argv[arg_index], "-ruge") == 0 )
         coarsen_type      = 1;
      else if ( strcmp(argv[arg_index], "-ruge1p") == 0 )
         coarsen_type      = 11;
      else if ( strcmp(argv[arg_index], "-ruge2b") == 0 )
         coarsen_type      = 2;
      else if ( strcmp(argv[arg_index], "-ruge3") == 0 )
         coarsen_type      = 3;
      else if ( strcmp(argv[arg_index], "-ruge3c") == 0 )
         coarsen_type      = 4;
      else if ( strcmp(argv[arg_index], "-rugerlx") == 0 )
         coarsen_type      = 5;
      else if ( strcmp(argv[arg_index], "-falgout") == 0 )
         coarsen_type      = 6;
      else if ( strcmp(argv[arg_index], "-coarse_enough") == 0 )
      {
         coarse_threshold = atoi(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-npre") == 0 )
      {
         npre = atoi(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-npost") == 0 )
      {
         npost = atoi(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-relax_type") == 0 )
      {
         relax_type[0] = atoi(argv[ ++arg_index ]);
         relax_type[1] = relax_type[0];
         relax_type[2] = relax_type[0];
      }
      else if ( strcmp(argv[arg_index], "-relax_wt") == 0 )
      {
         relax_wt = atof(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-interp_type") == 0 )
      {
         interp_type = atoi(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-theta") == 0 )
      {
         theta = atof(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-trunc_factor") == 0 )
      {
         trunc_factor = atof(argv[ ++arg_index ]);
      }
      else if ( strcmp(argv[arg_index], "-agg_trunc_factor") == 0 )
      {
         agg_trunc_factor = atof(argv[ ++arg_index ]);
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

   SolverType st( SolverType::BoomerAMG,
                  solver, 
                  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                  HYPRE_BoomerAMGDestroy );
   st.setNumIterationsFcn( HYPRE_BoomerAMGGetNumIterations );

   return st;
}

SolverType create_DiagScaling(void)
{
   return SolverType( SolverType::DiagScaling,
                      NULL, 
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
                      NULL );
}

SolverType create_precon( std::string name, int argc, char* argv[] )
{
   SolverType::Tags tag = SolverType::getTag( name );
   if ( tag == SolverType::BoomerAMG  ) return create_BoomerAMG( argc, argv, true );
   if ( tag == SolverType::DiagScaling) return create_DiagScaling();

   fprintf(stderr,"Precon not found in create_precon %s\n", name.c_str());
   exit(1);
}

SolverType create_GMRES( int argc, char* argv[], SolverType *precon )
{
   HYPRE_Solver solver;

   HYPRE_ParCSRGMRESCreate(hypre_MPI_COMM_WORLD, &solver);

   HYPRE_GMRESSetKDim(solver, solver_params.k);
   HYPRE_GMRESSetMaxIter(solver, solver_params.max_iters);
   HYPRE_GMRESSetTol(solver, solver_params.rtol);
   HYPRE_GMRESSetAbsoluteTol(solver, solver_params.atol);
   HYPRE_GMRESSetLogging(solver, 1 );
   HYPRE_GMRESSetPrintLevel(solver, solver_params.verbosity );

   if (precon != NULL)
      HYPRE_GMRESSetPrecond( solver,
                             (HYPRE_PtrToSolverFcn) precon->solve,
                             (HYPRE_PtrToSolverFcn) precon->setup,
                             (HYPRE_Solver) precon->object );

   SolverType st( SolverType::GMRES,
                      solver,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRGMRESSolve,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRGMRESSetup,
                      HYPRE_ParCSRGMRESDestroy,
                      precon );

   st.setNumIterationsFcn( HYPRE_GMRESGetNumIterations );

   return st;
}

SolverType create_FlexGMRES( int argc, char* argv[], SolverType *precon )
{
   HYPRE_Solver solver;

   HYPRE_ParCSRFlexGMRESCreate(hypre_MPI_COMM_WORLD, &solver);

   HYPRE_FlexGMRESSetKDim(solver, solver_params.k);
   HYPRE_FlexGMRESSetMaxIter(solver, solver_params.max_iters);
   HYPRE_FlexGMRESSetTol(solver, solver_params.rtol);
   HYPRE_FlexGMRESSetAbsoluteTol(solver, solver_params.atol);
   HYPRE_FlexGMRESSetLogging(solver, 1 );
   HYPRE_FlexGMRESSetPrintLevel(solver, solver_params.verbosity );

   if (precon != NULL)
      HYPRE_FlexGMRESSetPrecond( solver,
                             (HYPRE_PtrToSolverFcn) precon->solve,
                             (HYPRE_PtrToSolverFcn) precon->setup,
                             (HYPRE_Solver) precon->object );

   SolverType st( SolverType::FlexGMRES,
                      solver,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRFlexGMRESSolve,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRFlexGMRESSetup,
                      HYPRE_ParCSRFlexGMRESDestroy,
                      precon );

   st.setNumIterationsFcn( HYPRE_FlexGMRESGetNumIterations );

   return st;
}

SolverType create_BiCGSTAB( int argc, char* argv[], SolverType *precon )
{
   HYPRE_Solver solver;

   HYPRE_ParCSRBiCGSTABCreate(hypre_MPI_COMM_WORLD, &solver);

   HYPRE_BiCGSTABSetMaxIter(solver, solver_params.max_iters);
   HYPRE_BiCGSTABSetTol(solver, solver_params.rtol);
   HYPRE_BiCGSTABSetAbsoluteTol(solver, solver_params.atol);
   HYPRE_BiCGSTABSetLogging(solver, 1 );
   HYPRE_BiCGSTABSetPrintLevel(solver, solver_params.verbosity );

   if (precon != NULL)
      HYPRE_BiCGSTABSetPrecond( solver,
                             (HYPRE_PtrToSolverFcn) precon->solve,
                             (HYPRE_PtrToSolverFcn) precon->setup,
                             (HYPRE_Solver) precon->object );

   SolverType st( SolverType::BiCGSTAB,
                      solver,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRBiCGSTABSolve,
                      (HYPRE_PtrToSolverFcn) HYPRE_ParCSRBiCGSTABSetup,
                      HYPRE_ParCSRBiCGSTABDestroy,
                      precon );

   st.setNumIterationsFcn( HYPRE_BiCGSTABGetNumIterations );

   return st;
}

SolverType create_solver( std::string name, SolverType* precon, int argc, char* argv[] )
{
   SolverType::Tags tag = SolverType::getTag( name );
   if ( tag == SolverType::GMRES ) return create_GMRES( argc, argv, precon );
   if ( tag == SolverType::FlexGMRES ) return create_FlexGMRES( argc, argv, precon );
   if ( tag == SolverType::BiCGSTAB ) return create_BiCGSTAB( argc, argv, precon );

   fprintf(stderr,"Solver not found in create_solver %s\n", name.c_str());
   exit(1);
}

void print_usage(const char* prg)
{
   printf("\n");
   printf("Usage: %s [<options>]\n", prg);
   printf("\n");
   printf("  -A <str>            : matrix read from a single file (CSR format)\n");
   printf("  -b <str>            : read rhs from a single file (CSR format)\n");
   printf("  -precon <str>       : [ Diagonal, AMG ]\n");
   printf("  -solver <str>       : [ GMRES, FlexGMRES, BiCGSTAB ]\n");
   printf("  -rtol  <float>      : set solver relative convergence tolerance\n");
   printf("  -atol  <float>      : set solver absolute convergence tolerance\n");
   printf("  -max_iters  <int>   : set max iterations\n");
   printf("  -k          <int>   : set krylov space size for GMRES and family\n");
}

HYPRE_IJMatrix ReadCSRMatrix( const char *file_name )
{
   FILE* fp = fopen(file_name, "rb");
   if (fp == NULL)
   {
      fprintf(stderr,"Error opening file %s\n", file_name);
      exit(1);
   }

   int num_rows;
   fread( &num_rows, 1, sizeof(int), fp );

   std::vector<int> rowptr( num_rows+1 );

   fread( rowptr.data(), num_rows+1, sizeof(int), fp );

   int num_nonzeros = rowptr[ num_rows ];
   printf("%d %d\n", num_rows, num_nonzeros);

   std::vector<int> colidx( num_nonzeros );
   std::vector<double> values( num_nonzeros );

   fread( colidx.data(), num_nonzeros, sizeof(int), fp );
   fread( values.data(), num_nonzeros, sizeof(double), fp );

   fclose(fp);

   ilower = 0;
   iupper = num_rows - 1;
   local_num_rows = iupper - ilower + 1;
   global_num_rows = num_rows;

   HYPRE_IJMatrix A;

   HYPRE_IJMatrixCreate(hypre_MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);

   for (int i = ilower; i <= iupper; ++i)
   {
      int k0 = rowptr[i];
      int nnz = rowptr[i+1] - k0;
      double *data_i = values.data() + k0;
      int *cols_i = colidx.data() + k0;

      HYPRE_IJMatrixSetValues( A, 1, &nnz, &i, cols_i, data_i );
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);

   return A;
}

HYPRE_IJVector ReadVector( const char *file_name )
{
   FILE *fp = fopen(file_name, "r");
   if (fp == NULL)
   {
      fprintf(stderr,"Error opening file %s\n", file_name);
      exit(1);
   }

   int size;
   fread( &size, 1, sizeof(int), fp );

   std::vector<double> values(size);

   fread( values.data(), size, sizeof(double), fp );

   fclose(fp);

   if ( size != global_num_rows )
   {
      fprintf(stderr,"size != global_num_rows %d %d\n", size, global_num_rows);
      exit(1);
   }

   HYPRE_IJVector x;
   HYPRE_IJVectorCreate(hypre_MPI_COMM_WORLD, ilower, iupper, &x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);

   std::vector<int> indices(local_num_rows);
   for (int i = 0; i < local_num_rows; ++i)
      indices[i] = ilower + i;

   HYPRE_IJVectorSetValues( x, local_num_rows, indices.data(), values.data() + ilower );

   HYPRE_IJVectorAssemble(x);

   return x;
}

int main( int argc, char* argv[] )
{
   /* Initialize MPI */
   hypre_MPI_Init( &argc, &argv );

   hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &mpi_nprocs );
   hypre_MPI_Comm_rank(hypre_MPI_COMM_WORLD, &mpi_id );

   /* Initialize Hypre */
   HYPRE_Init(argc, argv);

   //nvtxDomainHandle_t domain = nvtxDomainCreateA("Domain_A");
   /*
      hypre_InitMemoryDebug(myid);
   */

   HYPRE_IJMatrix      ij_A = NULL;
   HYPRE_IJVector      ij_x = NULL;
   HYPRE_IJVector      ij_b = NULL;
   HYPRE_IJVector      ij_r = NULL;

   auto finalize = [&](void)
      {
         if (ij_A) HYPRE_IJMatrixDestroy( ij_A );
         if (ij_x) HYPRE_IJVectorDestroy( ij_x );
         if (ij_b) HYPRE_IJVectorDestroy( ij_b );
         if (ij_r) HYPRE_IJVectorDestroy( ij_r );

         MPI_Finalize();
      };

   if (mpi_nprocs > 1 )
   {
      fprintf(stderr,"Only 1 MPI procs is support at this point\n");
      finalize();
      return 0;
   }

   /*-----------------------------------------------------------
    * Initialize some stuff
    *-----------------------------------------------------------*/

   char *A_file = NULL, *b_file = NULL;

   std::string solver_name("gmres");
   std::string precon_name("diagscaling");

   bool omp_flag = false;
   int nsolves = 1;

   /*-----------------------------------------------------------
    * Parse command line
    *-----------------------------------------------------------*/

   for (int arg_index = 1; arg_index < argc; arg_index++ )
   {
      if ( strcmp(argv[arg_index], "-A") == 0 )
      {
         A_file = argv[++arg_index];
      }
      else if ( strcmp(argv[arg_index], "-b") == 0 )
      {
         b_file = argv[++arg_index];
      }
      else if ( strcmp(argv[arg_index], "-omp") == 0 )
      {
         omp_flag = true;
      }
      else if ( strcmp(argv[arg_index], "-solver") == 0 )
      {
         solver_name = argv[++arg_index];
      }
      else if ( strcmp(argv[arg_index], "-precon") == 0 )
      {
         precon_name = argv[++arg_index];
      }
      else if ( strcmp(argv[arg_index], "-max_iters") == 0 )
      {
         solver_params.max_iters = atoi(argv[++arg_index]);
      }
      else if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         print_usage( argv[0] );
         finalize();
         return 0;
      }
      else if ( strcmp(argv[arg_index], "-atol") == 0 )
      {
         solver_params.atol  = atof( argv[++arg_index] );
      }
      else if ( strcmp(argv[arg_index], "-rtol") == 0 )
      {
         solver_params.rtol  = atof( argv[++arg_index] );
      }
      else if ( strcmp(argv[arg_index], "-v") == 0 )
      {
         solver_params.verbosity ++;
      }
      else if ( strcmp(argv[arg_index], "-k") == 0 )
      {
         solver_params.k = atoi(argv[++arg_index]);
      }
      else if ( strcmp(argv[arg_index], "-nsolves") == 0 )
      {
         nsolves = atoi(argv[++arg_index]);
      }
   }

   if ( A_file )
      ij_A = ReadCSRMatrix( A_file );
   else
   {
      fprintf(stderr,"Input matrix not specified.\n");
      finalize();
      return 1;
   }

   if ( b_file )
      ij_b = ReadVector( b_file );
   else
   {
      fprintf(stderr,"Input Rhs not specified.\n");
      finalize();
      return 1;
   }

   // Create solution and residual vectors.
   {
      HYPRE_IJVectorCreate( hypre_MPI_COMM_WORLD, ilower, iupper, &ij_x );
      HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_x);

      std::vector<double> values( local_num_rows );
      for (int i = 0; i < local_num_rows; i++)
          values[i] = 0.;

      HYPRE_IJVectorSetValues( ij_x, local_num_rows, NULL, values.data() );

      HYPRE_IJVectorAssemble( ij_x );

      HYPRE_IJVectorCreate( hypre_MPI_COMM_WORLD, ilower, iupper, &ij_r );
      HYPRE_IJVectorSetObjectType(ij_r, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(ij_r);

      HYPRE_IJVectorSetValues( ij_r, local_num_rows, NULL, values.data() );

      HYPRE_IJVectorAssemble( ij_r );
   }

   if (omp_flag)
      HYPRE_IJMatrixSetOMPFlag( ij_A, 1 );

   HYPRE_ParCSRMatrix  parcsr_A = NULL;
   HYPRE_ParVector     parvec_b = NULL;
   HYPRE_ParVector     parvec_x = NULL;
   HYPRE_ParVector     parvec_r = NULL;

   /* Get the parcsr matrix/vector objects to use */
   HYPRE_IJMatrixGetObject( ij_A, (void**) &parcsr_A);
   HYPRE_IJVectorGetObject( ij_x, (void**) &parvec_x );
   HYPRE_IJVectorGetObject( ij_b, (void**) &parvec_b );
   HYPRE_IJVectorGetObject( ij_r, (void**) &parvec_r );

   {
      HYPRE_ParVectorInnerProd( parvec_b, parvec_b, &solver_params.normb );
      solver_params.normb = sqrt(solver_params.normb);

      if (mpi_id == 0)
         printf("norm(b)= %e\n", solver_params.normb);
   }

   auto compute_normr = [&](void)
   {
      hypre_ParVectorCopy( parvec_b, parvec_r );
      hypre_ParCSRMatrixMatvec( -1.0, parcsr_A, parvec_x, 1.0, parvec_r );

      double normr;
      HYPRE_ParVectorInnerProd( parvec_r, parvec_r, &normr );

      normr = sqrt(normr);

      return normr;
   };

   HYPRE_ParVectorSetConstantValues( parvec_x, 0.0 );
   double normr = compute_normr();
   if (mpi_id == 0)
      printf("Initial residual = %e %e\n", normr, normr / ( solver_params.normb > 0.0 ? solver_params.normb : 1.0 ));

   HYPRE_ParVectorSetConstantValues( parvec_x, 1.0 );
   normr = compute_normr();
   if (mpi_id == 0)
      printf("Residual with x=1 %e\n", normr);

   HYPRE_ParVectorSetConstantValues( parvec_x, 0.0 );

   std::cout << "Solver: " << solver_name << ' ' << SolverType::getTag(solver_name) << '\n';
   std::cout << "Precon: " << precon_name << ' ' << SolverType::getTag(precon_name) << '\n';

   auto precon = create_precon( precon_name, argc, argv );
   auto solver = create_solver( solver_name, &precon, argc, argv );

   {
      auto t_start = getTimestamp();

      solver.setup( solver.object,
                    (HYPRE_Matrix) parcsr_A,
                    (HYPRE_Vector) parvec_b,
                    (HYPRE_Vector) parvec_x );

      auto t_end = getTimestamp();

      printf("Setup time (ms): %f\n", getElapsedTime( t_start, t_end ));
   }

   for (int iter = 0; iter < nsolves; ++iter)
   {
      HYPRE_ParVectorSetConstantValues( parvec_x, 0.0 );

      auto t_start = getTimestamp();

      solver.solve( solver.object,
                    (HYPRE_Matrix) parcsr_A,
                    (HYPRE_Vector) parvec_b,
                    (HYPRE_Vector) parvec_x );

      auto t_end = getTimestamp();

      printf("Solve time (ms): %f %d\n", getElapsedTime( t_start, t_end ), solver.getNumIterations() );

      normr = compute_normr();
      if (mpi_id == 0)
         printf("Final actual residual = %e %e\n", normr, normr / ( solver_params.normb > 0.0 ? solver_params.normb : 1.0 ));
   }

   finalize();

   return (0);
}
