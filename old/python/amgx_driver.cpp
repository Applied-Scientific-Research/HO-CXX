/*--------------------------------------------------------------------------
 * Test driver for unstructured matrix interface (IJ_matrix interface).
 * Do `driver -help' for usage info.
 * This driver started from the driver for parcsr_linear_solvers, and it
 * works by first building a parcsr matrix as before and then "copying"
 * that matrix row-by-row into the IJMatrix interface. AJC 7/99.
 *--------------------------------------------------------------------------*/
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <ctime>

#include <vector>
#include <string>
#include <iostream>

#include "amgx_c.h"

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

#include <ratio>
#include <chrono>

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

/* print callback (could be customized) */
void print_callback(const char *msg, int length)
{
    printf("%s", msg);
}

int main( int argc, char* argv[] )
{
    auto finalize = [&](int exit_code)
       {
          AMGX_SAFE_CALL( AMGX_finalize_plugins() );
          AMGX_SAFE_CALL( AMGX_finalize() );
          if (exit_code) exit( exit_code );
       };

    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());

    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());

    /* get api and build info */
    {
        int major, minor;
        AMGX_get_api_version(&major, &minor);
        printf("amgx api version: %d.%d\n", major, minor);

        char *ver, *date, *time;
        AMGX_get_build_info_strings(&ver, &date, &time);
        printf("amgx build version: %s\nBuild date and time: %s %s\n", ver, date, time);
    }

    std::vector< std::string > args;
    for (int i = 0; i < argc; ++i)
       args.push_back( std::string( argv[i] ) );

    AMGX_Mode mode = AMGX_mode_hDDI;

    int num_solves = 1;
    int max_iters = 100;
    int verbosity = 0;
    int bs = 1;
    int kdim = 16;
    double atol = 1e-8, rtol = 0;

    std::string A_file, b_file;
    std::string solver("fgmres");
    std::string precon("none");
    std::string config("");
    std::string config_file;

    /*-----------------------------------------------------------
     * Parse command line
     *-----------------------------------------------------------*/

    for (int i = 1; i < args.size(); i++)
    {
      if ( args[i] == "-A" )
         A_file = args[++i];
      else if ( args[i] == "-b" )
         b_file = args[++i];
      else if ( args[i] == "-solver" )
         solver = args[++i];
      else if ( args[i] == "-precon" )
         precon = args[++i];
      else if ( args[i] == "-max_iters" )
         max_iters = atoi( args[++i].c_str() );
      else if ( args[i] == "-atol" )
         atol  = atof( args[++i].c_str() );
      else if ( args[i] == "-rtol" )
         rtol  = atof( args[++i].c_str() );
      else if ( args[i] == "-v" )
         verbosity ++;
      else if ( args[i] == "-k" )
         kdim = atoi( args[++i].c_str() );
      else if ( args[i] == "-bs" )
         bs = atoi( args[++i].c_str() );
      else if ( args[i] == "-num_solves" )
         num_solves = atoi( args[++i].c_str() );
      else if ( args[i] == "-config" )
         config = args[++i];
      else if ( args[i] == "-config_file" )
         config_file = args[++i];
      else if ( args[i] == "-mode" )
      {
         std::string smode = args[++i];
         if ( smode == "hDDI" )
            mode = AMGX_mode_hDDI;
         else if ( smode == "dDDI" )
            mode = AMGX_mode_dDDI;
         else
         {
            fprintf(stderr,"%s mode not ready\n", smode.c_str());
            finalize(1);
         }
      }
   }

   int num_rows = 0, num_cols = 0, num_nonzeros = 0;

   std::vector<int> A_rowptr, A_colidx;
   std::vector<double> A_values;
   std::vector<double> b_values;

   if ( A_file.size() )
   {
      FILE* fp = fopen( A_file.c_str(), "rb");
      if (fp == NULL)
      {
         fprintf(stderr,"Error opening file %s\n", A_file.c_str() );
         finalize(1);
      }

      fread( &num_rows, 1, sizeof(int), fp );

      A_rowptr.resize( num_rows+1 );

      fread( A_rowptr.data(), num_rows+1, sizeof(int), fp );

      num_nonzeros = A_rowptr[ num_rows ];
      printf("%d %d\n", num_rows, num_nonzeros);

      A_colidx.resize( num_nonzeros );
      A_values.resize( num_nonzeros );

      fread( A_colidx.data(), num_nonzeros, sizeof(int), fp );
      fread( A_values.data(), num_nonzeros, sizeof(double), fp );

      fclose(fp);
   }
   else
   {
      fprintf(stderr,"Input matrix not specified.\n");
      finalize(1);
   }

   if ( b_file.size() )
   {
      FILE *fp = fopen( b_file.c_str(), "r");
      if (fp == NULL)
      {
         fprintf(stderr,"Error opening file %s\n", b_file.c_str() );
         finalize(1);
      }

      int size;
      fread( &size, 1, sizeof(int), fp );
      if ( size != num_rows )
      {
         fprintf(stderr,"b size != num_rows %d %d\n", size, num_rows);
         finalize(1);
      }

      b_values.resize(size);

      fread( b_values.data(), size, sizeof(double), fp );

      fclose(fp);
   }
   else
   {
      fprintf(stderr,"Input Rhs not specified.\n");
      finalize(1);
   }

   AMGX_config_handle cfg = NULL;
   if ( config_file.empty() )
      AMGX_SAFE_CALL( AMGX_config_create( &cfg, config.c_str() ) )
   else
      AMGX_SAFE_CALL( AMGX_config_create_from_file_and_string(&cfg, config_file.c_str(), config.c_str() ) )

   AMGX_resources_handle rsrc = NULL;
   AMGX_SAFE_CALL( AMGX_resources_create_simple( &rsrc, cfg ) );

   AMGX_matrix_handle A_handle = NULL;
   AMGX_vector_handle b_handle = NULL, x_handle = NULL;

   AMGX_matrix_create(&A_handle, rsrc, mode);
   AMGX_vector_create(&x_handle, rsrc, mode);
   AMGX_vector_create(&b_handle, rsrc, mode);

   AMGX_matrix_upload_all( A_handle, num_rows, num_nonzeros, bs, bs,
                           A_rowptr.data(), A_colidx.data(), A_values.data(), NULL );
   AMGX_vector_upload( b_handle, num_rows, bs, b_values.data() );
   AMGX_vector_set_zero( x_handle, num_rows, bs );

   {
      int n, bx, by;
      AMGX_matrix_get_size( A_handle, &n, &bx, &by );
      printf("%d %d %d %d %d\n", num_rows, bs, n, bx, by);
      AMGX_vector_get_size( x_handle, &n, &bx );
      printf("%d %d %d %d\n", num_rows, bs, n, bx);
   }

   AMGX_solver_handle solver_handle = NULL;
   AMGX_solver_create( &solver_handle, rsrc, mode, cfg );

   /* solver setup */
   AMGX_solver_setup(solver_handle, A_handle);

   auto get_resid = [&](const char* msg)
   {
      std::vector<double> norm(bs);
      AMGX_solver_calculate_residual_norm( solver_handle, A_handle, b_handle, x_handle, norm.data() );
      printf("%s: ", msg);
      for (int i = 0; i < norm.size(); ++i)
         printf("%e%c", norm[i], (i == norm.size()-1 ? '\n' : ' ') );
      printf("\n");
   };

   get_resid("initial");

   /* solver solve */
   AMGX_solver_solve(solver_handle, b_handle, x_handle);

   AMGX_SOLVE_STATUS status;
   AMGX_solver_get_status(solver_handle, &status);
   int nit;
   double res;
   AMGX_solver_get_iterations_number(solver_handle, &nit);
   printf("status: %d, %d\n", nit, (int)status);

   get_resid("final");

   AMGX_solver_destroy( solver_handle );

   AMGX_vector_destroy( x_handle );
   AMGX_vector_destroy( b_handle );
   AMGX_matrix_destroy( A_handle );

   finalize(0);

   return 0;
}
