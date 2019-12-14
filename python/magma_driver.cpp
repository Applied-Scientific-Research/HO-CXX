// This is a simple standalone example. See README.txt

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "magma_v2.h"
#include "magmasparse.h"


// ------------------------------------------------------------
// This is an example how magma can be integrated into another software.
int main( int argc, char** argv )
{
   int m, n = 1;
   std::vector<int> rowptr, colidx;
   std::vector<double> values, rhs;

   if ( argc < 3 )
   {
      fprintf(stderr,"Must specify the file\n");
      return 1;
   }

   {
      FILE *fp = fopen(argv[1],"r");
      fread( &m, sizeof(int), 1, fp );

      rowptr.resize(m+1);
      fread( rowptr.data(), sizeof(int), m+1, fp );

      int nnz = rowptr[m];
      colidx.resize(nnz);
      values.resize(nnz);

      fread( colidx.data(), sizeof(int), nnz, fp );
      fread( values.data(), sizeof(double), nnz, fp );

      fclose(fp);
   }

   {
      FILE *fp = fopen(argv[2],"r");
      int i;
      fread( &i, sizeof(int), 1, fp );
      if (i != m ) {
         fprintf(stderr,"i!=m %d %d\n", i, m);
         return 1;
      }

      rhs.resize(m);
      fread( rhs.data(), sizeof(double), m, fp );

      fclose(fp);
   }
   std::cout << "Loaded problem: " << m << ", " << rowptr[m] << "\n";

   // Initialize MAGMA and create some LA structures.
   magma_init();
   magma_dopts opts;
   magma_queue_t queue;
   magma_queue_create( 0, &queue );
    
   magma_d_matrix A={Magma_CSR}, dA={Magma_CSR};
   magma_d_matrix b={Magma_CSR}, db={Magma_CSR};
   magma_d_matrix x={Magma_CSR}, dx={Magma_CSR};

   // Pass the system to MAGMA.
   magma_dcsrset( m, m, rowptr.data(), colidx.data(), values.data(), &A, queue );

   //write_d_csrtomtx( A, "a.mtx", queue );

#define SAFECALL(_cmd_) \
{ int _ierr = (_cmd_); \
  if(_ierr != MAGMA_SUCCESS) { \
    fprintf(stderr,"Error in %s\n", #_cmd_); \
    exit(_ierr); \
}}

   //{
   //   magma_d_matrix Ab = {Magma_BCSR};
   //   Ab.blocksize = 9;
   //   SAFECALL( magma_dmconvert( A, &Ab, Magma_CSR, Magma_BCSR, queue ) );
   //   std::cout << "Converted" << std::endl;
   //   SAFECALL( magma_dmfree( &Ab, queue ) );
   //}
    
   magma_dvset( m, 1, rhs.data(), &b, queue );
 
  // // Choose a solver, preconditioner, etc. - see documentation for options.
  // opts.solver_par.solver     = Magma_PBICGSTAB;
  // opts.solver_par.restart    = 16;
  // opts.solver_par.maxiter    = 100;
  // opts.solver_par.rtol       = 1e-10;
  // opts.solver_par.verbose    = 1;
  // //opts.precond_par.solver    = Magma_PARILU;
  // opts.precond_par.solver    = Magma_ILU;
  // //opts.precond_par.solver    = Magma_JACOBI;
  // //opts.precond_par.solver    = Magma_PARILU;
  // opts.precond_par.levels    = 4;
  // opts.precond_par.trisolver = Magma_CUSOLVE;

   int nmats = 0;
   magma_dparse_opts( argc-2, argv+2, &opts, &nmats, queue );

   std::map<int,std::string> solver_name;
   solver_name[Magma_PBICGSTABMERGE] = "Magma_PBICGSTABMERGE";
   solver_name[Magma_PGMRES] = "Magma_PGMRES";

   std::map<int,std::string> precond_name;
   precond_name[Magma_PARILU] = "Magma_PARILU";
   precond_name[Magma_ILU] = "Magma_ILU";
   precond_name[Magma_JACOBI] = "Magma_JACOBI";

   std::cout << "format : " << opts.output_format << "\n";
   std::cout << "solver : " << opts.solver_par.solver << " " << solver_name[opts.solver_par.solver] << "\n";
   std::cout << "  restart: " << opts.solver_par.restart << "\n";
   std::cout << "precond: " << opts.precond_par.solver << " " << precond_name[opts.precond_par.solver] << "\n";
   std::cout << "  levels " << opts.precond_par.levels << "\n";
   std::cout << "  fill   " << opts.precond_par.atol << "\n";
   std::cout << "  dtol   " << opts.precond_par.rtol << "\n";
   std::cout << "  sweeps " << opts.precond_par.sweeps << "\n";

   // Initialize the solver.
   magma_dsolverinfo_init( &opts.solver_par, &opts.precond_par, queue );
    
   // Copy the system to the device
   magma_dmtransfer( A, &dA, Magma_CPU, Magma_DEV, queue );
   magma_dmtransfer( b, &db, Magma_CPU, Magma_DEV, queue );

   // initialize an initial guess for the iteration vector
   magma_dvinit( &dx, Magma_DEV, b.num_rows, b.num_cols, 0.0, queue );

   std::cout << "Initialized device data\n";

   {
      auto t_start = magma_wtime();

      // Generate the preconditioner.
      magma_d_precondsetup( dA, db, &opts.solver_par, &opts.precond_par, queue );

      auto t_stop = magma_wtime();
      printf("Finished precond setup in %g (ms)\n", 1000.*(t_stop-t_start));

      printf("L: %d, %d, %d\n", opts.precond_par.L.num_rows, opts.precond_par.L.num_cols, opts.precond_par.L.nnz);
      printf("U: %d, %d, %d\n", opts.precond_par.U.num_rows, opts.precond_par.U.num_cols, opts.precond_par.U.nnz);

      //if (opts.solver_par.verbose)
      //   magma_dsolverinfo( &opts.solver_par, &opts.precond_par, queue );
   }

   // In case we only wanted to generate a preconditioner, we are done.
   // The preconditioner in the opts.precond_par structure - in this case an ILU.
   // The lower ILU(0) factor is in opts.precond_par.L and
   // the upper ILU(0) factor is in opts.precond_par.U (in this case on the device).

   {
      auto t_start = magma_wtime();

      // If we want to solve the problem, run:
      magma_d_solver( dA, db, &dx, &opts, queue );

      auto t_stop = magma_wtime();
      auto t_solve = 1000.*(t_stop-t_start);
      printf("Finished solver: info:       %d\n", opts.solver_par.info);
      printf("                 iterations: %d\n", opts.solver_par.numiter);
      printf("                 residual:   %.2e %.2e\n", opts.solver_par.init_res,
                                                         opts.solver_par.final_res);
      printf("                 runtime:    %g\n", t_solve);

      if (opts.solver_par.verbose)
         magma_dsolverinfo( &opts.solver_par, &opts.precond_par, queue );
   }

   // Then copy the solution back to the host...
   magma_dmtransfer( dx, &x, Magma_DEV, Magma_CPU, queue );

   std::vector<double> sol(m);
   // and back to the application code
   magma_dvcopy( x, &m, &n, sol.data(), queue );

   // Free the allocated memory...
   magma_dmfree( &dx, queue );
   magma_dmfree( &db, queue ); 
   magma_dmfree( &dA, queue );
   magma_dmfree( &b, queue );  // won't do anything as MAGMA does not own the data. 
   magma_dmfree( &A, queue );  // won't do anything as MAGMA does not own the data. 

   // and finalize MAGMA.
   magma_queue_destroy( queue );
   magma_finalize();
    
   // From here on, the application code may continue with the solution in sol...
   //for (int i = 0; i < std::min(m,20); ++i)
   //   printf("%.4f\n", sol[i]);

   return 0;
}
