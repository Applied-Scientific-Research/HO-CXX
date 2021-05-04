#!/bin/bash
## opts.solver_par.solver     = Magma_PBICGSTAB;
# opts.solver_par.restart    = 16;
# opts.solver_par.maxiter    = 100;
# opts.solver_par.rtol       = 1e-10;
# opts.solver_par.verbose    = 1;
# //opts.precond_par.solver    = Magma_PARILU;
# opts.precond_par.solver    = Magma_ILU;
# //opts.precond_par.solver    = Magma_JACOBI;
# //opts.precond_par.solver    = Magma_PARILU;
# opts.precond_par.levels    = 4;
# opts.precond_par.trisolver = Magma_CUSOLVE;

function run()
{
   echo "$@"
   eval "$@"
}

opts="
--solver PBICGSTAB
--atol 1e-10
--rtol 0
--verbose 100
--precond ILU
  --plevels 4
  --patol 4.0
  --prtol 1e-5
  --psweeps 5
"
MAGMA_DIR=$HOME/high-order/magma
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MAGMA_DIR}/lib
run ./magma_driver ../../HO-Fortran/Test_Package/a.out ../../HO-Fortran/Test_Package/b.out $opts $*
