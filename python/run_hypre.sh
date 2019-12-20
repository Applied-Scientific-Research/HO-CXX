#!/bin/bash

function run()
{
   echo "$@"
   eval "$@"
}

opts="
-max_iters 100
-solver FlexGMRES
-precon boomeramg
-atol 1e-10
-rtol 0.0
-pmis
-interp_type 7
-theta 0.6
-k 16
-relax_type 0
-relax_wt 0.533333
-npost 2
-trunc_factor 0.25
"

run ./hypre_driver -A ../Test_Package/a.out -b ../Test_Package/b.out $opts $*
