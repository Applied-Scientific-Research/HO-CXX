#!/bin/bash

if [[ ! -n $MKLROOT ]] ; then
   . $HOME/mkl/bin/mklvars.sh
fi

export LD_PRELOAD=${MKLROOT}/lib/intel64/libmkl_def.so:${MKLROOT}/lib/intel64/libmkl_avx2.so:${MKLROOT}/lib/intel64/libmkl_core.so:${MKLROOT}/lib/intel64/libmkl_intel_lp64.so:${MKLROOT}/lib/intel64/libmkl_gnu_thread.so:/usr/lib64/libgomp.so.1
