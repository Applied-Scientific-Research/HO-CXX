import numpy as np
import ctypes as ct
import ctypes.util
import sys, platform

def ct_pointer( arr ):

    assert arr.dtype == np.double or arr.dtype == np.intc

    dtype = None

    if arr.dtype == np.double:
        dtype = ct.c_double
    elif arr.dtype == np.intc:
        dtype = ct.c_int

    return arr.ctypes.data_as(ct.POINTER(dtype))

class sp_solver:

    """ Python wrapper for the C shared library lib"""
 
    def __init__( self, L, U, pr, pc ):

        lib_name = "libsp_dgetrs.so"
        try:
            lib = ct.CDLL(lib_name)
        except OSError:
            print("Unable to load the system C library")
            sys.exit()

        c_double_p = ct.POINTER(ct.c_double)
        c_int_p = ct.POINTER(ct.c_int)

        self.sp_dgetrs_wrapper = lib.sp_dgetrs
        self.sp_dgetrs_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int,
                                              ct.c_int, ct.c_int,
                                              c_int_p, c_int_p, c_double_p,
                                              c_double_p ]
        #self.sp_dgetrs_wrapper.restype = [ ct.c_int ]

        self.L = L
        self.U = U
        self.pr = pr
        self.pc = pc

        assert self.L.data.dtype == np.double
        assert self.L.indices.dtype == np.intc and L.indptr.dtype == np.intc
        assert self.U.data.dtype == np.double
        assert self.U.indices.dtype == np.intc and U.indptr.dtype == np.intc

        print("Loaded sp_dgetr library {}", lib.sp_dgetrs)

    def __call__( self, x, y = None ):
        return self.solve(x,y)

    def solve( self, x, y = None ):

        assert x.dtype == np.double
        assert self.sp_dgetrs_wrapper is not None

        mrows = self.L.shape[0]
        ncols = self.L.shape[1]
        l_nnz = self.L.nnz
        u_nnz = self.U.nnz

        if y is None:
            y = np.empty_like(x)

        # Permute the rhs vector.
        px = np.empty_like(x)
        px[ self.pr[:] ] = x[:]

        # (Pr)Ly = (Pr)x
        isUpper = 0
        isUnit  = 1

        info = self.sp_dgetrs_wrapper (
                         mrows, ncols, l_nnz, isUpper, isUnit,
                         ct_pointer( self.L.indptr ), ct_pointer( self.L.indices ), ct_pointer( self.L.data ),
                         ct_pointer( px ) )

        # Ux' = y
        isUpper = 1
        isUnit  = 0
        info = self.sp_dgetrs_wrapper (
                         mrows, ncols, u_nnz, isUpper, isUnit,
                         ct_pointer( self.U.indptr ), ct_pointer( self.U.indices ), ct_pointer( self.U.data ),
                         ct_pointer( px ) )

        # x = (Pc)x'
        y = px[ self.pc[:] ]

        return y
