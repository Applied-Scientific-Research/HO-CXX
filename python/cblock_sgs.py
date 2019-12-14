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

class cblock_sgs:

    """ Python wrapper for the C shared library lib"""
 
    def __init__( self, niters = 1 ):

        lib_name = "./libcblock_sgs_impl.so"
        # Find the library and load it
        #lib_path = ct.util.find_library("cblock_sgs_impl")
        #if not lib_path:
        #    print("Unable to find the specified library.")
        #    sys.exit()
        
        try:
            lib = ct.CDLL(lib_name)
        except OSError:
            print("Unable to load the system C library")
            sys.exit()

        c_double_p = ct.POINTER(ct.c_double)
        c_int_p = ct.POINTER(ct.c_int)

        self.cblock_sgs_impl_wrapper = lib.cblock_sgs_impl
        self.cblock_sgs_impl_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int,
                                             c_int_p, c_int_p, c_double_p,
                                             c_int_p, c_int_p, c_double_p,
                                             c_double_p, c_double_p, c_double_p ]
        self.cblock_sgs_impl_wrapper.restype = None

        self.niters = niters

        print("Loaded cblock_sgs library {}", lib.cblock_sgs_impl, self.niters)

    def __call__( self, L, Dinv, U, x, y = None ):

        assert L.data.dtype == np.double and x.dtype == np.double
        assert L.indices.dtype == np.intc and L.indptr.dtype == np.intc
        assert self.cblock_sgs_impl_wrapper is not None

        nrows  = L.shape[0]
        bs     = L.blocksize[0]
        nbrows = nrows // bs

        if y is None:
            y = np.empty_like(x)

        self.cblock_sgs_impl_wrapper( nrows, bs, self.niters,
                                 ct_pointer( L.indptr ), ct_pointer( L.indices ), ct_pointer( L.data ),
                                 ct_pointer( U.indptr ), ct_pointer( U.indices ), ct_pointer( U.data ),
                                 ct_pointer( Dinv.data ), ct_pointer( x ), ct_pointer( y ) )

        return y
