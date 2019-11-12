import numpy as np
import scipy.sparse as sp
import ctypes as ct
import ctypes.util
import sys, platform
import time
import enum

def timestamp():
    return time.time()

def np_pointer( arr ):

    assert arr.dtype == np.double or arr.dtype == np.float32 or arr.dtype == np.intc

    dtype = None

    if arr.dtype == np.double:
        dtype = ct.c_double
    elif arr.dtype == np.float32:
        dtype = ct.c_float
    elif arr.dtype == np.intc:
        dtype = ct.c_int

    return arr.ctypes.data_as(ct.POINTER(dtype))

def load_lib( lib_name ):
    try:
        lib = ct.CDLL(lib_name)
    except OSError:
        print("Unable to load the system C library: ", lib_name)
        #sys.exit()
        return None

    return lib

shlib = load_lib( "./libitsol_impl.so" )

class PreconTags(enum.IntEnum):
   ILUK = 1
   ARMS = 2
   ILUT = 3
   BILUK = 4
   DEFAULT = ILUK

class itsol:

    def __init__( self, A, precon = None ):

        assert shlib is not None

        self.P = None
        self.solve_func = None
        self.destroy_func = None

        if precon is None:
           precon = 'iluk[1]'

        A_csr = None
        if isinstance(A, sp.csr_matrix):
           A_csr = A
        else:
           A_csr = A.tocsr()

        build_func = shlib.itsol_create
        #print("found build function: ", build_func)

        precon_tag = PreconTags.DEFAULT
        if 'biluk' in precon:
           precon_tag = PreconTags.BILUK
        elif 'iluk' in precon:
           precon_tag = PreconTags.ILUK
        elif 'ilut' in precon:
           precon_tag = PreconTags.ILUT

        lfil = 50
        droptol = 0.0001

        if precon_tag == PreconTags.ILUK or precon_tag == PreconTags.BILUK: # ILUK[fill]

           lfil = 1
           first = precon.find('[')
           last  = precon.find(']')
           if first != -1 and last != -1:
              s = precon[first+1:last]
              if len(s) > 0:
                 lfil = int(s)

           print("iluk({})".format(lfil))

        elif precon_tag == PreconTags.ILUT: # ILUT[drop,fill]

           first = precon.find('[')
           last  = precon.find(']')
           if first != -1 and last != -1:
              opts = precon[first+1:last]
              if len(opts) > 0:
                 opts = opts.split(',')
                 for i,opt in enumerate(opts):
                    opt0 = opt.split('=')
                    if len(opt0) > 1:
                       [key,val] = opt0[:2]
                       if key == 'fill':
                          lfil = int(val)
                       elif key == 'drop':
                          droptol = float(val)
                    elif i == 0:
                       droptol = float(opt0[0])
                    elif i == 1:
                       lfil = int(opt0[0])

           print("ilut({},{})".format(droptol,lfil))

        build_func.argtypes = [ ct.c_int, ct.c_int, ct.c_double, # precon choice, fill, droptol
                                ct.c_int, ct.c_int, # rows, nnz
                                ct.POINTER(ct.c_int), # A_rowptr(rows+1)
                                ct.POINTER(ct.c_int), # A_colidx(nnz)
                                ct.POINTER(ct.c_double), # A_values(nnz)
                                ct.POINTER(ct.c_void_p) ] # &P

        destroy_func = shlib.itsol_destroy
        #print("found destroy function: ", destroy_func)

        destroy_func.argtypes = [ ct.c_void_p ]

        solve_func = shlib.itsol_solve
        #print("found solve function: ", solve_func)

        solve_func.argtypes = [ ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_void_p ]

        A_values = A.data
        A_rowptr = np.copy( A.indptr )
        A_colidx = np.copy( A.indices )

        A_rowptr[:] += 1
        A_colidx[:] += 1

        mrows = A.shape[0]
        nnz = A.nnz

        print("precon= ", precon, precon_tag)

        P = ct.c_void_p()

        build_func( precon_tag, lfil, droptol,
                    mrows, nnz,
                    np_pointer(A_rowptr), np_pointer(A_colidx), np_pointer(A_values),
                    ct.byref(P) )

        self.P = P
        self.solve_func = solve_func
        self.destroy_func = destroy_func

    #def __del__(self):
    #    if self.destroy_func is not None and self.P is not None:
    #        self.destroy_func( self.P )

    def __call__( self, b, x = None ):
        return self.solve(b,x)

    # Solve Mx = b
    def solve( self, b, x = None ):

        if x is None:
            x = np.empty_like(b)

        ierr = self.solve_func( np_pointer( b ),
                                np_pointer( x ),
                                self.P )

        return x
