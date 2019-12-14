import numpy as np
import scipy.sparse as sp
import ctypes as ct
import ctypes.util
import sys, platform
import time

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

shlib = load_lib( "BILUM/libilum.so" )
shlib = load_lib( "BILUM/libilum_mp.so" )

class bilum:

    def __init__( self, A, options = None ):

        assert shlib is not None

        build_func = shlib.bilum_
        print("found build function: ", build_func)

        build_func.restype = None # void
        build_func.argtypes = [ ct.POINTER(ct.c_int), # n
                                ct.POINTER(ct.c_double), # A(*)
                                ct.POINTER(ct.c_int), # A_colidx(*)
                                ct.POINTER(ct.c_int), # A_rowptr(n+1)
                                ct.POINTER(ct.c_double), # LU(*)
                                ct.POINTER(ct.c_int), # LU_colidx(*)
                                ct.POINTER(ct.c_int), # LU_rowptr(*)
                                ct.POINTER(ct.c_int), # iord(n)
                                ct.POINTER(ct.c_int), # riord(*)
                                ct.POINTER(ct.c_int), # nlast(*)
                                ct.POINTER(ct.c_int), # iwk(*)
                                ct.POINTER(ct.c_int), # iwk1(n)
                                ct.POINTER(ct.c_int), # iwork(n)
                                ct.POINTER(ct.c_int), # ipar(10)
                                ct.POINTER(ct.c_double), # droptol
                                ct.POINTER(ct.c_double), # A1(*)
                                ct.POINTER(ct.c_int), # A1_colidx(*)
                                ct.POINTER(ct.c_int), # A1_rowptr(n+1)
                               ]

        solve_func = shlib.bsolve_
        print("found solve function: ", solve_func)

        solve_func.restype = None # void
        solve_func.argtypes = [ ct.POINTER(ct.c_int), # n
                                ct.POINTER(ct.c_double), # b(n), input
                                ct.POINTER(ct.c_double), # x(n), output
                                ct.POINTER(ct.c_double), # rwk(n)
                                ct.POINTER(ct.c_double), # LU(*)
                                ct.POINTER(ct.c_int), # LU_colidx(*)
                                ct.POINTER(ct.c_int), # LU_rowptr(*)
                                ct.POINTER(ct.c_int), # perm(*)
                                ct.POINTER(ct.c_int), # nlast(nlvls)
                                ct.POINTER(ct.c_int), # nlvls
                                ct.POINTER(ct.c_int), # iwk(n)
                               ]

        # Convert to base-1
        A_data = A.data
        A_rowptr = np.copy(A.indptr)
        A_colidx = np.copy(A.indices)

        A_rowptr[:] += 1
        A_colidx[:] += 1

        mrows = A.shape[0]
        nnz = A.nnz

        nlevels = 20
        bs = 9
        fill = 10
        droptol = 1e-3
        nmax = 3*nnz

        print("options: ", options)
        if options is not None:
           first = options.find('[')
           last  = options.find(']')
           if first != -1 and last != -1:
              opts = options[first+1:last]
              print(opts)
              if len(opts) > 0:
                 opts = opts.split(',')
                 for i,opt in enumerate(opts):
                    opt0 = opt.split('=')
                    print(opt0)
                    if len(opt0) > 1:
                       [key,val] = opt0[0:2]
                       if key == 'fill':
                          fill = int(val)
                       elif key == 'drop' or key == 'droptol':
                          droptol = float(val)
                       elif key == 'levels':
                          nlevels = int(val)
                       elif key == 'bs' or key == 'blocksize':
                          bs = int(val)
                       else:
                          print("Unknown option {}: {} {}".format(i,key,val))

        print("fill: {} droptol: {} nlevels: {} bs: {}".format(fill, droptol, nlevels, bs))

        jlu = np.empty( (nlevels*nmax), dtype=np.intc)
        ilu = np.empty( (nlevels*(mrows+1)), dtype=np.intc)
        alu = np.empty( (nlevels*nmax), dtype=np.double)

        ilu[:] = -999

        j1 = np.empty( (nmax), dtype=np.intc)
        i1 = np.empty( (mrows+1), dtype=np.intc)
        a1 = np.empty( (nmax), dtype=np.double)

        iwk = np.empty( (mrows), dtype=np.intc)
        iwk1 = np.empty( (mrows), dtype=np.intc)
        iwork = np.empty( (3*nmax), dtype=np.intc)
        iord = np.empty( (mrows), dtype=np.intc)
        riord = np.empty( (nlevels*mrows), dtype=np.intc)
        nlast = np.empty( (nlevels), dtype=np.intc)

        ipar = np.zeros( (10), dtype=np.intc )

        ipar[0] = nlevels # num levels requested
        ipar[1] = bs
        ipar[2] = fill
        ipar[5] = a1.shape[0] # size of A1
        ipar[6] = alu.shape[0] # size of Alu
        ipar[7] = nlevels # max number of levels

        build_func( ct.byref( ct.c_int(mrows) ),
                    np_pointer(A_data),
                    np_pointer(A_colidx),
                    np_pointer(A_rowptr),
                    np_pointer(alu),
                    np_pointer(jlu),
                    np_pointer(ilu),
                    np_pointer(iord),
                    np_pointer(riord),
                    np_pointer(nlast),
                    np_pointer(iwk),
                    np_pointer(iwk1),
                    np_pointer(iwork),
                    np_pointer(ipar),
                    ct.byref( ct.c_double(droptol) ),
                    np_pointer(a1),
                    np_pointer(j1),
                    np_pointer(i1) )

        print("ipar: ", ipar[-1])
        print(nlast)

        n = mrows
        ilu_ptr = 0
        print("Level, Size, Red, Black, Perm")
        for lvl in range(nlevels):
            left = n - nlast[lvl]
            print(lvl, n, nlast[lvl], left, riord[lvl], ilu[ilu_ptr+n])
            n = left
            ilu_ptr += 2*nlast[lvl] + left + 3

        idx = ilu.tolist().index(-999)
        print("ilu index: ", idx, ilu[idx], ilu[idx-1])

        self.alu = alu
        self.jlu = jlu
        self.ilu = ilu
        self.perm = riord
        self.nlast = nlast

        self.solve_func = solve_func

    def __call__( self, b, x = None ):
        return self.solve(b,x)

    # Solve Mx = b
    def psolve( self, b, x = None ):

        return self.solve(b, x)

    # Solve Mx = b
    def solve( self, b, x = None ):

        n = b.shape[0]
        nlevels = self.nlast.shape[0]
        rwk = np.empty((n),dtype=np.double)
        iwk = np.empty((n),dtype=np.intc)

        if x is None:
            x = np.empty_like(b)

        self.solve_func( ct.byref( ct.c_int(n) ), # n
                         np_pointer( b ),
                         np_pointer( x ),
                         np_pointer( rwk ),
                         np_pointer( self.alu ),
                         np_pointer( self.jlu ),
                         np_pointer( self.ilu ),
                         np_pointer( self.perm ),
                         np_pointer( self.nlast ),
                         ct.byref( ct.c_int(nlevels) ), # n
                         np_pointer( iwk ) )

        return x
