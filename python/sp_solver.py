import numpy as np
import scipy.sparse as sp
import ctypes as ct
import ctypes.util
import sys, platform
import time

def timestamp():
    return time.time()

def ct_pointer( arr ):

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

class sp_solver:

    """ Python wrapper for the C shared library lib"""
 
    def __init__( self, L, U, pr, pc ):

        lib = load_lib( "libsp_solver_impl.so" )

        c_double_p = ct.POINTER(ct.c_double)
        c_float_p = ct.POINTER(ct.c_float)
        c_int_p = ct.POINTER(ct.c_int)

        self.mixed_precision = False

        assert L.data.dtype == U.data.dtype
        self.mixed_precision = (L.data.dtype == np.float32)

        self.sp_dtrsv_wrapper = lib.sp_dtrsv
        self.sp_dtrsv_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int,
                                              ct.c_int, ct.c_int, ct.c_int,
                                              c_int_p, c_int_p, c_double_p,
                                              c_double_p ]
        #self.sp_dtrsv_wrapper.restype = [ ct.c_int ]

        print("Loaded sp_dtrsv routine: ", lib.sp_dtrsv)

        self.sp_strsv_wrapper = lib.sp_strsv
        self.sp_strsv_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int,
                                              ct.c_int, ct.c_int, ct.c_int,
                                              c_int_p, c_int_p, c_float_p,
                                              c_float_p ]

        self.sp_strsv_bsr_wrapper = lib.sp_strsv_bsr
        self.sp_strsv_bsr_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int, ct.c_int, # bs
                                               ct.c_int, ct.c_int,
                                               c_int_p, c_int_p, c_float_p,
                                               c_float_p ]

        self.sp_dtrsv_bsr_wrapper = lib.sp_dtrsv_bsr
        self.sp_dtrsv_bsr_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int, ct.c_int, # bs
                                               ct.c_int, ct.c_int,
                                               c_int_p, c_int_p, c_double_p,
                                               c_double_p ]

        self.use_mkl = False

        #libmkl = load_lib( "libmkl_sparse_impl.so" )
        #if libmkl is not None:
        #    self.sp_mkl_dtrsv_wrapper = libmkl.sp_mkl_dtrsv
        #    self.sp_mkl_dtrsv_wrapper.argtypes = [ ct.c_int, ct.c_int, ct.c_int,
        #                                          ct.c_int, ct.c_int, ct.c_int,
        #                                          c_int_p, c_int_p, c_double_p,
        #                                          c_double_p, c_double_p ]

        #    self.use_mkl = True

        if self.use_mkl:
            print("Loaded sp_mkl_dtrsv routine: ", libmkl.sp_mkl_dtrsv)
        else:
            print("Disabled mkl routines")

        t_start = timestamp()

        self.L = L
        self.U = U

        print("ILU with dtype= ", self.L.dtype, self.mixed_precision)

        self.L_bsr = None
        self.U_bsr = None
        if True:
           bs = 4
           print(bs)
           L_bsr = sp.bsr_matrix( self.L, blocksize=(bs,bs))
           U_bsr = sp.bsr_matrix( self.U, blocksize=(bs,bs))
           print("L_bsr: {} {} {} {:.2f} {:.2f}".format(L_bsr.nnz, L_bsr.blocksize, L_bsr.indptr[-1], 100.*float(L_bsr.nnz)/(self.L.shape[0]**2), L_bsr.nnz/self.L.nnz))
           print("U_bsr: {} {} {} {:.2f} {:.2f}".format(U_bsr.nnz, U_bsr.blocksize, U_bsr.indptr[-1], 100.*float(U_bsr.nnz)/(self.U.shape[0]**2), U_bsr.nnz/self.U.nnz))
           #print(M_L_bsr.data[0])
           #print(linalg.inv(M_L_bsr.data[0]))
           #print(M_U_bsr.data[0])
           #print(linalg.inv(M_U_bsr.data[0]))

           self.L_bsr = L_bsr
           self.U_bsr = U_bsr

        self.pr = pr
        self.pc = pc

        #assert self.L.data.dtype == np.double
        #assert self.L.indices.dtype == np.intc and L.indptr.dtype == np.intc
        #assert self.U.data.dtype == np.double
        #assert self.U.indices.dtype == np.intc and U.indptr.dtype == np.intc

        # Test of the columns are sorted (in rows)
        self.L_is_sorted = True
        #for i in range(self.L.shape[0]):
        #    #k0 = self.L.indptr[i]
        #    k1 = self.L.indptr[i+1]
        #    #row_is_sorted = np.all(self.L.indices[k0:k1-1] < self.L.indices[k0+1:k1])
        #    #if not (row_is_sorted and self.L.indices[k1-1] == i):
        #    if not (self.L.indices[k1-1] == i):
        #        self.L_is_sorted = False
        #        break

        self.U_is_sorted = True
        #for i in range(self.U.shape[0]):
        #    k0 = self.U.indptr[i]
        #    #k1 = self.U.indptr[i+1]
        #    #row_is_sorted = np.all(self.U.indices[k0:k1-1] < self.U.indices[k0+1:k1])
        #    #if not (row_is_sorted and self.U.indices[k0] == i):
        #    if not (self.U.indices[k0] == i):
        #        self.U_is_sorted = False
        #        break

        #t_finish = timestamp()
        #print("Finished analyzing L/U matrices: ", self.L_is_sorted, self.U_is_sorted, 1000.*(t_finish-t_start))

        #L_levels = self.find_parallelism(self.L)
        #U_levels = self.find_parallelism(self.U)

    def find_parallelism(self, A):

        t_start = timestamp()

        assert isinstance(A,sp.csr_matrix)

        mrows = A.shape[0] # csr format
        nnz   = A.nnz

        A_csc = sp.csc_matrix(A)

        A_nnzs = np.empty( (mrows), dtype='i' )
        #A_column_indices = []
        for i in range(mrows):
            k      = A.indptr[i]
            k_stop = A.indptr[i+1]

            A_nnzs[i] = k_stop - k

            #cols = A.indices[ k:k_stop ]
            #A_column_indices.append( list(cols) )

        row_indices = np.arange(mrows)

        maxiters = 10000
        maxiters = mrows
        iters = 0

        levels = []

        non_zeros = nnz
        active_rows = mrows

        while (iters < maxiters) and active_rows > 0:

            # Find any rows with only 1 dependency (that could be the diagonal).
            one_connection = (A_nnzs == 1)

            rows_to_remove = row_indices[ one_connection ]

            #print("{}: rows with 1 connection: {}".format(iters, rows_to_remove.shape[0]))

            levels.append( rows_to_remove )

            # Clear the current candidates
            for row in rows_to_remove:

                # For each row 'i' to be removed, sweep down the 'i' column of the CSC matrix
                # to find which other rows are connected and remove them.

                k_begin = A_csc.indptr[row  ]
                k_end   = A_csc.indptr[row+1]

                for k in range(k_begin, k_end):
                    other_row = A_csc.indices[k]
                    A_nnzs[ other_row ] -= 1
                    #A_column_indices[other_row].remove( row )
                    #k += 1
                    non_zeros -= 1

            #new_nnz = 0
            #for i in range(mrows):
            #    if A_nnzs[i] != len(A_column_indices[i]):
            #        print("mismatch ", i, A_nnzs[i], len(A_column_indices[i]))
            #        sys.exit(2)
            #    new_nnz += len(A_column_indices[i])

            active_rows = np.count_nonzero(A_nnzs)
            #print(iters, rows_to_remove.shape[0], non_zeros, active_rows)

            iters += 1

        nlevels = len(levels)
        nrows = np.empty( (nlevels), dtype='i' )
        for i,lvl in enumerate(levels):
            nrows[i] = len(lvl)

        t_stop = timestamp()
        t_ms = 1000.*(t_stop - t_start)
        if active_rows > 0:
            print("Failed to complete analysis in {:.2f} ms: {} levels with min/max {} / {}".format( t_ms, nlevels, nrows.min(), nrows.max()))
            return []
        else:
            print("Finished analysis in {:.2f} ms: {} levels with min/max {} / {}".format( t_ms, nlevels, nrows.min(), nrows.max()))
            return levels

    def __call__( self, x, y = None ):
        return self.solve(x,y)

    def solve( self, x, y = None ):

        assert x.dtype == np.double
        if y is not None:
            assert x.dtype == y.dtype

        assert self.sp_dtrsv_wrapper is not None
        if self.mixed_precision:
            assert self.sp_strsv_wrapper is not None
        if self.use_mkl:
            assert self.sp_mkl_dtrsv_wrapper is not None

        mrows = self.L.shape[0]
        ncols = self.L.shape[1]
        l_nnz = self.L.nnz
        u_nnz = self.U.nnz

        px = None
        if self.mixed_precision:
            px = np.empty( (mrows), dtype=np.float32)
        else:
            px = np.empty_like(x)

        # Permute the rhs vector.
        px[ self.pr[:] ] = x[:]

        pz = None
        if self.use_mkl:
            pz = np.empty_like(px)

        func = self.sp_dtrsv_wrapper
        if self.mixed_precision:
            func = self.sp_strsv_wrapper

        bsr_func = self.sp_dtrsv_bsr_wrapper
        if self.mixed_precision:
            bsr_func = self.sp_strsv_bsr_wrapper

        # (Pr)Ly = (Pr)x
        isUpper = 0
        isUnit  = 1
        isSorted = 0
        if self.L_is_sorted:
            isSorted = 1

        if self.use_mkl:
            info = self.sp_mkl_dtrsv_wrapper (
                         mrows, ncols, l_nnz, isUpper, isUnit, isSorted,
                         ct_pointer( self.L.indptr ), ct_pointer( self.L.indices ), ct_pointer( self.L.data ),
                         ct_pointer( px ), ct_pointer( pz ) )
        else:
            if self.L_bsr is not None:
                info = bsr_func (
                         mrows, ncols, self.L_bsr.nnz, self.L_bsr.blocksize[0],
                         isUpper, isUnit,
                         ct_pointer( self.L_bsr.indptr ), ct_pointer( self.L_bsr.indices ), ct_pointer( self.L_bsr.data ),
                         ct_pointer( px ) )
            else:
                #px_tmp = np.copy(px)
                info = func (
                         mrows, ncols, l_nnz, isUpper, isUnit, isSorted,
                         ct_pointer( self.L.indptr ), ct_pointer( self.L.indices ), ct_pointer( self.L.data ),
                         ct_pointer( px ) )

                #info = bsr_func (
                #         mrows, ncols, self.L_bsr.nnz, self.L_bsr.blocksize[0],
                #         isUpper, isUnit,
                #         ct_pointer( self.L_bsr.indptr ), ct_pointer( self.L_bsr.indices ), ct_pointer( self.L_bsr.data ),
                #         ct_pointer( px_tmp ) )

                #print(np.linalg.norm( px - px_tmp ))

        # Ux' = y
        isUpper = 1
        isUnit  = 0
        isSorted = 0
        if self.U_is_sorted:
            isSorted = 1

        if self.use_mkl:
            info = self.sp_mkl_dtrsv_wrapper (
                         mrows, ncols, u_nnz, isUpper, isUnit, isSorted,
                         ct_pointer( self.U.indptr ), ct_pointer( self.U.indices ), ct_pointer( self.U.data ),
                         ct_pointer( pz ), ct_pointer( px ) )
        else:
            if self.U_bsr is not None:
                info = bsr_func (
                         mrows, ncols, l_nnz, self.U_bsr.blocksize[0],
                         isUpper, isUnit,
                         ct_pointer( self.U_bsr.indptr ), ct_pointer( self.U_bsr.indices ), ct_pointer( self.U_bsr.data ),
                         ct_pointer( px ) )
            else:
                info = func (
                         mrows, ncols, u_nnz, isUpper, isUnit, isSorted,
                         ct_pointer( self.U.indptr ), ct_pointer( self.U.indices ), ct_pointer( self.U.data ),
                         ct_pointer( px ) )

        # x = (Pc)x'
        if y is None:
            y = np.empty_like(x)

        y[:] = px[ self.pc[:] ]

        return y
