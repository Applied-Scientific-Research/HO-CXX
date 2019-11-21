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

shlib = load_lib( "./libindset_impl.so" )

def find_matrix_trisolve_parallelism(A):
    """
    Search the triangular CSR input matrix and return a sequence of rows that can be solved concurrently.
    """

    assert shlib is not None

    t_start = timestamp()

    assert (isinstance(A, sp.csr_matrix) or isinstance(A,sp.bsr_matrix))

    mrows = A.shape[0] # csr format
    nnz   = A.nnz

    lib_func = shlib.find_independent_sets

    lib_func.argtypes = [ ct.c_int, ct.c_int, # rows, nnz
                          ct.POINTER(ct.c_int), # A_rowptr(rows+1)
                          ct.POINTER(ct.c_int), # A_colidx(nnz)
                          ct.POINTER(ct.c_int) ]

    row_level = np.full( (mrows), dtype='i', fill_value=mrows )

    num_levels = lib_func( mrows, nnz, np_pointer(A.indptr), np_pointer(A.indices), np_pointer(row_level) )

    print(num_levels)

    ordering = np.argsort( row_level, kind='stable' )

    def index(arr, val):
       'Locate the leftmost value exactly equal to x'
       i = np.searchsorted( arr, val )
       return i

    row_level_ordered = row_level[ ordering ]

    level_offsets = np.zeros( (num_levels+1), dtype='i' )
    level_size = np.zeros( (num_levels), dtype='i' )
    for level in range(num_levels):
        i = index( row_level_ordered, level+1 )
        level_offsets[level+1] = i
        level_size[level] = -(i - level_offsets[level])
        if level < 10:
            print('level ', level, i, -level_size[level])

    ordered_level_size_index = np.argsort( level_size, kind='stable' )
    for k in range(min(10,len(ordered_level_size_index.tolist()))):
        j = ordered_level_size_index[k]
        print(k, j, abs(level_size[j]))

    return num_levels, ordering, level_offsets
