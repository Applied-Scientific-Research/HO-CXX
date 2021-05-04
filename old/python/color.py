import numpy as np
import scipy.sparse as sp
import time

def timestamp():
    return time.time()

def color_matrix( Ai ):

    print("Starting matrix coloring")

    time_start = timestamp()

    A = None
    if not (isinstance(Ai, sp.csr_matrix) or isinstance(Ai, sp.bsr_matrix)):
        A = sp.csr_matrix(Ai)
    else:
        A = Ai

    assert A.shape[0] == A.shape[1]

    n = A.shape[0]

    if isinstance(A, sp.bsr_matrix):
        assert A.blocksize[0] == A.blocksize[1]
        n = n // A.blocksize[0]

    num_colors = 0

    # initial the vertex colors so that all will get lumped with the last row
    # to be processed. These will be updated during the analysis.
    colors = np.full( (n), fill_value=(n-1), dtype='i' )

    mark = np.full( (n), fill_value=n, dtype='i' )

    for row in range(n):

        #for k in range( A.indptr[row], A.indptr[row+1] ):
        #    col = A.indices[k]
        #    mark[colors[col]] = row
        mark[colors[ A.indices[ A.indptr[row]:A.indptr[row+1] ]]] = row

        # Find the color of this vertex
        row_color = 0
        while row_color < num_colors and mark[row_color] == row:
            row_color += 1

        # if we didn't hit an existing color, create a new one.
        if row_color == num_colors: 
            num_colors += 1

        colors[row] = row_color

    time_stop = timestamp()

    print("Finished matrix coloring: colors= {} time= {}".format( num_colors, (time_stop-time_start)) )

    return num_colors, colors

def coloring_order( A ):

    num_colors, row_color = color_matrix(A)

    #print(row_color[:20])

    # Sort the rows by color.
    ordering = np.argsort( row_color, kind='stable' )

    #print(ordering[:20])

    #for i in range(20):
    #    print(i, ordering[i], row_color[ ordering[i] ])

    # find the offsets into the ordering for each color.

    #import bisect

    def index(arr, val):
       'Locate the leftmost value exactly equal to x'
       i = np.searchsorted( arr, val )
       return i
       #i = bisect_left(a, x)
       #if i != len(a) and a[i] == x:
       #    return i
       #raise ValueError

    row_color_ordered = row_color[ ordering ]

    color_offsets = np.zeros( (num_colors+1), dtype='i' )
    color_size = np.zeros( (num_colors), dtype='i' )
    for color in range(num_colors):
        i = index( row_color_ordered, color+1 )
        color_offsets[color+1] = i
        color_size[color] = -(i - color_offsets[color])
        #print('color ', color, i, i-color_offsets[color])

    ordered_color_size_index = np.argsort( color_size, kind='stable' )
    for k in range(min(10,len(ordered_color_size_index.tolist()))):
        j = ordered_color_size_index[k]
        print(k, j, abs(color_size[j]))

    return num_colors, ordering, color_offsets
