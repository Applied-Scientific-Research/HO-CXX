#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>
#include <type_traits>

#ifdef WITH_MKL
#include "mkl_spblas.h"

extern "C"
int sp_mkl_dtrsv( const int mrows, const int ncols, const int nnz,
                  const int isUpper, const int isUnit, const int isSorted,
                  int A_indptr[], int A_indices[], double A_data[],
                  double x[], double y[] )
{
    static_assert( std::is_same<int,MKL_INT>::value, "Type mismatch in MKL interface");

    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t csrA;

    // Create handle with matrix stored in CSR format
    auto status = mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO,
                                    mrows,  // number of rows
                                    ncols,  // number of cols
                                    A_indptr,
                                    A_indptr+1,
                                    A_indices,
                                    A_data );
    if ( status != SPARSE_STATUS_SUCCESS )
    {
        fprintf(stderr,"Error in mkl_sparse_d_create_csr %d\n", status);
        exit(1);
    }

    // Create matrix descriptor
    struct matrix_descr descrA;

    descrA.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrA.mode = ( isUpper ) ? SPARSE_FILL_MODE_UPPER : SPARSE_FILL_MODE_LOWER;
    descrA.diag = ( isUnit ) ? SPARSE_DIAG_UNIT : SPARSE_DIAG_NON_UNIT;

    // Compute y = alpha * A^{-1} * x
    const double alpha = 1.0;
    status = mkl_sparse_d_trsv ( SPARSE_OPERATION_NON_TRANSPOSE,
                             alpha,
                             csrA, descrA,
                             x, y );

    // Release matrix structure
    mkl_sparse_destroy ( csrA );

    if ( status != SPARSE_STATUS_SUCCESS )
    {
        fprintf(stderr,"Error in mkl_sparse_d_trsv %d\n", status);
        exit(-2);
    }

    return 0;
}

#endif
