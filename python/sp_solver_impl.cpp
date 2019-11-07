#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace details
{

template <typename T>
int sp_ztrsv( const int mrows, const int ncols, const int nnz,
               const bool isUpper, const bool isUnit, const bool isSorted,
               int A_indptr[], int A_indices[], T A_data[],
               T x[] )
{
    if (isUpper)
    {
        // Start at the bottom and work backwards to solve Ux = y

        for (int i = mrows-1; i >= 0; --i)
        {
            const int k0 = A_indptr[i];
            const int k1 = A_indptr[i+1];

            // The diagonal is at the first term.
            if (isSorted)
            {
                const int kdiag = k0;
                double xsum = x[i]; // b

                for (int k = k0+1; k < k1; ++k)
                {
                    const int j = A_indices[k];
                    xsum -= A_data[k] * x[j];
                }

                if (isUnit)
                    x[i] = xsum;
                else
                    x[i] = xsum / A_data[kdiag];
            }
            else
            {
                int kdiag = mrows;
                double xsum = x[i]; // b

                for (int k = k0; k < k1; ++k)
                {
                    const int j = A_indices[k];
                    if (j == i) // diagonal
                        kdiag = k;
                    else
                        xsum -= A_data[k] * x[j];
                }

                //x[i] = xsum / A_data[kdiag];
                if (isUnit)
                    x[i] = xsum;
                else
                    x[i] = xsum / A_data[kdiag];
            }
        }
    }
    else
    {
        for (int i = 0; i < mrows; ++i)
        {
            const int k0 = A_indptr[i];
            const int k1 = A_indptr[i+1];

            // The diagonal is at the last term.
            if (isSorted)
            {
                const int kdiag = k1-1;
                double xsum = x[i]; // b

                for (int k = k0; k < (k1-1); ++k)
                {
                    const int j = A_indices[k];
                    xsum -= A_data[k] * x[j];
                }

                //x[i] = xsum / A_data[kdiag];
                if (isUnit)
                    x[i] = xsum;
                else
                    x[i] = xsum / A_data[kdiag];
            }
            else
            {
                int kdiag = mrows;
                double xsum = x[i]; // b

                for (int k = k0; k < k1; ++k)
                {
                    const int j = A_indices[k];
                    if (j == i) // diagonal
                        kdiag = k;
                    else
                        xsum -= A_data[k] * x[j];
                }

                //x[i] = xsum / A_data[kdiag];
                if (isUnit)
                    x[i] = xsum;
                else
                    x[i] = xsum / A_data[kdiag];
            }
        }
    }
}

} // namespace-details

extern "C"
{

int sp_dtrsv( const int mrows, const int ncols, const int nnz,
               const int isUpper, const int isUnit, const int isSorted,
               int A_indptr[], int A_indices[], double A_data[],
               double x[] )
{
    if (isUpper)
        if (isUnit)
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 1, 1, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 1, 1, 0, A_indptr, A_indices, A_data, x );
        else
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 1, 0, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 1, 0, 0, A_indptr, A_indices, A_data, x );
    else
        if (isUnit)
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 0, 1, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 0, 1, 0, A_indptr, A_indices, A_data, x );
        else
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 0, 0, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 0, 0, 0, A_indptr, A_indices, A_data, x );
}

int sp_strsv( const int mrows, const int ncols, const int nnz,
               const int isUpper, const int isUnit, const int isSorted,
               int A_indptr[], int A_indices[], float A_data[],
               float x[] )
{
    if (isUpper)
        if (isUnit)
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 1, 1, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 1, 1, 0, A_indptr, A_indices, A_data, x );
        else
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 1, 0, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 1, 0, 0, A_indptr, A_indices, A_data, x );
    else
        if (isUnit)
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 0, 1, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 0, 1, 0, A_indptr, A_indices, A_data, x );
        else
            if (isSorted) return details::sp_ztrsv( mrows, ncols, nnz, 0, 0, 1, A_indptr, A_indices, A_data, x );
            else          return details::sp_ztrsv( mrows, ncols, nnz, 0, 0, 0, A_indptr, A_indices, A_data, x );
}

} // extern-C
