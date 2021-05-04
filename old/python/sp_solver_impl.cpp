#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

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
                T xsum = x[i]; // b

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
                T xsum = x[i]; // b

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
                T xsum = x[i]; // b

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
                T xsum = x[i]; // b

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

#define A(i,j) ( _A[ (i)*m + (j) ] )

// y := alpha Ax + beta y
template <typename T>
void gemv( const int _m, const T alpha, const T _A[], const T x[], const T beta, T y[] )
{
    const int m = 4;

    for (int i = 0; i < m; ++i)
    {
        T yi = 0.0;

        #pragma ivdep
        for (int j = 0; j < m; ++j)
            yi += A(i,j) * x[j];

        if (beta == T(0))
            y[i] = alpha * yi;
        else
            y[i] = alpha * yi + beta * y[i];
    }

    return;
}

// Solve Lx = y or Ux = y
template <typename T>
void trsv( const bool isUpper, const bool isUnit, const int _m, const T _A[], T x[] )
{
    const int m = 4;

   // A is Upper triangular; solve Ux = b
   if (isUpper)
      for (int j = m-1; j >= 0; --j)
      {
         T x_j = x[j];
         if (x_j != T(0))
         {
            if ( not(isUnit) ) { x[j] /= A(j,j); x_j = x[j]; }

            #pragma ivdep
            for (int i = 0; i < j; ++i)
                x[i] -= x_j * A(i,j);
         }
      }
   // A is Lower triangular; solve Lx = b
   else
      for (int j = 0; j < m; ++j)
      {
         T x_j = x[j];
         if (x_j != T(0))
         {
            if ( not(isUnit) ) x_j = ( x[j] /= A(j,j) );

            #pragma ivdep
            for (int i = j+1; i < m; ++i)
               x[i] -= x_j * A(i,j);
         }
      }

    return;
}

#undef A

template <typename T>
int sp_ztrsv_bsr( const int mrows, const int ncols, const int nnz, const int bs,
                  const bool isUpper, const bool isUnit,
                  const int A_indptr[], const int A_indices[], const T A_data[],
                  T x[] )
{
    if ( mrows % bs ) { fprintf(stderr,"mrows %% bs = %d %d %d\n", mrows % bs, mrows, bs); exit(1); }

    const int brows = mrows / bs;
    const int bs2 = bs*bs;

    //printf("ztrsv_bsr: %s %d %d %d %d %d %d %d %d\n", typeid(T).name(), isUpper, isUnit, mrows, ncols, nnz, bs, brows, A_indptr[brows]);

    for (int ib = 0; ib < brows; ++ib)
    {
        const int i = isUpper ? ((brows-1) - ib) : ib;

        const int k0 = A_indptr[i];
        const int k1 = A_indptr[i+1];

        // Factor out the off-diagonal blocks first.
        int kdiag = k0;
        T *xi = x + i * bs;

        for (int k = k0; k < k1; ++k)
        {
            const int j = A_indices[k];
            if (j == i)
                kdiag = k;
            else
            {
                const T *Aij = A_data + k * bs2;
                const T *xj  = x + j * bs;

                // x_i = b_i - Sum_j{ A_i,j * x_j }
                gemv( bs, T(-1), Aij, xj, T(1), xi );
            }
        }

        // Now solve Ax = y on the diagonal block.

        const T* Aii = A_data + kdiag * bs2;
        trsv( isUpper, isUnit, bs, Aii, xi );
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

int sp_strsv_bsr( const int mrows, const int ncols, const int nnz, const int bs,
               const int isUpper, const int isUnit,
               int A_indptr[], int A_indices[], float A_data[],
               float x[] )
{
    if (isUpper)
        if (isUnit)
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 1, 1, A_indptr, A_indices, A_data, x );
        else
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 1, 0, A_indptr, A_indices, A_data, x );
    else
        if (isUnit)
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 0, 1, A_indptr, A_indices, A_data, x );
        else
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 0, 0, A_indptr, A_indices, A_data, x );
}

int sp_dtrsv_bsr( const int mrows, const int ncols, const int nnz, const int bs,
               const int isUpper, const int isUnit,
               int A_indptr[], int A_indices[], double A_data[],
               double x[] )
{
    if (isUpper)
        if (isUnit)
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 1, 1, A_indptr, A_indices, A_data, x );
        else
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 1, 0, A_indptr, A_indices, A_data, x );
    else
        if (isUnit)
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 0, 1, A_indptr, A_indices, A_data, x );
        else
            return details::sp_ztrsv_bsr( mrows, ncols, nnz, bs, 0, 0, A_indptr, A_indices, A_data, x );
}

} // extern-C
