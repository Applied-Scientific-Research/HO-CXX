#include <cmath>
#include <cstdio>
#include <cstdlib>

extern "C"
int sp_dgetrs( const int mrows, const int ncols, const int nnz,
               const int IsUpper, const int IsUnit,
               int A_indptr[], int A_indices[], double A_data[],
               double x[] )
{
    if (IsUpper)
    {
        // Start at the bottom and work backwards to solve Ux = y

        for (int i = mrows-1; i >= 0; --i)
        {
            int kdiag = mrows;
            double xsum = x[i]; // b

            const int k0 = A_indptr[i];
            const int k1 = A_indptr[i+1];

            for (int k = k0; k < k1; ++k)
            {
                const int j = A_indices[k];
                if (j == i) // diagonal
                    kdiag = k;
                else if (j < i) // ensure it's a upper terms
                    { kdiag = -k; break; }
                else
                    xsum -= A_data[k] * x[j];
            }

            if (kdiag == mrows)
                { fprintf(stderr,"No diagonal found on row %d\n", i); exit(-1); }
            else if (kdiag < 0)
                { fprintf(stderr,"Lower diagonal term found in Upper matrix on row %d %d\n", i, A_indices[-kdiag]); exit(-1); }
            else
                x[i] = xsum / A_data[kdiag];
        }
    }
    else
    {
        for (int i = 0; i < mrows; ++i)
        {
            int kdiag = mrows;
            double xsum = x[i]; // b

            const int k0 = A_indptr[i];
            const int k1 = A_indptr[i+1];

            for (int k = k0; k < k1; ++k)
            {
                const int j = A_indices[k];
                if (j == i) // diagonal
                    kdiag = k;
                else if (j > i) // ensure it's a lower term
                    { kdiag = -k; break; }
                else
                    xsum -= A_data[k] * x[j];
            }

            if (kdiag == mrows)
                { fprintf(stderr,"No diagonal found on row %d\n", i); exit(-1); }
            else if (kdiag < 0)
                { fprintf(stderr,"Upper diagonal term found in Lower matrix on row %d %d\n", i, A_indices[-kdiag]); exit(-1); }
            else
                x[i] = xsum / A_data[kdiag];
        }
    }
}
