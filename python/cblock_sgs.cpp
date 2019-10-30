#include <cmath>

namespace details
{

// y := alpha*A*x + beta*y
template <int _bs>
struct dgemv
{
    enum : int { bs = _bs };

    void _dgemv( const int n, const double alpha, double A[], double x[], const double beta, double y[] )
    {
        for (int i = 0; i < n; ++i)
        {
            double yi = 0.0;
            double *Ai = A + i*n;

            #pragma ivdep
            for (int j = 0; j < n; ++j)
                yi += Ai[j] * x[j];
    
            if (beta == 0.0)
                y[i] = alpha * yi;
            else
                y[i] = alpha * yi + beta * y[i];
        }
    
        return;
    }

    void operator()( const int n, const double alpha, double A[], double x[], const double beta, double y[] )
    {
        if (bs)
            this->_dgemv(bs, alpha, A, x, beta, y );
        else
            this->_dgemv( n, alpha, A, x, beta, y );
    }

};

template <class Functor>
void cblock_sgs_impl( const int nrows, const int bs, const int niters,
                      int L_indptr[], int L_indices[], double L_data[],
                      int U_indptr[], int U_indices[], double U_data[],
                      double Dinv_data[], double x[], double y[], Functor dgemv )
{
    const int nbrows = nrows / bs;
    const int bs2 = bs*bs;

    double yprime[bs];

    for (int iter = 0; iter < niters; iter++)
    {

    // Forward sweep: (D + L) y^(k+1/2) = x - U y^(k)

    for (int ib = 0; ib < nbrows; ib++)
    {
        double *yib = y + ib * bs;
        double *xib = x + ib * bs;

        // y' = x - Uy
        for (int j = 0; j < bs; ++j)
            yprime[j] = xib[j];

            // y^0 = 0 so we skip the U*y term.
        if (iter != 0)
        {
            const int k0 = U_indptr[ib];
            const int k1 = U_indptr[ib+1];

            for (int k = k0; k < k1; ++k)
            {
                const int j = U_indices[k];
                dgemv( bs, -1.0, U_data + (k*bs2), y + j*bs, 1.0, yprime );
            }
        }

        // Dy^(k+1/2) = [ x - Uy ] - Ly
        {
            const int k0 = L_indptr[ib];
            const int k1 = L_indptr[ib+1];

            for (int k = k0; k < k1; ++k)
            {
                const int j = L_indices[k];
                dgemv( bs, -1.0, L_data + (k*bs2), y + j*bs, 1.0, yprime );
            }
        }

        dgemv( bs, 1.0, Dinv_data + (ib*bs2), yprime, 0.0, yib );
    }

    // Backward sweep: (D + U) y^(k+1) = x - L y^(k+1/2) )

    for (int ib = nbrows-1; ib >= 0; ib--)
    {
        double *yib = y + ib * bs;
        double *xib = x + ib * bs;

        // y' = x - Ly
        for (int j = 0; j < bs; ++j)
            yprime[j] = xib[j];

        {
            const int k0 = L_indptr[ib];
            const int k1 = L_indptr[ib+1];

            for (int k = k0; k < k1; ++k)
            {
                const int j = L_indices[k];
                dgemv( bs, -1.0, L_data + (k*bs2), y + j*bs, 1.0, yprime );
            }
        }

        {
            const int k0 = U_indptr[ib];
            const int k1 = U_indptr[ib+1];

            // D y^(k+1) = [x - Ly^(k+1/2)] - Uy^(k+1) ... = y' - Uy
            for (int k = k0; k < k1; ++k)
            {
                const int j = U_indices[k];
                dgemv( bs, -1.0, U_data + (k*bs2), y + j*bs, 1.0, yprime );
            }
        }

        // y^(k+1) = D^(-1) y'
        dgemv( bs, 1.0, Dinv_data + (ib*bs2), yprime, 0.0, yib );
    }

    } // iters

    return;
}

} // namespace

extern "C"
void cblock_sgs_impl( const int nrows, const int bs, const int niters,
                      int L_indptr[], int L_indices[], double L_data[],
                      int U_indptr[], int U_indices[], double U_data[],
                      double Dinv_data[], double x[], double y[] )
{
    if (bs == 9)
        details::cblock_sgs_impl( nrows, 9, niters,
                      L_indptr, L_indices, L_data,
                      U_indptr, U_indices, U_data,
                      Dinv_data, x, y, details::dgemv<9>() );
    else
        details::cblock_sgs_impl( nrows, bs, niters,
                      L_indptr, L_indices, L_data,
                      U_indptr, U_indices, U_data,
                      Dinv_data, x, y, details::dgemv<0>() );

    return;
}
