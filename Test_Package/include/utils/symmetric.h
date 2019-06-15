#ifndef __has_utils_symmetric_h
#define __has_utils_symmetric_h

#include <splib/formats.h>
#include <splib/matrix_utils.h>

#include <utils/timer.h>

namespace utils
{

namespace details
{

template <class Matrix, typename ValueType>
size_t _symmetry (Matrix& A, ValueType &symmetry, splib::csr_format)
{
   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   value_type norm_a = 0.0;
   value_type norm_amat = 0.0;

   size_t num_missed = 0;

   for (index_type i = 0; i < A.num_rows(); ++i)
   {
      value_type a_sum_j = 0.0;
      value_type amat_sum_j = 0.0;

      for (index_type k = A.m_rowptr[i]; k < A.m_rowptr[i+1]; k++)
      {
         index_type j = A.m_colidx[k];
         value_type a_ij = A.m_values[k];
         value_type a_ji(0);

         index_type idx = A.find(j,i);
         if (idx != A.num_values())
            a_ji = A.m_values[idx];
         else
            num_missed++;

         amat_sum_j += fabs(a_ij - a_ji);
         a_sum_j += fabs(a_ij);
      }

      norm_a = std::max(norm_a, a_sum_j);
      norm_amat = std::max(norm_amat, amat_sum_j);
   }

   return num_missed;
}
template <class Matrix>
size_t _symmetry_pattern (Matrix& A, splib::csr_format)
{
   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   size_t num_missed = 0;

   for (index_type i = 0; i < A.num_rows(); ++i)
   {
      for (index_type k = A.m_rowptr[i]; k < A.m_rowptr[i+1]; k++)
      {
         index_type j = A.m_colidx[k];

         index_type idx = A.find(j,i);
         if (idx == A.num_values())
            num_missed++;
      }
   }

   return num_missed;
}

template <class Matrix, typename ValueType>
size_t _symmetry (Matrix& A, ValueType &symmetry, splib::coo_format)
{
   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   splib::CSR<index_type,value_type> A_csr;
   A_csr = A;

   return _symmetry(A_csr, symmetry, splib::csr_format());
}
template <class Matrix>
size_t _symmetry_pattern (Matrix& A, splib::coo_format)
{
   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   size_t num_missed = 0;

   for (index_type k = 0; k < A.num_values(); ++k)
   {
      const index_type i = A.m_rowidx[k];
      const index_type j = A.m_colidx[k];
      index_type idx = A.find(j,i);
      if (idx == A.num_values())
         num_missed++;
   }

   return num_missed;
}

}

// Test the symmetric of the matrix. || A - A^T || / ||A||
// ... using the ||A||_inf = max_i Sum_j |a_ij|

template <class Matrix>
bool test_symmetric (const Matrix& A, typename Matrix::value_type tol = 1.0e-10)
{
   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   value_type a_norm;

   size_t num_missed = details::_symmetry (A, a_norm, typename Matrix::format_type() );
   printf("||A - A^T|| / ||A|| (inf-norm) = %e %d\n", a_norm, num_missed);

   return (a_norm < tol);
}
template <class Matrix>
bool test_symmetric_pattern (const Matrix& A)
{
   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   return (details::_symmetry_pattern (A, typename Matrix::format_type()) == 0);
}

} // namespace utils

#endif
