#ifndef __multiply_h
#define __multiply_h

#include <splib/vector.h>
#include <splib/base_matrix.h>

// Template function for Y := A * X or A^T * X

namespace splib
{

template <typename Matrix, typename MatrixOrVector1, typename MatrixOrVector2>
void multiply (const Matrix&          A,
               const MatrixOrVector1& X,
                     MatrixOrVector2& Y,
               const bool transposeA = false);

// z := alpha Ax + beta y
template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
void
gemv (const typename Vector3::value_type	alpha,
      const Matrix&	A,
      const Vector1&	x,
      const typename Vector3::value_type	beta,
      const Vector2&	y,
            Vector3&	z);

void multiply (BaseMatrix *A, BaseVector *x, BaseVector *y);
void gemv (double alpha, BaseMatrix *A, BaseVector *x, double beta, BaseVector *y, BaseVector *z);

} // namespace splib

#endif
