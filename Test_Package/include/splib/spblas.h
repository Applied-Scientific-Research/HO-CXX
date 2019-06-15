#ifndef __has_spblas_h
#define __has_spblas_h

namespace splib
{
namespace spblas
{

/* Vector object interfaces */

// y := alpha * x + y
template <typename Vector1,typename Vector2>
void
axpy (const typename Vector1::value_type alpha,
      const Vector1&   x,
            Vector2&   y);

// z := alpha * x + beta * y
template <typename Vector1, typename Vector2, typename Vector3>
void
axpby (const typename Vector1::value_type alpha,
       const typename Vector2::value_type beta,
       const Vector1&   x,
       const Vector2&   y,
             Vector3&   z);

// z := x * y
template <typename Vector1, typename Vector2, typename Vector3>
void
xmy (const Vector1& x, const Vector2& y, Vector3& z);

// y := alpha * x
template <typename Vector1, typename Vector2>
void
scale (const typename Vector1::value_type alpha, const Vector1& x, Vector2& y);

// y := x
template <typename Vector1, typename Vector2>
void
copy (const Vector1& x, Vector2& y);

// x := alpha
template <typename Vector>
void
fill (const typename Vector::value_type alpha, Vector& x);

// y := 1/x
template <typename Vector1, typename Vector2>
void
recip (const Vector1& x, Vector2& y);

// Inner product (i.e., x^T * y)
template <typename Vector1, typename Vector2>
typename Vector1::value_type
dot (const Vector1& x, const Vector2& y);

// ||x||_2
template <typename Vector>
typename Vector::value_type
norm (const Vector& x);

// ||x||_2
template <typename Vector>
typename Vector::value_type
norm2 (const Vector& x);

// ||x||_inf
template <typename Vector>
typename Vector::value_type
norm_inf (const Vector& x);

// max(x_i)
template <typename Vector>
typename Vector::value_type
max(const Vector& x);

// min(x_i)
template <typename Vector>
typename Vector::value_type
min(const Vector& x);

// fill with linear sequence.
template <typename Vector>
void
sequence (Vector& x, typename Vector::value_type x0 = typename Vector::value_type(0));

#if 0
// x[] = y[ p[] ]
template <typename Vector1, typename Vector2>
void
gather (const Vector1& y, const Vector2& p, Vector1& x);

// y[ p[] ] = x[]
template <typename Vector1, typename Vector2>
void
scatter (const Vector1& x, const Vector2& p, Vector1& y);
#endif

// permute the vector x by p ... x = x[ p[:] ]
template <typename Vector1, typename Vector2>
void
permute (const Vector1& p, Vector2& x);

/*
// invert the permutation: y[ p[:] ] = x
template <typename Vector1, typename Vector2, typename vector3>
void
inverse_permute (const Vector1& p, const Vector2& x, Vector3& y);
*/

} // namespace spblas
} // namespace splib

#endif
