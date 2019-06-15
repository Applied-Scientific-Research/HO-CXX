#ifndef __has_splib_copy_h
#define __has_splib_copy_h

namespace splib
{

template <typename Matrix1, typename Matrix2>
void matrix_copy (const Matrix1 &A, Matrix2 &B);

template <typename Vector1, typename Vector2>
void vector_copy (const Vector1 &x, Vector2 &y);

} // namespace splib

#endif
