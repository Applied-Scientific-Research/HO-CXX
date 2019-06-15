#ifndef __sort_h
#define __sort_h

namespace splib
{

// Vector-only functions ...

// Default sort ... ascending order
template <typename Vector>
void sort (Vector &x);

// sort keys[] and shuffle values[] at the same time.
template <typename Vector1, typename Vector2>
void sort_by_key (Vector1 &keys, Vector2 &values);

// sort keys[] and shuffle values1 and values2[] at the same time.
template <typename Vector1, typename Vector2, typename Vector3>
void sort_by_key (Vector1 &keys, Vector2 &values2, Vector3 &values3);

} // namespace splib

//#include <splib/details/sort.inl>

#endif
