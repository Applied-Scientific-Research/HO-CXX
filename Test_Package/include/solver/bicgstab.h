#ifndef __solver_bicgstab_h
#define __solver_bicgstab_h

#include <solver/solvers.h>
#include <solver/base_solver.h>
#include <splib/base_matrix.h>
#include <splib/vector.h>
#include <precon/identity.h>

namespace solver
{

template <typename Matrix>//, typename Precon = precon::Identity >
struct BiCGstab
   : public BaseSolver
{
   enum { SolverTag = SolverTags::BiCGstab };
   enum { ValueTag  = Matrix::ValueTag };
   enum { FormatTag = Matrix::FormatTag };

   typedef Matrix			matrix_type;
   typedef typename Matrix::value_type	value_type;
   typedef BiCGstab<Matrix>		solver_type;

   typedef BaseSolver			Parent;

#ifdef USE_LIBNUMA
//#error 'BiCGstab not numa ready'
#warning 'BiCGstab not numa ready'
//   typedef splib::Vector<value_type,parallel::page_aligned_allocator<value_type> > vector_type;
//   typedef splib::Vector<value_type,parallel::numa_interleaved_allocator<value_type> > interleaved_vector_type;
//#else
#endif
   typedef splib::Vector<value_type> vector_type;
//   typedef vector_type interleaved_vector_type;
//#endif

//   typedef typename interleaved_vector_type::allocator_type _allocator;
//   typedef typename _allocator::template rebind<float>::other _mixed_allocator;
//   typedef splib::Vector<float, _mixed_allocator> mixed_interleaved_vector_type;

   const Matrix *A_ptr;
   const precon::BasePrecon *baseM;
   //vector_type r, p, q, s, t, z;
   vector_type r, p, q, s;// z;
   vector_type phat, shat;
   vector_type rtilde;
   //interleaved_vector_type p, z;

   BiCGstab (void);
   ~BiCGstab ();

   std::string name (void) const;
   int tag (void) const;
   int value_tag (void) const;
   int format_tag (void) const;

   int build (const splib::BaseMatrix *baseA);
   int build (const Matrix &A);

   //template <typename Vector1, typename Vector2>
   //int solve (const Matrix &A, Vector1 &x, const Vector2 &b, precon::BasePrecon *M);
   int solve (const Matrix &A, vector_type &x, const vector_type &b, precon::BasePrecon *M);
   //...int solve (splib::Vector<double> &x, const splib::Vector<double> &b, /*const*/ precon::BasePrecon *M);
};

} // namespace solver

#endif
