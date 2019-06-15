#ifndef __solver_pcg_h
#define __solver_pcg_h

#include <solver/solvers.h>
#include <solver/base_solver.h>
#include <splib/base_matrix.h>
#include <splib/vector.h>
#include <precon/identity.h>

namespace solver
{

template <typename Matrix>//, typename Precon = precon::Identity >
//template <typename ValueType>//, typename Precon = precon::Identity >
struct PCG
   : public BaseSolver
{
   //enum { SolverTag = SolverTags::PCG };
//   enum { PreconTag = Precon::PreconTag };

   typedef Matrix			matrix_type;
   typedef typename Matrix::value_type	value_type;
   //typedef ValueType			value_type;
//   typedef PCG<Matrix,Precon>		solver_type;
   typedef PCG<Matrix>			solver_type;
   //typedef PCG<ValueType>		solver_type;
//   typedef Precon			precon_type;

   typedef BaseSolver			Parent;

   enum { ValueTag = Matrix::ValueTag };
   enum { FormatTag = Matrix::FormatTag };

   //typedef Parent::BaseMatrix		BaseMatrix;
   //typedef Parent::BaseVector		BaseVector;
   //typedef Parent::BasePrecon		BasePrecon;

   typedef splib::Vector<value_type> vector_type;
   typedef vector_type interleaved_vector_type;
#if 0
#ifdef USE_LIBNUMA
   //typedef splib::Vector<value_type,parallel::page_aligned_allocator<value_type> > vector_type;
   typedef splib::Vector<value_type,parallel::numa_interleaved_allocator<value_type> > interleaved_vector_type;
#else
   //typedef splib::Vector<value_type> vector_type;
   typedef vector_type interleaved_vector_type;
#endif
#endif

   typedef typename interleaved_vector_type::allocator_type _allocator;
   typedef typename _allocator::template rebind<float>::other _mixed_allocator;
   typedef splib::Vector<float, _mixed_allocator> mixed_interleaved_vector_type;

   //const Matrix *A_ptr;
   //const Precon *M_ptr;
   //const precon::BasePrecon *baseM;
   vector_type r, q;
   interleaved_vector_type p, z;
   //mixed_interleaved_vector_type z_mixed;

   PCG (void);
   ~PCG ();

   std::string name (void) const;
   int tag (void) const;
   int value_tag (void) const;
   int format_tag (void) const;
   //int precon_tag (void) const;

   int build (const splib::BaseMatrix *baseA);
   //template <typename OtherMatrix> int build (const OtherMatrix &A);
   int build (const Matrix &A);

/*   int setPrecon (const precon::BasePrecon *baseM);
   template <typename Precon> int setPrecon (const Precon &M);

   precon::BasePrecon* getPrecon (void) const;*/

   //template <typename Matrix, typename Vector1, typename Vector2, typename Precon>
   //template <typename Vector1, typename Vector2, typename Precon>
   //int solve (const Matrix& A, Vector1 &x, const Vector2& b, const Precon& M);
   //template <typename Vector1, typename Vector2, typename Precon>
   //int solve (Vector1 &x, const Vector2& b, const Precon &M);

   //int solve (const BaseMatrix *A, BaseVector *x, const BaseVector *b, const BasePrecon *M);
   //int solve (splib::BaseVector *x, const splib::BaseVector *b);

   //template <typename Matrix, typename Vector1, typename Vector2>
   //template <typename Vector1, typename Vector2>
   //int solve (const Matrix& A, Vector1 &x, const Vector2& b, precon::BasePrecon *baseM);
   int solve (const Matrix& A, vector_type &x, const vector_type& b, precon::BasePrecon *baseM);
   //int solve (splib::BaseMatrix *baseA,
   //           splib::BaseVector  *basex,
   //           splib::BaseVector  *baseb,
   //           precon::BasePrecon *baseM = NULL);

   //...int solve (splib::Vector<double> &x, const splib::Vector<double> &b, /*const*/ precon::BasePrecon *M);
};

} // namespace solver

#endif
