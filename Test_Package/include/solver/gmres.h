#ifndef __solver_gmres_h
#define __solver_gmres_h

#include <solver/solvers.h>
#include <solver/base_solver.h>
#include <splib/base_matrix.h>
#include <splib/vector.h>
#include <splib/dense_matrix.h>
#include <precon/identity.h>

namespace solver
{

template <typename Matrix>
struct GMRES
   : public BaseSolver
{
   //enum { SolverTag = SolverTags::GMRES };

   typedef Matrix			matrix_type;
   typedef typename Matrix::value_type	value_type;
   typedef GMRES<Matrix>		solver_type;

   typedef BaseSolver			Parent;

#ifdef USE_LIBNUMA
//#error 'GMRES not numa ready yet'
#warning 'GMRES not numa ready yet'
#endif
   typedef splib::Vector<value_type> vector_type;
   typedef vector_type interleaved_vector_type;

   typedef typename interleaved_vector_type::allocator_type _allocator;
   typedef typename _allocator::template rebind<float>::other _mixed_allocator;
   typedef splib::Vector<float, _mixed_allocator> mixed_interleaved_vector_type;

   typedef splib::DenseMatrix<value_type> dense_matrix_type;

   const Matrix *A_ptr;
   const precon::BasePrecon *baseM;
   vector_type r, w;
   //interleaved_vector_type p, z;
   //vector_type p, z;

#define __gmres_use_vector_of_vectors

   // Create a bunch of vectors ... these persist through the whole process.
#ifdef __gmres_use_vector_of_vectors
#warning 'using vector-of-vectors in fgmres'
   std::vector<vector_type> v, z;
#else
   dense_matrix_type v,z;
#endif
   vector_type cs, sn, s;
   dense_matrix_type h; // Hessenberg upper triangular matrix.

   int restart;

   GMRES (void);
   ~GMRES ();

   std::string name (void) const;
   int tag (void) const;
   int value_tag (void) const;
   int format_tag (void) const;

   int build (const splib::BaseMatrix *baseA);
   int build (const Matrix &A);

   //...int solve (splib::Vector<double> &x, const splib::Vector<double> &b, /*const*/ precon::BasePrecon *M);
   //int solve (splib::Vector<float> &x, const splib::Vector<float> &b, /*const*/ precon::BasePrecon *M);
   //template <typename Vector1, typename Vector2>
   //int solve (const Matrix &A, Vector1 &x, const Vector2 &b, precon::BasePrecon *baseM);
   int solve (const Matrix &A, vector_type &x, const vector_type &b, precon::BasePrecon *baseM);
};

} // namespace solver

#endif
