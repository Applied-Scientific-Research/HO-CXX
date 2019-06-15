#ifndef __has_precon_krylov_h
#define __has_precon_krylov_h

#include <splib/formats.h>

#include <precon/precons.h>
#include <precon/base_precon.h>
#include <precon/relaxation.h>

#include <solver/gmres.h>
#include <solver/pcg.h>

#include <string>

namespace precon
{

template <typename Matrix>
struct Krylov
   : public BasePrecon
{
   enum { PreconTag = PreconTags::Krylov };
   enum { FormatTag = Matrix::FormatTag };
   enum { ValueTag = Matrix::ValueTag };

   typedef BasePrecon			Parent;

   typedef solver::GMRES<Matrix>	SolverType;
   //typedef solver::PCG<Matrix>	SolverType;
   typedef precon::Relaxation<Matrix>	InnerPreconType;
   //typedef precon::Relaxation<Matrix,PreconTags::SGS>	InnerPreconType;

   typedef Matrix				matrix_type;

   typedef typename matrix_type::value_type	value_type;
   typedef typename matrix_type::index_type	index_type;

   //splib::BaseMatrix	*baseA;
   matrix_type		*A_ptr;
   bool			A_is_allocated;
   //solver::BaseSolver	*base_solver;
   SolverType		solver;
   InnerPreconType	precon;

   Krylov (void);
   Krylov (const BasePreconOptions &options);
   ~Krylov();

   std::string name (const bool full_name = true) const;
   int tag (void) const;
   int value_tag (void) const;
   int format_tag (void) const;

   //int build (const matrix_type &A);
   template <typename OtherMatrix>
   int build (const OtherMatrix &A);
   //int build (const Matrix &A);

   int build (const splib::BaseMatrix *baseA);

   // Apply N relaxation steps. y = M^(-1)x
   template <typename Vector1, typename Vector2>
   int solve (const Vector1& x, Vector2& y) ;//const;

   //int solve (const splib::Vector<value_type>& x, splib::Vector<value_type>& y) ;//const;
   //int solve (const splib::Vector<double>& x, splib::Vector<double>& y) ;//const;
   //int solve (const splib::Vector<value_type>& x, splib::Vector<value_type>& y) ;//const;
   //int solve (const splib::Vector<double>& x, splib::Vector<float >& y) ;//const;
   //int solve (splib::BaseVector *basex, splib::BaseVector *basey) ;//const;
};

} // namespace precon

#endif
