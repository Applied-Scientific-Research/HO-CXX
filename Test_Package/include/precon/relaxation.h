#ifndef __has_precon_relaxation_h
#define __has_precon_relaxation_h

#include <splib/formats.h>
#include <splib/csr_matrix.h>
#include <splib/ell_matrix.h>

#include <precon/precons.h>
#include <precon/base_precon.h>
#include <precon/color.h>

#include <utils/type_traits.h>

#include <string>

namespace precon
{

template <typename Matrix, PreconTags::Tag Tag = PreconTags::Jacobi>
struct Relaxation
   : public BasePrecon
{
   enum { PreconTag = Tag };
   enum { FormatTag = Matrix::FormatTag };
   enum { ValueTag = Matrix::ValueTag };

   typedef BasePrecon			Parent;

   // Test for valid relaxation enum ...
   typedef typename utils::enable_if< PreconTag == (int)PreconTags::Jacobi ||
                                      PreconTag == (int)PreconTags::SGS, bool>::type __type;

   typedef Matrix				matrix_type;

   typedef typename matrix_type::value_type	value_type;
   typedef typename matrix_type::index_type	index_type;

   typedef splib::Vector<value_type>		vector_type;

   typedef splib::csr_matrix<value_type>	csr_matrix;

   value_type		omega;

   matrix_type*		A_ptr;	// May be a copy or a view.
   bool			A_is_allocated;
   csr_matrix		A_csr;	// Only need a CSR of the input matrix for SGS.

   vector_type		dinv;
   vector_type		r;
   vector_type		x_mixed, y_mixed;

   typedef point_relax::ColoredMatrix<value_type> ColoredMatrix;

   ColoredMatrix	colored_matrix;
   bool			color_matrix;
   value_type		truncation_rate;
   int			max_colors;
   int			color_threshold;

   Relaxation (void);
   Relaxation (const BasePreconOptions &options);

   ~Relaxation();

   std::string name (const bool full_name = true) const;
   int tag (void) const;
   int value_tag (void) const;
   int format_tag (void) const;

   //int build (const matrix_type &A);
   template <typename OtherMatrix>
   int build (const OtherMatrix &A);

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
