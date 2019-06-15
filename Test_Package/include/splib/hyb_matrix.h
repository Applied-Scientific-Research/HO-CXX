#ifndef __has_hyb_matrix_h
#define __has_hyb_matrix_h

#include <splib/base_matrix.h>
#include <splib/formats.h>
#include <splib/vector.h>
#include <parallel/allocator.h>

#include <splib/ell_matrix.h>
#include <splib/coo_matrix.h>

namespace splib
{

template <typename ValueType,
          typename Allocator = parallel::aligned_allocator<ValueType> >
struct hyb_matrix :
   public BaseMatrix
{
   enum { FormatTag = FormatTags::HYB };
   enum { ValueTag = ValueTags::TYPE<ValueType>::tag };

   typedef BaseMatrix			Parent;

   typedef hyb_matrix<ValueType>	self_type;
   typedef Allocator			allocator_type;

   typedef splib::hyb_format	format_type;
   typedef ValueType		value_type;
   typedef Parent::index_type	index_type;

   typedef ell_matrix<ValueType>	ell_matrix_type;
   typedef coo_matrix<ValueType>	coo_matrix_type;

   ell_matrix_type	ell;
   coo_matrix_type	coo;

   template <typename OtherValueType,
             typename OtherAllocator = typename Allocator::template rebind<OtherValueType>::other>
   struct rebind
   {
      typedef hyb_matrix<OtherValueType,OtherAllocator> other;
   };

   // Construct empty matrix
   hyb_matrix (void);
   hyb_matrix (const index_type _num_rows,
               const index_type _num_columns,
               const index_type _num_values_coo,
               const index_type _num_values_per_row_ell);

   // Reallocate the internal storage ...
   void resize (const index_type _num_rows,
                const index_type _num_columns,
                const index_type _num_values_coo,
                const index_type _num_values_per_row_ell);

   int format_tag (void) const;
   int value_tag (void) const;
   void sort(void);
};

} // splib

#endif
