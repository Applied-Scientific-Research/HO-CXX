#ifndef __has_ell_matrix_h
#define __has_ell_matrix_h

#include <splib/base_matrix.h>
#include <splib/formats.h>
#include <splib/vector.h>
#include <parallel/allocator.h>
#include <parallel/numa.h>

namespace splib
{

// ELLPACK (ELL) format

template <typename ValueType,
          typename Allocator = parallel::aligned_allocator<ValueType> >
struct ell_matrix :
   public BaseMatrix
{
   enum { FormatTag = FormatTags::ELL };
   enum { ValueTag = ValueTags::TYPE<ValueType>::tag };

   typedef BaseMatrix Parent;

   typedef ell_matrix<ValueType>	self_type;
   typedef Allocator					allocator_type;

   typedef splib::ell_format	format_type;
   typedef ValueType		value_type;
   typedef Parent::index_type	index_type;

   enum { ignore_invalid_index = 0 };
   const static index_type invalid_index = static_cast<index_type>(-1);

   parallel::ParallelPropertiesType parallel_properties;

   index_type num_values_per_row;
   index_type array_pitch;

   // Sparse storage vectors.
   typedef typename allocator_type::template rebind<index_type>::other index_allocator_type;
   Vector<value_type,       allocator_type > values;
   Vector<index_type, index_allocator_type > colidx;

   template <typename OtherValueType,
             typename OtherAllocator = typename Allocator::template rebind<OtherValueType>::other>
   struct rebind
   {
      typedef ell_matrix<OtherValueType,OtherAllocator> other;
   };

   // Construct empty matrix
   ell_matrix (void);
   ell_matrix (const index_type _num_rows,
               const index_type _num_columns,
               const index_type _num_values_per_row);

   // Reallocate the internal storage ...
   void resize (const index_type _num_rows,
                const index_type _num_columns,
                const index_type _num_values_per_row);

   int format_tag (void) const;
   int value_tag (void) const;
   void sort(void);

   inline index_type pitch(void) const
   {
      return (this->array_pitch) ? this->array_pitch : this->num_rows;
   }

   inline index_type index (const index_type &i, const index_type &j) const
   {
      return i + this->pitch() * j;
   }
};

} // splib

#endif
