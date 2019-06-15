#ifndef __has_dia_matrix_h
#define __has_dia_matrix_h

#include <splib/base_matrix.h>
#include <splib/formats.h>
#include <splib/vector.h>
#include <parallel/allocator.h>
#include <parallel/numa.h>

namespace splib
{

// DIAgonal format

template <typename ValueType,
          typename Allocator = parallel::aligned_allocator<ValueType> >
struct dia_matrix :
   public BaseMatrix
{
   enum { FormatTag = FormatTags::DIA };
   enum { ValueTag = ValueTags::TYPE<ValueType>::tag };

   typedef BaseMatrix				Parent;

   typedef dia_matrix<ValueType,Allocator>	self_type;
   typedef Allocator				allocator_type;

   //typedef typename Parent::format_type	format_type;
   //typedef typename Parent::value_type	value_type;
   typedef splib::dia_format			format_type;
   typedef ValueType				value_type;
   typedef typename Parent::index_type		index_type;

   template <typename OtherValueType,
             typename OtherAllocator = typename Allocator::template rebind<OtherValueType>::other >
   struct rebind
   {
      typedef dia_matrix<OtherValueType,OtherAllocator> other;
   };

   typedef splib::Vector<value_type,allocator_type>	vector_type;

   // Sparse storage vectors.
   int				num_diagonals;
   splib::Vector<int>		offsets;
   //splib::Vector<vector_type>	diagonals;
   std::vector<vector_type>	diagonals;
   //splib::Vector<value_type>	values;
   vector_type			values;

   parallel::ParallelPropertiesType parallel_properties;

   dia_matrix (void);

   // Make a CSR view from existing data ...
   dia_matrix (const index_type _num_rows,
               const index_type _num_columns,
               const index_type _num_diagonals,
                     index_type *_offsets,
                     value_type *_values);

   void
   resize (const index_type _num_rows, const index_type _num_columns, const index_type _num_diagonals);

   // Reset view from existing data ...
   void
   set_view (const index_type _num_rows,
             const index_type _num_columns,
             const index_type _num_diagonals,
                 index_type *_offsets, value_type *_values);

   int format_tag (void) const;
   int value_tag (void) const;
   //std::string name (void) const;
   //void sort (void);
};

} // splib

#endif
