#ifndef __has_coo_matrix_h
#define __has_coo_matrix_h

#include <splib/formats.h>
#include <splib/base_matrix.h>
#include <splib/vector.h>
#include <parallel/allocator.h>
#include <parallel/numa.h>

namespace splib
{

// Sparse Coordinate (COO) format

template <typename ValueType,
          typename Allocator = parallel::aligned_allocator<ValueType> >
struct coo_matrix :
   //public splib::BaseMatrix<ValueType,splib::coo_format,IndexType>
   public splib::BaseMatrix
{
   enum { FormatTag = FormatTags::COO };
   enum { ValueTag = ValueTags::TYPE<ValueType>::tag };

   //typedef BaseMatrix<ValueType,splib::coo_format,IndexType> Parent;
   typedef BaseMatrix	Parent;

   //typedef coo_matrix<ValueType,Allocator,IndexType>	self_type;
   typedef coo_matrix<ValueType,Allocator>	self_type;
   typedef Allocator				allocator_type;

   //typedef typename Parent::format_type		format_type;
   //typedef typename Parent::value_type		value_type;
   typedef splib::coo_format		format_type;
   typedef ValueType			value_type;
   typedef typename Parent::index_type	index_type;

   template <typename OtherValueType,
             typename OtherAllocator = typename Allocator::template rebind<OtherValueType>::other >
   struct rebind
   {
      typedef coo_matrix<OtherValueType,OtherAllocator> other;
   };

   // Sparse storage vectors.
   splib::Vector<value_type,allocator_type>	values;
   splib::Vector<index_type,typename allocator_type::template rebind<index_type>::other >	colidx;
   splib::Vector<index_type,typename allocator_type::template rebind<index_type>::other >	rowidx;

   parallel::ParallelPropertiesType parallel_properties;

   coo_matrix (void);

   // Construct CSR with existing data ...
   template <typename IndexVector, typename ValueVector>
   coo_matrix (const index_type _num_rows,
               const index_type _num_columns,
               const index_type _num_values,
                     IndexVector& _rowidx, IndexVector& _colidx, ValueVector& _values);

   void resize (const index_type _num_rows, const index_type _num_columns, const index_type _num_values);

   int format_tag (void) const;
   int value_tag (void) const;
   //std::string name (void) const;
   void sort(void);

}; // coo_matrix

} // splib

#endif
