#ifndef __has_csr_matrix_h
#define __has_csr_matrix_h

#include <splib/base_matrix.h>
#include <splib/formats.h>
#include <splib/vector.h>
#include <parallel/allocator.h>
#include <parallel/numa.h>

namespace splib
{

// Compressed Sparse Row (CSR) format

template <typename ValueType,
          typename Allocator = parallel::aligned_allocator<ValueType> >
struct csr_matrix :
   public BaseMatrix
{
   enum { FormatTag = FormatTags::CSR };
   enum { ValueTag = ValueTags::TYPE<ValueType>::tag };

   //enum { ChunkSize = 4 };
   //enum { ChunkSize = 8 };
   //enum { ChunkSize = __LEVEL1_DCACHE_LINESIZE / sizeof(ValueType) };
   enum { ChunkSize = 1 };

   //typedef BaseMatrix<ValueType,splib::csr_format,IndexType> Parent;
   typedef BaseMatrix Parent;

   //typedef csr_matrix<ValueType,Allocator,IndexType>	self_type;
   typedef csr_matrix<ValueType>	self_type;
   typedef Allocator					allocator_type;

   typedef splib::csr_format	format_type;
   typedef ValueType		value_type;
   typedef Parent::index_type	index_type;

   // Sparse storage vectors.
   Vector<value_type,allocator_type>	values;
   Vector<index_type,typename allocator_type::template rebind<index_type>::other >	colidx;
   Vector<index_type,typename allocator_type::template rebind<index_type>::other >	rowptr;

#ifdef _OPENMP
//#define __EnableNonuniformRowPartitioning
#ifdef  __EnableNonuniformRowPartitioning
   Vector<index_type> row_partition;
#endif
#endif

   parallel::ParallelPropertiesType parallel_properties;

   template <typename OtherValueType,
             typename OtherAllocator = typename Allocator::template rebind<OtherValueType>::other>
   struct rebind
   {
      typedef csr_matrix<OtherValueType,OtherAllocator> other;
   };

   // Construct empty matrix
   csr_matrix (void);

   // Construct CSR with existing data ...
   template <typename IndexVector, typename ValueVector>
   csr_matrix (const index_type _num_rows,
               const index_type _num_columns,
               const index_type _num_values,
                     IndexVector& _rowptr, IndexVector& _colidx, ValueVector& _values);

   // Make a CSR view from existing data ...
   csr_matrix (const index_type _num_rows,
               const index_type _num_columns,
               const index_type _num_values,
                     index_type *_rowptr, index_type *_colidx, value_type *_values);

   // Reallocate the internal storage ...
   void resize (const index_type _num_rows, const index_type _num_columns, const index_type _num_values);

   // Reset a CSR view from existing data ...
   void set_view (const index_type _num_rows,
                  const index_type _num_columns,
                  const index_type _num_values,
                        index_type *_rowptr, index_type *_colidx, value_type *_values);

   int format_tag (void) const;
   int value_tag (void) const;
   //std::string name (void) const;
   void sort(void);
}; // csr_matrix

} // splib

#endif
