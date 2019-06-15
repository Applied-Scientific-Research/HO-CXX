#ifndef __has_splib_dense_matrix_h
#define __has_splib_dense_matrix_h

#include <splib/formats.h>
#include <splib/base_matrix.h>
#include <splib/vector.h>

namespace splib
{

template <typename ValueType, typename Orientation = splib::column_major>
struct DenseMatrix
   : public BaseMatrix
{
   enum { FormatTag = FormatTags::DENSE };
   enum { ValueTag = ValueTags::TYPE<ValueType>::tag };

   typedef Orientation						orientation;
   typedef BaseMatrix						Parent;
   typedef DenseMatrix<ValueType>				self_type;
   typedef splib::Vector<ValueType>				vector_type;
   typedef splib::dense_format					format_type;
   //typedef IndexType						index_type;
   typedef int							index_type;
   typedef ValueType						value_type;

   typedef size_t	size_type;

   vector_type	values;

   DenseMatrix (void);
   DenseMatrix (index_type __m, index_type __n, const value_type __value = value_type(0));

   explicit
   DenseMatrix (index_type __m, index_type __n, const value_type *__p);

   template <typename MatrixType>
   DenseMatrix(MatrixType& A);

   void resize (const index_type& __m, const index_type& __n, const value_type __value = value_type());

   const value_type& operator() (const index_type& i, const index_type& j) const;
         value_type& operator() (const index_type& i, const index_type& j);

   int format_tag (void) const;
   int value_tag (void) const;

#if 0
   // Define column and row-vector view types
   typedef splib::StrideVector<value_type>	vector_view_type;

   vector_view_type _column_vector_view (index_type& j, splib::column_major) const
      {
         value_type *p = const_cast<value_type*>( &this->operator()(0,j) );
         return vector_view_type( p, this->num_rows, 1 );
      }
   vector_view_type _column_vector_view (index_type& j, splib::row_major) const
      {
         value_type *p = const_cast<value_type*>( &this->operator()(0,j) );
         return vector_view_type( p, this->num_rows, this->num_columns );
      }

   vector_view_type _row_vector_view (index_type i, splib::column_major) const
      {
         value_type *p = const_cast<value_type*>( &this->operator()(i,0) );
         return vector_view_type( p, this->num_columns, this->num_rows );
      }
   vector_view_type _row_vector_view (index_type i, splib::row_major) const
      {
         value_type *p = const_cast<value_type*>( &this->operator()(i,0) );
         return vector_view_type( p, this->num_columns, 1 );
      }


   vector_view_type column_vector_view (index_type j) const
      { return this->_column_vector_view (j,orientation()); }
   vector_view_type row_vector_view (index_type i) const
      { return this->_row_vector_view (i,orientation()); }
#endif
};

} // namespace splib


#endif
