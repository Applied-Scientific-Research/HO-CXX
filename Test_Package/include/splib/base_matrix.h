#ifndef __has_BaseMatrix_h
#define __has_BaseMatrix_h

//#include <splib/formats.h>
#include <splib/attributes.h>
//#include <parallel/numa.h>

namespace splib
{

// Base class
struct BaseMatrix
{
   typedef int	index_type;

   index_type	num_rows;
   index_type	num_columns;
   index_type	num_values;

   AttributeType attributes;
   //parallel::ParallelPropertiesType parallel_properties;

   BaseMatrix (void);
   BaseMatrix (const index_type _num_rows, const index_type _num_columns, const index_type _num_values);

   virtual ~BaseMatrix(void);

   // All matrices have a format and value-type tag for run-time checks.
   virtual int format_tag (void) const = 0;
   virtual int value_tag (void) const = 0;
   //virtual std::string name (void) const = 0;

   // All matrices have a sort method but it does different thangs for each format (possibily nothing at all) ...
   virtual void sort (void);
};

} // splib

   #define __decode_matrix_cast__(__newtype__, __pointer__) \
         dynamic_cast<__newtype__ *>(__pointer__);

   #define __decode_matrix_format__(__baseA, __value, __A, __func) \
   { \
      if      (__baseA->format_tag() == splib::FormatTags::CSR) \
      { \
         typedef splib::csr_matrix< __value > __matrix; \
         __matrix *__A = __decode_matrix_cast__(__matrix, __baseA); \
         __func; \
      } \
      else if (__baseA->format_tag() == splib::FormatTags::COO) \
      { \
         typedef splib::coo_matrix< __value > __matrix; \
         __matrix *__A = __decode_matrix_cast__(__matrix, __baseA); \
         __func; \
      } \
      else if (__baseA->format_tag() == splib::FormatTags::DIA) \
      { \
         typedef splib::dia_matrix< __value > __matrix; \
         __matrix *__A = __decode_matrix_cast__(__matrix, __baseA); \
         __func; \
      } \
      else if (__baseA->format_tag() == splib::FormatTags::ELL) \
      { \
         typedef splib::ell_matrix< __value > __matrix; \
         __matrix *__A = __decode_matrix_cast__(__matrix, __baseA); \
         __func; \
      } \
      else if (__baseA->format_tag() == splib::FormatTags::HYB) \
      { \
         typedef splib::hyb_matrix< __value > __matrix; \
         __matrix *__A = __decode_matrix_cast__(__matrix, __baseA); \
         __func; \
      } \
      else \
      { \
         fprintf(stderr,"Invalid A->format_tag() %d\n", __baseA->format_tag()); \
         exit(-1); \
      } \
   }
   #define __decode_matrix__(__baseA, __A, __func) \
   { \
      if      (__baseA->value_tag() == splib::ValueTags::DOUBLE) \
      { \
         __decode_matrix_format__(__baseA, double, __A, __func); \
      } \
      else if (__baseA->value_tag() == splib::ValueTags::FLOAT) \
      { \
         __decode_matrix_format__(__baseA, float, __A, __func); \
      } \
      else \
      { \
         fprintf(stderr,"Invalid A->value_tag() %d\n", __baseA->value_tag()); \
         exit(-1); \
      } \
   }

#endif
