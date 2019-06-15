#ifndef __has_formats_h
#define __has_formats_h

#include <string>

namespace splib
{

struct row_major {};
struct column_major {};

struct matrix_format {};

struct sparse_format {};
struct dense_format {};

struct coo_format : public sparse_format, public matrix_format {};
struct csr_format : public sparse_format, public matrix_format {};
struct dia_format : public sparse_format, public matrix_format {};
struct ell_format : public sparse_format, public matrix_format {};
struct hyb_format : public sparse_format, public matrix_format {};

struct vector_format : public dense_format {};

struct partitioned_matrix_format : public sparse_format {};
//struct partitioned_vector_format : public vector_format {};
struct partitioned_vector_format : public dense_format {};

namespace FormatTags { enum Formats {CSR = 1, COO, DIA, ELL, HYB, SPARSE, DENSE, VECTOR}; };
namespace ValueTags
{
   enum Values {BOOL = 1, INT, FLOAT, DOUBLE};
   template <typename T> struct TYPE {};
   template <> struct TYPE<bool>    { enum { tag = BOOL	}; };
   template <> struct TYPE<int>     { enum { tag = INT	}; };
   template <> struct TYPE<float>   { enum { tag = FLOAT	}; };
   template <> struct TYPE<double>  { enum { tag = DOUBLE	}; };
}

std::string format_name (FormatTags::Formats);
std::string format_name (csr_format);
std::string format_name (coo_format);
std::string format_name (dia_format);
std::string format_name (ell_format);
std::string format_name (hyb_format);
std::string format_name (vector_format);
std::string format_name (dense_format);
std::string format_name (sparse_format);

int format_tag (csr_format);
int format_tag (coo_format);
int format_tag (dia_format);
int format_tag (ell_format);
int format_tag (hyb_format);
int format_tag (vector_format);
int format_tag (dense_format);
int format_tag (sparse_format);

} // splib

#endif
