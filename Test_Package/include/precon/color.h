#ifndef __has_precon_color_h
#define __has_precon_color_h

#include <splib/csr_matrix.h>
#include <splib/vector.h>

namespace precon
{
namespace point_relax
{
   enum { seq_threshold = 50 };

   // Special block CSR matrix in split (L + D + U) format.
   template <typename ValueType>
   struct ColoredMatrix
   {
      enum { ValueTag = splib::ValueTags::TYPE< ValueType >::tag };
      enum { FormatTag = splib::FormatTags::SPARSE };

      typedef ValueType					value_type;
      typedef int					index_type;

      typedef splib::csr_matrix<value_type>		matrix_type;
      typedef splib::Vector<value_type>			vector_type;
      typedef splib::Vector<index_type>			permute_type;

      //const Matrix *_A;
      const matrix_type *_A;
      int num_colors;
      int max_block_size;
      //std::vector<int> offsets;
      splib::Vector<int> offsets;
      matrix_type *L, *U;
      vector_type dinv;
      permute_type	p, pinv; // permutation (and inverse permutation) vector
      //std::vector<int> is_seq; // should this set be executed sequentially (lex order)
      splib::Vector<int> is_seq; // should this set be executed sequentially (lex order)
      bool is_sparsified;

      ColoredMatrix(void);
      ~ColoredMatrix();
   };

} // namespace point_relax
} // namespace precon

#endif
