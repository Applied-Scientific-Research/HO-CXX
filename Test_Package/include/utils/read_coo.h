#ifndef __utils_read_coo_h
#define __utils_read_coo_h

#include <splib/coo_matrix.h>

namespace utils
{

template <typename Matrix>
int read_coo (const char *filename, Matrix& _A)
   {
      typedef typename Matrix::value_type value_type;
      typedef typename Matrix::index_type index_type;

      typedef splib::coo_matrix<index_type,value_type> coo_matrix;

      coo_matrix A;

      FILE *f = fopen(filename,"r");
      if (f == NULL)
      {
         fprintf(stderr,"error opening %s", filename);
         return 1;
      }

      int num_rows, num_columns, num_values;

      if (fread(&num_rows, sizeof(int), 1, f) != 1)
      {
         fprintf(stderr,"error reading num_rows");
         return 1;
      }
      if (fread(&num_columns, sizeof(int), 1, f) != 1)
      {
         fprintf(stderr,"error reading num_columns");
         return 1;
      }
      if (fread(&num_values, sizeof(int), 1, f) != 1)
      {
         fprintf(stderr,"error reading num_values");
         return 1;
      }

      if (__DEBUG) printf("%s: num_rows,num_columns,num_values=%d,%d,%d\n", filename, num_rows, num_columns, num_values);

      A.resize(num_rows,num_columns,num_values);

      if (fread(&A.rowidx[0], sizeof(int), num_values, f) != num_values)
      {
         fprintf(stderr,"error reading rowidx");
         return 1;
      }
      if (fread(&A.colidx[0], sizeof(int), num_values, f) != num_values)
      {
         fprintf(stderr,"error reading colidx");
         return 1;
      }
      if (fread(&A.values[0], sizeof(double), num_values, f) != num_values)
      {
         fprintf(stderr,"error reading values");
         return 1;
      }

      fclose(f);

      splib::copy(A,_A);

      return 0;
   }

} // utils

#endif
