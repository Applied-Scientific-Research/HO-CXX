#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <typeinfo>

#include <vector>
#include <map>
#include <algorithm>
#include <memory>

namespace details
{

template <typename IndexType, typename Vector>
int find_independent_sets ( const IndexType n, const IndexType nnz,
                            IndexType rowptr[], IndexType colidx[],
                            Vector& level_ptr, Vector& row_sets )
{
   // create a CSC struct ...
   std::vector<int> colcnt(n,0);
   std::vector<int> colptr(n+1), rowidx(nnz);

   for (int i = 0; i < n; ++i)
      for (int k = rowptr[i]; k < rowptr[i+1]; ++k)
         colcnt[ colidx[k] ]++;

   //for (int i = 0; i < 10; ++i)
   //   printf("colcnt[%d]= %d\n", i, colcnt[i]);

   colptr[0] = 0;
   for (int i = 0; i < n; i++)
   {
      colptr[i+1] = colptr[i] + colcnt[i];
      colcnt[i] = 0;
   }

   for (int i = 0; i < n; ++i)
      for (int k = rowptr[i]; k < rowptr[i+1]; ++k)
      {
         int col = colidx[k];
         int kk = colptr[col] + colcnt[col];
         rowidx[ kk ] = i;
         ++colcnt[col];
      }

   //for (int j = 0; j < 10; ++j)
   //   for (int k = colptr[j]; k < colptr[j+1]; ++k)
   //      printf("col %d, %d\n", j, rowidx[k]);

   level_ptr.resize(1,0);
   row_sets.resize(n,0);

   std::vector<int> rownnz(n);

   for (int i = 0; i < n; ++i)
   {
      rownnz[i] = rowptr[i+1] - rowptr[i];
      row_sets[i] = i;
   }

   auto comp_sort = [&]( const int& i, const int& j )
      {
         //int row_i = row_sets[i];
         //int row_j = row_sets[j];
         //return rownnz[row_i] < rownnz[row_j];
         return rownnz[i] < rownnz[j];
      };

   auto comp_upper_bound = [&]( const int val, const int& item )
      {
         //int row = row_sets[item];
         return val < rownnz[item];
      };

   //for (int i = 0; i < 10; ++i)
   //   printf("%d, %d, %d\n", i, rowcnt[i].nnz, rowcnt[i].row);

   int rows_remaining = n;

   int level = 0;
   while (level < n and rows_remaining > 0)
   {
      int first = level_ptr[level];

      //for (int i = first; i < std::min(n,first+10); ++i)
      //   printf("before: %d, %d\n", row_sets[i], rownnz[row_sets[i]]);

      std::stable_sort( &row_sets[first], row_sets.data() + n, comp_sort );

      //for (int i = first; i < std::min(n,first+10); ++i)
      //   printf("top: %d, %d\n", row_sets[i], rownnz[row_sets[i]]);

      // Find the rows with the smallest count. ... 0 means a diagonal term.
      int next = -1;
      int nitems = 0;

      if ( rownnz[ row_sets[first] ] > 1 ) // this row isn't independent so process in order seq.
      {
         //next = first+1;
         fprintf(stderr,"No independent rows found at level %d\n", level);
         exit(1);
      }
      else if ( rownnz[ row_sets[first] ] != 1 ) // this row isn't independent so process in order seq.
      {
         //next = first+1;
         fprintf(stderr,"Unexpected dependency count at level %d\n", level);
         for (int i = first; i < std::min(n,first+10); ++i)
            printf("top: %d, %d\n", row_sets[i], rownnz[row_sets[i]]);
         exit(2);
      }
      else
      {
         // Find the length of this sequence.
         nitems = std::upper_bound( &row_sets[first], row_sets.data() + n, 1, comp_upper_bound ) - &row_sets[first];
         next = first + nitems;
      }

      //if ( next - first > 1 )
      //   printf("level %d, first %d, next %d, len %d\n", level, first, next, nitems);

      // Delete the dependence of these rows from the graph.
      for (int i = first; i < next; ++i)
      {
         // For each column on this row, break the edge.
         int row = row_sets[i];
         //printf("elem %d\n", row);
         if ( rownnz[row] > 1 ) { fprintf(stderr,"Error in nnz %d %d\n", row, rownnz[row]); exit(2); }

         const int col = row; // must be the diagonal block.
         for (int jk = colptr[col]; jk < colptr[col+1]; ++jk)
         {
            int other_row = rowidx[jk];
            //printf("col %d, %d, %d\n", col, other_row, rownnz[other_row]);
            --rownnz[other_row];
         }
      }

      level_ptr.push_back( next );

      rows_remaining -= nitems;
      level++;
   }

   const int nlevels = level_ptr.size()-1;
   printf("nlevels = %d\n", nlevels);

   for (int i = 0; i < std::min(20,nlevels); i++) {
      int len = level_ptr[i+1]-level_ptr[i];
      printf("level: %d count: %d\n", i, len);
      printf("[");
      len = std::min(20,len);
      for (int j = 0; j < len; ++j)
         printf("%d%s", row_sets[level_ptr[i]+j], (j != len-1) ? ", " : "]\n");
   }

   return nlevels;
}

} // details

extern "C"
{

int find_independent_sets(
               const int mrows, const int nnz,
               int rowptr[], int colidx[],
               int rowseq[] )
{
   std::vector<int> level_ptr, row_sets;
   int nlevels = details::find_independent_sets( mrows, nnz, rowptr, colidx, level_ptr, row_sets );

   for (int level = 0; level < nlevels; level++)
   {
      for (int k = level_ptr[level]; k < level; ++k)
      {
         const int row = row_sets[k];
         rowseq[row] = level;
      }
   }

   return nlevels;
}

} // extern-C
