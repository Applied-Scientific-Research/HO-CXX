#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <typeinfo>

#include <vector>
#include <map>
#include <algorithm>
#include <memory>

extern "C"
{
#include "globheads.h"
#include "defs.h"
#include "protos.h"

}

namespace details
{

   struct IndependentSets
   {
      std::vector<int> level_ptr, row_sets;

      IndependentSets()
      {
         printf("Created IndependentSets %x\n", this);
      }
      ~IndependentSets()
      {
         printf("Destroyed IndependentSets %x\n", this);
      }
   };
   struct BlockLUIndependentSets
   {
      IndependentSets L, U;
      int bs;

      BlockLUIndependentSets() : bs(0)
      {
         printf("Created BlockLUIndependentSets %x\n", this);
      }
      ~BlockLUIndependentSets()
      {
         printf("Destroyed BlockLUIndependentSets %x\n", this);
      }

      // y := alpha Ax + beta y
      template <typename T>
      void gemv( const int m, const T alpha, const T _A[], const T x[], const T beta, T y[] )
      {
#define A(i,j) ( _A[ (i) + (j)*m ] ) // column-major
          for (int i = 0; i < m; ++i)
          {
              T yi = 0.0;

              #pragma ivdep
              for (int j = 0; j < m; ++j)
                  yi += A(i,j) * x[j];

              if (beta == T(0))
                  y[i] = alpha * yi;
              else
                  y[i] = alpha * yi + beta * y[i];
          }
#undef A

          return;
      }

#ifndef _OPENMP
      template <int _bs>
      int _solve_bs( double *b, double *x, VBILUSpar *LU, const int bs_rt = 0 )
      {
         const int bs = (_bs > 0) ? _bs : bs_rt;
         const int n = LU->n;
         double xtmp[bs];

         vbsptr L = LU->L;
         vbsptr U = LU->U;
         BData *D = LU->D;

         /* Block L solve */
         for (int i = 0; i < n; i++)
         {
            int nBs = bs * i;
            double *x_i = x + nBs;
            for (int k = 0; k < bs; k++)
               x_i[k] = b[nBs+k];

            int nzcount = L->nzcount[i];
            int *ja = L->ja[i];
            double **ba = L->ba[i];
            for (int j = 0; j < nzcount; j++)
            {
               int jcol = ja[j];
               double *data = ba[j];
               double *x_j = x + jcol * bs;
               gemv( bs, -1.0, data, x_j, 1.0, x_i ); 
            }
         }

         /* Block -- U solve */
         for (int i = n-1; i >= 0; i--)
         {
            int nzcount = U->nzcount[i];
            double *x_i = x + i*bs;

            int *ja = U->ja[i];
            double **ba = U->ba[i];

            for (int j = 0; j < nzcount; j++ )
            {
               int jcol = ja[j];
               double *x_j = x + jcol*bs;
               double *data = ba[j];
               gemv( bs, -1.0, data, x_j, 1.0, x_i ); 
            }

            double *data = D[i];
            gemv( bs, 1.0, data, x_i, 0.0, xtmp );

            for (int k = 0; k < bs; ++k)
               x_i[k] = xtmp[k];
         }

         return 0;
      }
#else
      template <int _bs>
      int _solve_bs( double *b, double *x, VBILUSpar *LU, const int bs_rt = 0 )
      {
         const int min_rows = 5;

         const int bs = (_bs > 0) ? _bs : bs_rt;
         const int n = LU->n;

         vbsptr L = LU->L;
         vbsptr U = LU->U;
         BData *D = LU->D;

         #pragma omp parallel
         {
            double xtmp[bs];

         /* Block L solve */
         const int num_levels_L = this->L.level_ptr.size()-1;

         auto forward_op = [&](const int i) {
               int nBs = bs * i;
               double *x_i = x + nBs;
               for (int k = 0; k < bs; k++)
                  x_i[k] = b[nBs+k];

               int nzcount = L->nzcount[i];
               int *ja = L->ja[i];
               double **ba = L->ba[i];
               for (int j = 0; j < nzcount; j++)
               {
                  int jcol = ja[j];
                  double *data = ba[j];
                  double *x_j = x + jcol * bs;
                  gemv( bs, -1.0, data, x_j, 1.0, x_i ); 
               }
            };

         for (int level = 0; level < num_levels_L; ++level)
         {
            const int nrows = this->L.level_ptr[level+1] - this->L.level_ptr[level];
            if ( nrows == min_rows )
            {
               #pragma omp master
               for (int ii = this->L.level_ptr[level]; ii < this->L.level_ptr[level+1]; ++ii)
                  forward_op( this->L.row_sets[ ii ] );
            }
            else
            {
               #pragma omp barrier

            #pragma omp for
            for (int ii = this->L.level_ptr[level]; ii < this->L.level_ptr[level+1]; ++ii)
            {
               const int i = this->L.row_sets[ii];

               forward_op(i);
               /*int nBs = bs * i;
               double *x_i = x + nBs;
               for (int k = 0; k < bs; k++)
                  x_i[k] = b[nBs+k];

               int nzcount = L->nzcount[i];
               int *ja = L->ja[i];
               double **ba = L->ba[i];
               for (int j = 0; j < nzcount; j++)
               {
                  int jcol = ja[j];
                  double *data = ba[j];
                  double *x_j = x + jcol * bs;
                  gemv( bs, -1.0, data, x_j, 1.0, x_i ); 
               }*/
            }

            }
         }

         /* Block -- U solve */
         const int num_levels_U = this->U.level_ptr.size()-1;

         auto backward_op = [&]( const int i ) {
               int nzcount = U->nzcount[i];
               double *x_i = x + i*bs;

               int *ja = U->ja[i];
               double **ba = U->ba[i];

               for (int j = 0; j < nzcount; j++ )
               {
                  int jcol = ja[j];
                  double *x_j = x + jcol*bs;
                  double *data = ba[j];
                  gemv( bs, -1.0, data, x_j, 1.0, x_i ); 
               }

               double *data = D[i];
               gemv( bs, 1.0, data, x_i, 0.0, xtmp );

               for (int k = 0; k < bs; ++k)
                  x_i[k] = xtmp[k];
            };

         for (int level = 0; level < num_levels_U; ++level)
         {
            const int nrows = this->U.level_ptr[level+1] - this->U.level_ptr[level];
            if ( nrows < min_rows )
            {
               #pragma omp master
               for (int ii = this->U.level_ptr[level]; ii < this->U.level_ptr[level+1]; ++ii)
                  backward_op( this->U.row_sets[ii] );
            }
            else
            {

            #pragma omp barrier

            #pragma omp for
            for (int ii = this->U.level_ptr[level]; ii < this->U.level_ptr[level+1]; ++ii)
            {
               const int i = this->U.row_sets[ii];

               backward_op(i);
               /*int nzcount = U->nzcount[i];
               double *x_i = x + i*bs;

               int *ja = U->ja[i];
               double **ba = U->ba[i];

               for (int j = 0; j < nzcount; j++ )
               {
                  int jcol = ja[j];
                  double *x_j = x + jcol*bs;
                  double *data = ba[j];
                  gemv( bs, -1.0, data, x_j, 1.0, x_i ); 
               }

               double *data = D[i];
               gemv( bs, 1.0, data, x_i, 0.0, xtmp );

               for (int k = 0; k < bs; ++k)
                  x_i[k] = xtmp[k];*/
            }

            }
         }

         } // end parallel

         return 0;
      }
#endif

      int solve_constant_bs( double *b, double *x, VBILUSpar *LU )
      {
         if (this->bs == 9)
            return this->_solve_bs<9>( b, x, LU );
         else
            return this->_solve_bs<0>( b, x, LU, this->bs );
      }

      int solve( double *b, double *x, VBILUSpar *LU )
      {
         const int n = LU->n;
         const int* bsz = LU->bsz;

         double alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0;
         int inc = 1;
         char* transA = "n";

         vbsptr L = LU->L;
         vbsptr U = LU->U;
         BData *D = LU->D;

         /* Block L solve */
         for (int i = 0; i < n; i++)
         {
            int dim = B_DIM(bsz,i);
            int nBs = bsz[i];
            for (int j = 0; j < dim; j++)
               x[nBs+j] = b[nBs+j];

            int nzcount = L->nzcount[i];
            int *ja = L->ja[i];
            double **ba = L->ba[i];
            for (int j = 0; j < nzcount; j++)
            {
               int icol = ja[j];
               int sz = B_DIM(bsz,icol);
               double *data = ba[j];
               DGEMV( transA, dim, sz, alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc ); 
            }
         }

         /* Block -- U solve */
#if 0
         for (int i = n-1; i >= 0; i--)
         {
            int dim = B_DIM(bsz,i);
            int nzcount = U->nzcount[i];
            int nBs = bsz[i];
            int *ja = U->ja[i];
            double **ba = U->ba[i];

            for (int j = 0; j < nzcount; j++ )
            {
               int icol = ja[j];
               int sz = B_DIM(bsz,icol);
               double *data = ba[j];
               DGEMV( transA, dim, sz, alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc ); 
            }

            double *data = D[i];
	    DGEMV( transA, dim, dim, alpha2, data, dim, x+nBs, inc, beta2, LU->bf, inc ); 

            for (int bi = 0; bi < dim; bi++)
               x[nBs+bi] = LU->bf[bi];
         }
#else
         const int num_levels_U = this->U.level_ptr.size()-1;
         const int bs = 9;

         #pragma omp parallel
         for (int level = 0; level < num_levels_U; level++)
         {
            //printf("%d, %d, %d\n", level, this->U.level_ptr[level], this->U.level_ptr[level+1]);

            #pragma omp for
            for (int ii = this->U.level_ptr[level]; ii < this->U.level_ptr[level+1]; ++ii)
            {
               const int i = this->U.row_sets[ii];

               int dim = B_DIM(bsz,i);
               int nzcount = U->nzcount[i];
               int nBs = bsz[i];
               int *ja = U->ja[i];
               double **ba = U->ba[i];

               for (int j = 0; j < nzcount; j++ )
               {
                  int icol = ja[j];
                  int sz = B_DIM(bsz,icol);
                  double *data = ba[j];
                  DGEMV( transA, dim, sz, alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc ); 
               }

               double *data = D[i];
               double tmp[bs];
	       //DGEMV( transA, dim, dim, alpha2, data, dim, x+nBs, inc, beta2, LU->bf, inc ); 
	       DGEMV( transA, dim, dim, alpha2, data, dim, x+nBs, inc, beta2, tmp, inc ); 

               for (int bi = 0; bi < dim; bi++)
                  //x[nBs+bi] = LU->bf[bi];
                  x[nBs+bi] = tmp[bi];
            }
         }
#endif

         return 0;
      }
   };

   std::map< void *, std::shared_ptr< BlockLUIndependentSets > > blocklu_ptr_map;

   typedef enum : int
   {
      ILUK = 1,
      ARMS2= 2,
      ILUT = 3,
      ILUC = 4,
      BILUK= 5,
      BILUT= 6,
   }
   PreconTag;

   PreconTag getPreconTag (int precon)
   {
      switch(precon)
      {
         case(ILUK ): return ILUK;
         case(ARMS2): return ARMS2;
         case(ILUT ): return ILUT;
         case(ILUC ): return ILUC;
         case(BILUK): return BILUK;
         case(BILUT): return BILUT;
         default:
         {
            fprintf(stderr,"Unknown tag %d\n", precon);
            exit(2);
         }
      }
   }

   int find_independet_sets ( VBSparMat* A, std::vector<int>& level_ptr, std::vector<int>& row_sets )
   {
      const int n = A->n;
      const int nnz0 = nnzVBMat(A);
      const int nnz = nnz0 + n;
      //printf("n, nnz, nnz0 %d %d %d\n", n, nnz, nnz0);

      // Convert to normal CSR struct ...
      // ... and ... create a CSC struct ...
      std::vector<int> rowptr(n+1), colidx(nnz);
      std::vector<int> colcnt(n,0);
      std::vector<int> colptr(n+1), rowidx(nnz);

      rowptr[0] = 0;
      for (int i = 0; i < n; ++i)
      {
         rowptr[i+1] = rowptr[i] + ( A->nzcount[i]+1 );

         colcnt[i]++;
         colidx[ rowptr[i] ] = i;
         for (int k = 0; k < A->nzcount[i]; ++k)
         {
            int col = A->ja[i][k];
            colcnt[ col ]++;
            colidx[ rowptr[i]+(k+1) ] = col;
         }
      }

      struct coo_pair {
         int row,col;
         coo_pair(int i, int j) : row(i), col(j) {}
      };

      //std::vector< coo_pair > coo;

      //for (int i = 0; i < n; ++i)
      //{
      //   rowptr[i+1] = rowptr[i] + ( A->nzcount[i]+1 );

      //   coo.push_back( coo_pair(i,i) );

      //   for (int k = 0; k < A->nzcount[i]; ++k)
      //   {
      //      int j = A->ja[i][k];
      //      coo.push_back( coo_pair(i,j) );
      //   }
      //}

      //printf("coo: %d\n", coo.size());

      //for (int i = 0; i < 10; ++i)
      //   for (int k = rowptr[i]; k < rowptr[i+1]; ++k)
      //      printf("row %d, %d\n", i, colidx[k]);

      //for (int i = 0; i < 10; ++i)
      //   printf("colcnt[%d]= %d\n", i, colcnt[i]);

      colptr[0] = 0;
      for (int i = 0; i < n; i++)
      {
         colptr[i+1] = colptr[i] + colcnt[i];
         colcnt[i] = 0;
      }
      //printf("colptr[n]= %d\n", colptr[n]);

      //for (int i = 0; i < 10; ++i)
      //   printf("colptr[%d]= %d %d\n", i, colptr[i], colptr[i+1]-colptr[i]);

      for (int i = 0; i < n; ++i)
      {
         for (int k = rowptr[i]; k < rowptr[i+1]; ++k)
         {
            int col = colidx[k];
            int kk = colptr[col] + colcnt[col];
            rowidx[ kk ] = i;
            ++colcnt[col];
         }
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

         //for (int i = first; i < min(n,first+10); ++i)
         //   printf("before: %d, %d\n", row_sets[i], rownnz[row_sets[i]]);

         std::stable_sort( &row_sets[first], row_sets.data() + n, comp_sort );

         //for (int i = first; i < min(n,first+10); ++i)
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
            for (int i = first; i < min(n,first+10); ++i)
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

      std::vector<int> level_sizes(nlevels), level_index(nlevels);
      for (int i = 0; i < nlevels; i++) {
         level_sizes[i] = level_ptr[i+1] - level_ptr[i];
         level_index[i] = i;
      }

      auto sort_sizes = [&] ( const int& i, const int& j ) { return level_sizes[i] > level_sizes[j]; };

      std::stable_sort( level_index.begin(), level_index.end(), sort_sizes );

      //for (int i = 0; i < min(20,level_ptr.size()-1); i++) {
      for (int i = 0; i < min(20,nlevels); i++)
      {
         int level = level_index[i];
         int len = level_ptr[level+1]-level_ptr[level];
         printf("level: %d count: %d\n", level, len);
         printf("[");
         len = min(20,len);
         for (int j = 0; j < len; ++j)
            printf("%d%s", row_sets[level_ptr[level]+j], (j != len-1) ? ", " : "]\n");
      }
   }

   SPreptr create_arms2( csptr spmat, int lfil0, double tol0, int max_levels = 20)
   {
      armsMat *arms_p = new armsMat;
      setup_arms( arms_p );

      int lfil[7];
      double droptol[7];
      int ipar[18];

      int diagscal = 1;
      double tolind = 0.05; //0.7;

      const char *names[] = { "Max levels",
                        "Reordering",
                        "Bsize",
                        "Verbose" };

      for (int i = 0; i < 18; ++i)
         ipar[i] = 0;

      ipar[0] = max_levels; // max levels to consider (hardwired at 10 max)

      //ipar[1] = 0; // level reordering
      //ipar[1] = 1; // ddPQ
      ipar[1] = 2;

      ipar[2] = 9; //Bsize ... smallest size allowed fr last schur comp.
      ipar[3] = 1; // verbose

      /*-------------------- interlevel methods */  
      ipar[10] = 0;       /* Always do permutations - currently not used  */    
      ipar[11] = 0;       /* ILUT or ILUTP - currently only ILUT is implemented */
      ipar[12] = diagscal;  /* diagonal row scaling before PILUT 0:no 1:yes */
      ipar[13] = diagscal;  /* diagonal column scaling before PILUT 0:no 1:yes */
      /*-------------------- last level methods */  
      ipar[14] = 1;       /* Always do permutations at last level */
      ipar[15] = 1;       /* ILUTP for last level(0 = ILUT at last level) */
      ipar[16] = diagscal;  /* diagonal row scaling  0:no 1:yes */
      ipar[17] = diagscal;  /* diagonal column scaling  0:no 1:yes */

      /*--------- dropcoef (droptol[k] = tol0*dropcoef[k]) ----- */
      droptol[0] = 1.6;     /* dropcoef for L of B block */
      droptol[1] = 1.6;     /* dropcoef for U of B block */
      droptol[2] = 1.6;     /* dropcoef for L\ inv F */
      droptol[3] = 1.6;     /* dropcoef for U\ inv E */
      droptol[4] = 0.004;   /* dropcoef for forming schur comple. */
      droptol[5] = 0.004;   /* dropcoef for last level L */
      droptol[6] = 0.004;   /* dropcoef for last level U */

      droptol[4] = droptol[5] = droptol[6] = 1./1.6;
      droptol[4] = droptol[5] = droptol[6] = 1.;

      for (int i = 0; i < 4; i++)
         printf("%s : %d\n", names[i], ipar[i]);

      int nnz = nnz_cs(spmat);
      int n = spmat->n;

      for (int i = 0; i < 7; ++i)
      {
         lfil[i] = lfil0 * int(nnz/n);
         droptol[i] *= tol0;
         printf("%d %f %d %f\n", lfil[i], droptol[i], lfil0, tol0);
      }

      int ierr = arms2( spmat, ipar, droptol, lfil, tolind, arms_p, stderr );
      if (ierr != 0)
      {
         fprintf(stderr,"ARMS2 error %d\n", ierr);
         exit(1);
      }

      int lu_nnz = nnz_arms(arms_p, stderr);

      SPre *P = new SPre;

      P->ARMS = arms_p;
      P->precon = preconARMS;

      printf("built ARMS2 P %x fill %f\n", P, float(lu_nnz) / nnz);

      return P;
   }
   SPreptr create_iluk( csptr spmat, int lfil )
   {
      ILUSpar *lu = new ILUSpar;

      printf("starting iluk(%d)\n", lfil);
      int ierr = ilukC( lfil, spmat, lu, stderr );
      if (ierr != 0)
      {
         fprintf(stderr,"ILUK error %d\n", ierr);
         exit(1);
      }

      auto fill = (double) nnz_ilu(lu) / nnz_cs(spmat);

      SPre *P = new SPre;

      P->ILU = lu;
      P->precon = preconILU;

      printf("built ILUK %x fill-factor= %f\n", P, fill);

      condestLU( lu, stderr );

      return P;
   }
   SPreptr create_ilut( csptr spmat, int lfil, double tol )
   {
      ILUSpar *lu = new ILUSpar;

      printf("starting ilut(%d,%f)\n", lfil, tol);
      int ierr = ilut( spmat, lu, lfil, tol, stderr );
      if (ierr != 0)
      {
         fprintf(stderr,"ILUT error %d\n", ierr);
         exit(1);
      }

      auto fill = (double) nnz_ilu(lu) / nnz_cs(spmat);

      SPre *P = new SPre;

      P->ILU = lu;
      P->precon = preconILU;

      printf("built ILUT %x fill-factor= %f\n", P, fill);

      condestLU( lu, stderr );

      return P;
   }
   SPreptr create_iluc( csptr spmat, int lfil, double tol )
   {
      LDUmat *ldumat = new LDUmat; // for input
      {
         int ierr = CSClumC( spmat, ldumat, 0 );
         if (ierr) {
            fprintf(stderr,"CSClumC error %d\n", ierr);
            exit(1);
         }
      }

      ILUSpar *lu = new ILUSpar;

      const int dropmeth = 0; // standard method

      printf("starting iluc(%d,%f)\n", lfil, tol);
      int ierr = ilutc( ldumat, lu, lfil, tol, dropmeth, stderr );
      if (ierr != 0)
      {
         fprintf(stderr,"ILUTC error %d\n", ierr);
         exit(1);
      }

      auto fill = (double) nnz_ilu(lu) / nnz_cs(spmat);

      SPre *P = new SPre;

      P->ILU = lu;
      P->precon = preconILU;

      printf("built ILUC %x fill-factor= %f\n", P, fill);

      condestLU( lu, stderr );

      //cleanILU ( ldumat ); 
      //delete ldumat;

      return P;
   }
   SPreptr create_vbilu( csptr spmat, int lfil, double tol, const bool useLevelDropMethod = true)
   {
      int nBlocks;
      int *blockDim;
      int *perm;

      {
         double t_h, t_a;
         double eps = 0.8;

         int ierr = init_blocks( spmat, &nBlocks, &blockDim, &perm, eps, &t_h, &t_a );
         if (ierr != 0) {
            fprintf(stderr,"Error in init_blocks %d\n", ierr);
            exit(2);
         }

         printf("nBlocks %d\n", nBlocks);
      }

      auto *vbmat = new VBSparMat;

      {
         int ierr = csrvbsrC( 1, nBlocks, blockDim, spmat, vbmat );
         if (ierr != 0) {
            fprintf(stderr,"Error in csrvbsr %d\n", ierr);
            exit(2);
         }

         printf("block rows: %d %d\n", spmat->n, vbmat->n);
         auto nnz = nnz_cs(spmat);
         auto bnnz = nnzVBMat(vbmat);
         printf("block vals: %d %d %f\n", nnz, bnnz, (double) nnz/ bnnz);
      }

      // Test if all are the same blocksize.
      bool A_all_equal = true;
      int max_blk_sz = 0;
      for (int i = 0; i < vbmat->n; ++i) {
         A_all_equal &= ( B_DIM(vbmat->bsz, 0) == B_DIM(vbmat->bsz,i) );
         max_blk_sz = max( max_blk_sz, B_DIM(vbmat->bsz,i) );
      }
      printf("max_blk_sz: %d %d\n", A_all_equal, max_blk_sz);

      auto *lu = new VBILUSpar;

      if ( useLevelDropMethod )
      {
         printf("starting vbiluk(%d)\n", lfil);
         int ierr = vbilukC( lfil, vbmat, lu, stderr );
         if (ierr != 0)
         {
            fprintf(stderr,"VBILUK error %d\n", ierr);
            exit(1);
         }
      }
      else
      {
         std::vector<double*> wp( vbmat->n );
         std::vector<double> w( vbmat->n * max_blk_sz * max_blk_sz );
         for (int i = 0; i < vbmat->n; ++i)
            wp[i] = w.data() + i * max_blk_sz * max_blk_sz;

         printf("starting vbilut(%d,%.2e)\n", lfil, tol);
         int ierr = vbilutC( vbmat, lu, lfil, tol, wp.data(), stderr );
         if (ierr != 0)
         {
            fprintf(stderr,"VBILUT error %d\n", ierr);
            exit(1);
         }
      }

      auto fill = (double) nnz_vbilu(lu) / nnz_cs(spmat);

      SPre *P = new SPre;

      P->VBILU = lu;
      P->precon = preconVBR;

      bool permuted = false;
      for (int i = 0; i < spmat->n; ++i)
         if ( perm[i] != i ) {
            permuted = true;
            break;
         }

      bool LU_all_equal = true;
      for (int i = 0; i < lu->n; ++i)
         LU_all_equal &= ( B_DIM(lu->bsz, 0) == B_DIM(lu->bsz,i) );

      printf("built VBILUx %x fill-factor= %f permuted=%d diag= %d equal-blocks= %d %d %d\n", P, fill, permuted, lu->DiagOpt, A_all_equal, LU_all_equal, B_DIM(lu->bsz, 0));

      {
         auto ptr = std::make_shared< BlockLUIndependentSets >();
         blocklu_ptr_map[P] = ptr;

         find_independet_sets ( lu->L, ptr->L.level_ptr, ptr->L.row_sets );
         find_independet_sets ( lu->U, ptr->U.level_ptr, ptr->U.row_sets );

         if (LU_all_equal)
            ptr->bs = B_DIM(lu->bsz, 0);

         P->Ptype = 2;
      }

      VBcondestC( lu, stderr );

      cleanVBMat( vbmat );
      free( perm );
      free( blockDim );

      return P;
   }
}

extern "C"
{

int itsol_create(
               const int precon, // choice of precon
               const int lfil,
               const double droptol,
               const int max_levels,
               const int mrows, const int nnz,
               int A_rowptr[], int A_colidx[], double A_values[],
               void **v_P )
{
   using namespace details;

   PreconTag tag = getPreconTag(precon);

   printf("Inside itsol_create %d %d %d %d %d %f %d\n", precon, tag, mrows, nnz, lfil, droptol, max_levels);

   int ierr = 0;

   SparMat *spmat = new SparMat;
   if ( (ierr = CSRcs( mrows, A_values, A_colidx, A_rowptr, spmat, 0 ) ) != 0)
   {
      fprintf(stderr,"Error building SparMat from CSR format %d\n", ierr);
      return ierr;
   }

   //printf("after CSRcs\n");

   SPreptr P = NULL;
   switch(tag)
   {
      case(ARMS2) : { P = create_arms2( spmat, lfil, droptol, max_levels ); break; }
      case(ILUK)  : { P = create_iluk( spmat, lfil ); break; }
      case(ILUT)  : { P = create_ilut( spmat, lfil, droptol ); break; }
      case(ILUC)  : { P = create_iluc( spmat, lfil, droptol ); break; }
      case(BILUK) : { P = create_vbilu( spmat, lfil, droptol ); break; }
      case(BILUT) : { P = create_vbilu( spmat, lfil, droptol, false ); break; }
   }

   cleanCS( spmat );

   *v_P = (void*)P;
   return 0;
}

int itsol_destroy( SPre *P )
{
   //printf("destroy P %x\n", P);
   if (P)
   {
      if (P->ARMS) cleanARMS( P->ARMS );
      if (P->ILU) cleanILU( P->ILU );
      if (P->VBILU)
      {
         using details::blocklu_ptr_map;
         auto it = blocklu_ptr_map.find(P);
         if ( it != blocklu_ptr_map.end() )
            blocklu_ptr_map.erase( it );

         cleanVBILU( P->VBILU );
      }
      delete P;
      return 0;
   }
   else
      return 1;
}

int itsol_solve( double *b, double *x, SPre *P )
{
//   printf("solve P %x %d\n", P, P->Ptype);
/*   if (P->VBILU)
   {
      using details::blocklu_ptr_map;
      auto it = blocklu_ptr_map.find(P);
      if ( it != blocklu_ptr_map.end() )
         if ( it->second->bs )
            return it->second->solve_constant_bs( b, x, P->VBILU );
         else
            return it->second->solve( b, x, P->VBILU );
   }*/

   return P->precon( b, x, P );
}

} // extern-C
