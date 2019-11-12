#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <typeinfo>

extern "C"
{
#include "globheads.h"
#include "defs.h"
#include "protos.h"
#include "ios.h"

void set_arms_pars(io_t* io,  int Dscale, int *ipar, double *tolcoef, int *lfil);

}

namespace details
{

   typedef enum : int {
      ILUK = 1,
      ARMS2,
      ILUT,
      BILUK
   }
   PreconTag;

   PreconTag getPreconTag (int precon)
   {
      switch(precon)
      {
         case(ILUK ): return ILUK;
         case(ARMS2): return ARMS2;
         case(ILUT ): return ILUT;
         case(BILUK): return BILUK;
         default:
         {
            fprintf(stderr,"Unknown tag %d\n", precon);
            exit(2);
         }
      }
   }

   io_t init_io(size_t ndim, size_t nnz)
   {
      int lfil_arr[7];
      double droptol[7], dropcoef[7];
      int ipar[18];
      int diagscal = 1;
      double tolind = 0.7;

      io_t io;
      memset( &io, 0, sizeof(io) ); // clear all terms.

      /*-------------------- set parameters for arms */
      io.nparam   = 1;
      io.im       = 60;
      io.maxits   = 200;
      io.tol      = 1e-8;
      io.lfil0    = 50;
      io.lfilInc  = 1;
      io.tol0     = 0.001;
      io.tolMul   = 0.01;
      io.fill_lev = 1;
      io.perm_type= 0;
      io.Bsize    = 30;
      io.nnz      = nnz;
      io.ndim     = ndim;

      return io;
   }

   SPreptr create_arms2( csptr spmat, io_t& io )
   {
      armsMat *arms_p = new armsMat;
      setup_arms( arms_p );

      int lfil_arr[7];
      double droptol[7], dropcoef[7];
      int ipar[18];

      int diagscal = 1;
      double tolind = 0.7;

      /*-------------------- set parameters for arms */
   
      set_arms_pars( &io, diagscal, ipar, dropcoef, lfil_arr );

      int ierr = arms2( spmat, ipar, droptol, lfil_arr, tolind, arms_p, stderr );
      if (ierr != 0)
      {
         fprintf(stderr,"ARMS2 error %d\n", ierr);
         exit(1);
      }

      SPre *P = new SPre;

      P->ARMS = arms_p;
      P->precon = preconARMS;

      //printf("built ARMS2 P %x\n", P);

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

      return P;
   }
   SPreptr create_vbiluk( csptr spmat, int lfil )
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

      auto *lu = new VBILUSpar;

      printf("starting vbiluk(%d)\n", lfil);
      int ierr = vbilukC( lfil, vbmat, lu, stderr );
      if (ierr != 0)
      {
         fprintf(stderr,"VBILUK error %d\n", ierr);
         exit(1);
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

      printf("built VBILUK %x fill-factor= %f permuted=%d\n", P, fill, permuted);

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
               const int mrows, const int nnz,
               int A_rowptr[], int A_colidx[], double A_values[],
               void **v_P )
{
   using namespace details;

   PreconTag tag = getPreconTag(precon);

   printf("Inside itsol_create %d %d %d %d %d %f\n", precon, tag, mrows, nnz, lfil, droptol);

   int ierr = 0;

   SparMat *spmat = new SparMat;
   if ( (ierr = CSRcs( mrows, A_values, A_colidx, A_rowptr, spmat, 0 ) ) != 0)
   {
      fprintf(stderr,"Error building SparMat from CSR format %d\n", ierr);
      return ierr;
   }

   //printf("after CSRcs\n");

   auto io = init_io( mrows, nnz );

   /*-------------------- set parameters for arms */
   io.nparam   = 1;
   io.im       = 60;
   io.maxits   = 200;
   io.tol      = 1e-8;
   io.lfil0    = 50;
   io.lfilInc  = 1;
   io.tol0     = 0.001;
   io.tolMul   = 0.01;
   io.fill_lev = lfil;
   io.perm_type= 0;
   io.Bsize    = 30;

   SPreptr P = NULL;
   switch(tag)
   {
      case(ARMS2) : { P = create_arms2( spmat, io ); break; }
      case(ILUK)  : { P = create_iluk( spmat, lfil ); break; }
      case(ILUT)  : { P = create_ilut( spmat, lfil, droptol ); break; }
      case(BILUK) : { P = create_vbiluk( spmat, lfil ); break; }
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
      if (P->VBILU) cleanVBILU( P->VBILU );
      delete P;
      return 0;
   }
   else
      return 1;
}

int itsol_solve( double *b, double *x, SPre *P )
{
   //printf("solve P %x\n", P);
   return P->precon( b, x, P );
}

} // extern-C
