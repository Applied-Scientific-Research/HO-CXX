#ifndef __linsolver_hpp
#define __linsolver_hpp

namespace HighOrderFEM
{
namespace LinearSolver
{

struct LinearSolverTag
{
   enum Tag : int
   {
      APLLES = 1,
      AMGCL,
      EIGEN,
      HYPRE,
      AMGX
   };

   Tag tag;
   int id;

   static int getTag( Tag etag )
   {
      switch( etag )
      {
         case APLLES: return (int) APLLES;
         case AMGCL : return (int) AMGCL;
         case EIGEN : return (int) EIGEN;
         case HYPRE : return (int) HYPRE;
         case AMGX  : return (int) AMGX;
      }

      // Not possible.
      return 0;
   }

   static Tag getTag( const int id )
   {
      switch( id )
      {
         case (int)APLLES: return APLLES;
         case (int)AMGCL : return AMGCL;
         case (int)EIGEN : return EIGEN;
         case (int)HYPRE : return HYPRE;
         case (int)AMGX  : return AMGX;
      }

      // Error.
      fprintf(stderr,"Unknown solver id %d ... shutting down.\n", id);
      exit(1);
   };

   LinearSolverTag( Tag tag ) : tag(tag), id( getTag(tag) ) {}
};

struct BaseLinearSolver
{
   typedef double value_type;
   typedef int    index_type;

   LinearSolverTag tag;

   double abstol, reltol;
   int maxiters;
   index_type num_rows;

   BaseLinearSolver( LinearSolverTag::Tag tag )
      : abstol(1e-10), reltol(0), maxiters(200),
        tag(tag),
        num_rows(0)
   {
      std::cout << "BaseLinearSolver: " << this << std::endl;
   }

   ~BaseLinearSolver()
   {
      std::cout << "~BaseLinearSolver: " << this << std::endl;
   }

   virtual int build( const index_type num_rows, index_type* rowptr, index_type* colidx, value_type* values ) = 0;
   virtual int solve ( const value_type* fp, value_type* up, const value_type* u0 = NULL ) = 0;
};

} // LinearSolver
} // HO

#ifdef ENABLE_EIGEN
# include "eigen.hpp"
#endif

#ifdef ENABLE_AMGCL
# include "amgcl.hpp"
#endif

#endif
