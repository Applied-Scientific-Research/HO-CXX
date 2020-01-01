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

   //!< Maximum # of iterations allowed by the iterative solver.
   int maxiters;

   /*! Target tolerance normalized by the RHS norm:
    * @f[
    * \frac{|| Ax - f ||}{||f||} < relTol
    * @f]
    */
   double reltol;

   /*! Target absolute tolerance.
    * @f[
    * || Ax - f || < absTol
    * @f]
    */
   double abstol;

   //!< Krylov subspace dimension for GMRES (etc.)
   int restart_k;

   //!< Control the level of reporting. 0=none, 1=minimal, 2=very.
   int verbosity;

   index_type num_rows;

   BaseLinearSolver( LinearSolverTag::Tag tag )
      : abstol(1e-10), reltol(0), maxiters(200), verbosity(0), restart_k(16),
        tag(tag),
        num_rows(0)
   {
      std::cout << "BaseLinearSolver: " << this << std::endl;
   }

   ~BaseLinearSolver()
   {
      std::cout << "~BaseLinearSolver: " << this << std::endl;
   }

   void setMaxIterations(const int max_iters) { this->maxiters = max_iters; }
   int  getMaxIterations(void) const { return this->maxiters; }

   void setRelativeTolerance(const double reltol) { this->reltol = reltol; }
   void setAbsoluteTolerance(const double abstol) { this->abstol = abstol; }

   double getRelativeTolerance(void) const { return this->reltol; }
   double getAbsoluteTolerance(void) const { return this->abstol; }

   void setVerbosity(const int level) { this->verbosity = level; }
   int  getVerbosity(void) const { return this->verbosity; }

   virtual int build( const index_type num_rows, const std::vector<index_type>& rowptr,
                                                 const std::vector<index_type>& colidx,
                                                 const std::vector<value_type>& values ) = 0;

   virtual int solve ( const std::vector<value_type>& f, std::vector<value_type>& u ) = 0;
   //virtual int solve ( const std::vector<value_type>& f, std::vector<value_type>& u, std::vector<value_type>& u0 ) = 0;
};

} // LinearSolver
} // HO

#ifdef ENABLE_EIGEN
# include "eigen.hpp"
#endif

#ifdef ENABLE_AMGCL
# include "amgcl.hpp"
#endif

#ifdef ENABLE_HYPRE
# include "hypre.hpp"
#endif

#ifdef ENABLE_APLLES
# include "aplles.hpp"
#endif

#endif
