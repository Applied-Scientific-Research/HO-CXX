#include <cstdio>
#include <cstdlib>

#include <iostream>

#include "linsolvers/linsolver.hpp"

// Putting these in decending order based on efficiency / availability.

#ifdef ENABLE_HYPRE
# include "hypre.hpp"
# ifndef DEFAULT_LINEAR_SOLVER
#  define DEFAULT_LINEAR_SOLVER HYPRE
# endif
#endif

#ifdef ENABLE_APLLES
# include "aplles.hpp"
# ifndef DEFAULT_LINEAR_SOLVER
#  define DEFAULT_LINEAR_SOLVER APLLES
# endif
#endif

#ifdef ENABLE_EIGEN
# include "eigen.hpp"
# ifndef DEFAULT_LINEAR_SOLVER
#  define DEFAULT_LINEAR_SOLVER EIGEN
# endif
#endif

#ifdef ENABLE_AMGCL
# include "amgcl.hpp"
# ifndef DEFAULT_LINEAR_SOLVER
#  define DEFAULT_LINEAR_SOLVER AMGCL
# endif
#endif

namespace HighOrderFEM
{
namespace LinearSolver
{

namespace details
{

LinearSolverTag::Tag default_tag = LinearSolverTag::DEFAULT_LINEAR_SOLVER;

} // namespace

std::string LinearSolverTag::name( void ) const
   {
      switch( this->tag )
      {
         case APLLES: return "APLLES";
         case AMGCL : return "AMGCL";
         case EIGEN : return "EIGEN";
         case HYPRE : return "HYPRE";
         case AMGX  : return "AMGX";
      }
   }

LinearSolverTag::Tag getLinearSolverTag( const int id )
   {
      switch( id )
      {
         case (int)LinearSolverTag::APLLES: return LinearSolverTag::APLLES;
         case (int)LinearSolverTag::AMGCL : return LinearSolverTag::AMGCL;
         case (int)LinearSolverTag::EIGEN : return LinearSolverTag::EIGEN;
         case (int)LinearSolverTag::HYPRE : return LinearSolverTag::HYPRE;
         case (int)LinearSolverTag::AMGX  : return LinearSolverTag::AMGX;
      }

      // Error.
      fprintf(stderr,"Unknown solver id %d ... shutting down.\n", id);
      exit(1);
   }

std::shared_ptr<BaseLinearSolver>
  CreateLinearSolver( LinearSolverTag::Tag tag )
{
   switch( tag )
   {
#ifdef ENABLE_APLLES
      case ( LinearSolverTag::APLLES ) :
      {
         auto sptr = std::make_shared< ApllesSolver >();
         return std::dynamic_pointer_cast< BaseLinearSolver >( sptr );
      }
#endif
#ifdef ENABLE_HYPRE
      case ( LinearSolverTag::HYPRE ) :
      {
         auto sptr = std::make_shared< HypreSolver >();
         return std::dynamic_pointer_cast< BaseLinearSolver >( sptr );
      }
#endif
#ifdef ENABLE_AMGCL
      case ( LinearSolverTag::AMGCL ) :
      {
         auto sptr = std::make_shared< AMGCLSolver >();
         return std::dynamic_pointer_cast< BaseLinearSolver >( sptr );
      }
#endif
#ifdef ENABLE_EIGEN
      case ( LinearSolverTag::EIGEN ) :
      {
         auto sptr = std::make_shared< EigenSolver >();
         return std::dynamic_pointer_cast< BaseLinearSolver >( sptr );
      }
#endif
      default:
      {
         fprintf(stderr,"Requested linear solver %d not available. Reverting to default: %d\n", (int)tag, (int)details::default_tag );
         return CreateLinearSolver( details::default_tag );
      }
   }
}

std::shared_ptr<BaseLinearSolver>
  CreateLinearSolver( void )
{
   return CreateLinearSolver( details::default_tag );
}

std::vector< std::shared_ptr<BaseLinearSolver> >
  CreateAllLinearSolvers(void)
{
   std::vector< std::shared_ptr<BaseLinearSolver> > solvers;

#ifdef ENABLE_APLLES
      {
         auto sptr = std::make_shared< ApllesSolver >();
         solvers.push_back( std::dynamic_pointer_cast< BaseLinearSolver >( sptr ) );
      }
#endif
#ifdef ENABLE_HYPRE
      {
         auto sptr = std::make_shared< HypreSolver >();
         solvers.push_back( std::dynamic_pointer_cast< BaseLinearSolver >( sptr ) );
      }
#endif
#ifdef ENABLE_AMGCL
      {
         auto sptr = std::make_shared< AMGCLSolver >();
         solvers.push_back( std::dynamic_pointer_cast< BaseLinearSolver >( sptr ) );
      }
#endif
#ifdef ENABLE_EIGEN
      {
         auto sptr = std::make_shared< EigenSolver >();
         solvers.push_back( std::dynamic_pointer_cast< BaseLinearSolver >( sptr ) );
      }
#endif

   return solvers;
}

} // LinearSolver
} // HO
