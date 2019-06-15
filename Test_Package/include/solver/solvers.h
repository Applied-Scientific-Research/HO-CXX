#ifndef _solver_solvers_h
#define _solver_solvers_h

#include <string>

namespace solver
{

namespace SolverTags
{
   //enum Tag {PCG = 1, BiCGstab, GMRES, FGMRES, AMG };
   enum Tag {PCG = 1, BiCGstab, GMRES, RGMRES, FGMRES };
};

SolverTags::Tag solver_tag (int tag);
std::string solver_name (SolverTags::Tag tag);
std::string solver_name (int tag);

} // namespace solver

#endif
