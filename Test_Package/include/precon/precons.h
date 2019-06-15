#ifndef __precons_h
#define __precons_h

#include <string>

namespace precon
{

namespace PreconTags
{
   enum Tag { Identity = 1, Jacobi, SGS, AMG, Krylov, Diagonal = Jacobi };
};

PreconTags::Tag precon_tag (int tag);
std::string precon_name (PreconTags::Tag tag);
std::string precon_name (int tag);

} // namespace precon

#endif
