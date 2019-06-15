#ifndef __precon_identity_h
#define __precon_identity_h

#include <precon/precons.h>
#include <precon/base_precon.h>
#include <splib/vector.h>

#include <string>

namespace precon
{

struct Identity
   : public BasePrecon
{
   enum { PreconTag = PreconTags::Identity };

   typedef BasePrecon	Parent;

   Identity (void);
   Identity (const BasePreconOptions &options);
   ~Identity ();

   std::string name (const bool full_name = true) const;
   int tag (void) const;
   int format_tag (void) const;
   int value_tag (void) const;

   template <typename Matrix>
   int build (const Matrix& A);

   int build (const splib::BaseMatrix *baseA);

   template <typename Vector1, typename Vector2>
   int solve (const Vector1& x, Vector2& z) /*const*/;
   //template <typename Vector>
   //int solve (const Vector& x, Vector& z) /*const*/;

   //int solve (const splib::Vector<double>& x, splib::Vector<double>& z) /*const*/;
   //int solve (const splib::Vector<float >& x, splib::Vector<float >& z) /*const*/;
   //int solve (splib::BaseVector *basex, splib::BaseVector *basez) /*const*/;
};

} // namespace precon

#endif
