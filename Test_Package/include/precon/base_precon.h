#ifndef __base_precon_h
#define __base_precon_h

#include <precon/precons.h>
//#include <utils/xml.h>
#include <splib/base_matrix.h>
#include <splib/vector.h>

#include <string>

namespace precon
{

// Run-time options for the BasePrecon
struct BasePreconOptions
{
   int max_iters;
   bool mixed_precision;

   BasePreconOptions (void);
   virtual ~BasePreconOptions();
};

struct BasePrecon
   : public BasePreconOptions
{
   PreconTags::Tag precon_tag;

   BasePrecon (PreconTags::Tag tag);
   BasePrecon (PreconTags::Tag tag, const BasePreconOptions& base_options);
   virtual ~BasePrecon();

   virtual std::string name (const bool full_name = true) const = 0;
   virtual int tag (void) const = 0;
   virtual int value_tag (void) const = 0;
   virtual int format_tag (void) const = 0;

   // Empty build method by default
   //template <typename Matrix>
   //int build (const Matrix& A);
   virtual int build (const splib::BaseMatrix *baseA) = 0;
   //virtual int build (const splib::BaseMatrix *baseA);

   //template <typename Vector1, typename Vector2>
   //virtual int solve (const Vector1& x, Vector2& z) /*const*/;

   //virtual int solve (const splib::Vector<double>& x, splib::Vector<double>& z) /*const*/ = 0;
   //virtual int solve (const splib::Vector<float>& x, splib::Vector<float >& z) /*const*/ = 0;
   //virtual int solve (const splib::BaseVector *basex, splib::BaseVector *basez) = 0;
   //virtual int solve (splib::BaseVector *basex, splib::BaseVector *basez) = 0;
   virtual int solve (const splib::BaseVector *basex, splib::BaseVector *basez);
};

//int solve (BasePrecon *baseM, splib::BaseVector *basex, splib::BaseVector *basez);

} // namespace precon

#endif
