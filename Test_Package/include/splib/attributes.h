#ifndef __has_attributes_h
#define __has_attributes_h

#include <string>

// Matrix attributes

namespace splib
{

namespace AttributeTags
{
   // Bits for attributes ...
   enum AttributeTag
        { Symmetric	= 0,
          Pattern	= 1,
          Sparse	= 2,
          Sorted	= 4,
          NUMAed	= 5 };

   //Attribute attr_get_enum (int tag);
   std::string attr_name (AttributeTag tag);
}

// Binary attributes ... on or off only.
class AttributeType
{
private:
   unsigned int value;

public:
   AttributeType(void);

   bool is_enabled (AttributeTags::AttributeTag attr) const;
   void enable (AttributeTags::AttributeTag attr);
   void disable (AttributeTags::AttributeTag attr);

   std::string to_string();
};

bool attr_is_enabled (const AttributeType &obj, AttributeTags::AttributeTag attr);
AttributeType attr_enable (const AttributeType &obj, AttributeTags::AttributeTag attr);
AttributeType attr_disable (const AttributeType &obj, AttributeTags::AttributeTag attr);

} // splib

#endif
