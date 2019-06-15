#ifndef __has_vector_h
#define __has_vector_h

#include <parallel/allocator.h>
#include <splib/formats.h>

namespace splib
{

struct BaseVector
{
   //typedef splib::vector_format	format_type;

   //virtual ~BaseVector(void) = 0;

   //virtual size_t size (void) const;
   //virtual void* getPointer (void) const;
   //virtual void set_view (void *ptr, size_t n);
   virtual size_t size(void) const = 0;
   virtual bool empty(void) const = 0;
   virtual int value_tag(void) const = 0;
};

template <typename T, typename Allocator = parallel::aligned_allocator<T> >
struct Vector
   : public BaseVector
{
   typedef splib::vector_format			format_type;
   typedef Allocator				allocator_type;
   typedef BaseVector				Parent;

   enum { ValueTag = ValueTags::TYPE<T>::tag };

   // Take from parent ... this should be the allocator.
   typedef typename allocator_type::value_type		value_type;
   typedef typename allocator_type::size_type		size_type;
   typedef typename allocator_type::difference_type	difference_type;
   typedef typename allocator_type::pointer		pointer;
   typedef typename allocator_type::const_pointer	const_pointer;
   typedef typename allocator_type::reference		reference;
   typedef typename allocator_type::const_reference	const_reference;

   // Use pointer as iterator
   typedef pointer		iterator;
   typedef const_pointer	const_iterator;

   allocator_type	m_allocator;
   pointer		m_pointer; // pointer to the first element.
   size_type		m_size;
   bool			m_allocated;

   AttributeType	attributes;
   parallel::ParallelPropertiesType parallel_properties;

   template <typename OtherT,
             typename OtherAlloc = typename allocator_type::template rebind<OtherT>::other >
   struct rebind
   {
      typedef Vector<OtherT,OtherAlloc> other;
   };

   Vector (void);

   // Create a new array of length n and initialize with __value (or the default)
   explicit
   Vector (size_type __n, const value_type& __value = value_type());
   //explicit
   //Vector (const int __n, const value_type& __value = value_type());

   // Create a new array from existing memory and length n ... no allocation and no initialization.
   explicit
   Vector (const pointer __p, size_type __n);

   // Create a new array that's a copy of __a ... new memory is formed!
   //explicit
   Vector (const Vector& __a);
   //Vector (Vector& __a);

   // Create an array from an iterator range.
   template <typename It>
   Vector (It __begin, It __end);

   ~Vector (void);

   // ... use the default copy constructor
   template <typename _OtherVectorType>
   Vector& operator= (const _OtherVectorType& __a);
   //Vector& operator= (_OtherVectorType& __a);

   // Vector inquires
   size_type size(void) const;
   size_type max_size(void) const;
   bool empty(void) const;
   int value_tag(void) const;

   // Adjust size ...
   void resize (size_type __n, value_type __value = value_type());

   void clear (void);

   // Explicitly set the internal data.
   //explicit
   void set_view (const pointer __p, size_type __n);

   pointer get_pointer (void) const;

   // STL information
   iterator       begin()       { return iterator(this->m_pointer); };
   const_iterator begin() const { return const_iterator(this->m_pointer); };
   iterator       end()         { return iterator(this->m_pointer + this->m_size); };
   const_iterator end()   const { return const_iterator(this->m_pointer + this->m_size); };

   // Random access
   inline reference operator[](size_type i) const { return m_pointer[i]; };
   inline reference operator[](size_type i)       { return m_pointer[i]; };
   inline reference operator[](int i) const { return m_pointer[i]; };
   inline reference operator[](int i)       { return m_pointer[i]; };
};

   #define __decode_vector__(__baseObject, __Object, __func) \
   { \
      if      (__baseObject->value_tag() == ValueTags::DOUBLE) \
      { \
         typedef Vector< double > __vector; \
         __vector *__Object = dynamic_cast<__vector *>(__baseObject); \
         __func; \
      } \
      else if (__baseObject->value_tag() == ValueTags::FLOAT) \
      { \
         typedef Vector< float > __vector; \
         __vector *__Object = dynamic_cast<__vector *>(__baseObject); \
         __func; \
      } \
      else \
      { \
         fprintf(stderr,"Invalid Object->value_tag() %d\n", __baseObject->value_tag()); \
         exit(-1); \
      } \
   }

} // splib

#endif
