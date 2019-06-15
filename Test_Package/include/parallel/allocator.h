#ifndef __has_parallel_allocator_h
#define __has_parallel_allocator_h

#include <parallel/parallel.h>
#include <parallel/numa.h>

#include <stddef.h>

#ifndef __LEVEL1_DCACHE_LINESIZE
#define __LEVEL1_DCACHE_LINESIZE (64)
#endif
#ifndef __PAGE_SIZE
#define __PAGE_SIZE (4096)
#endif

namespace parallel
{
   // Certain alignement sizes ...
   enum { __Align_m128	= 16 }; // sse/mmx
   enum { __Align_cache	= __LEVEL1_DCACHE_LINESIZE }; // probably 64
   enum { __Align_page	= __PAGE_SIZE }; // probably 4096

#ifdef USE_LIBNUMA
//#warning 'Default alignment = page'
   enum { __Align_default = __Align_page };
#else
//#warning 'Default alignment = cache'
   enum { __Align_default = __Align_cache };
#endif

   // Allocator for aligned data.
   //template <typename T, size_t __Alignment = sizeof(void*) >
   //template <typename T, size_t __Alignment = __Align_cache >
   template <typename T, size_t __Alignment = __Align_default >
   struct aligned_allocator
   {
      //public:

      enum { Alignment		= __Alignment };

      // The following will be the same for virtually all allocators.
      typedef T*		pointer;
      typedef const T*		const_pointer;
      typedef T&		reference;
      typedef const T&		const_reference;
      typedef T			value_type;
      typedef std::size_t	size_type;
      typedef ptrdiff_t		difference_type;

      // The following must be the same for all allocators.
      template <typename U>
      struct rebind
      {
         typedef aligned_allocator<U, Alignment> other;
      };

      T* address (T& r) const;
      const T* address (const T& s) const;

      size_t max_size() const;

      bool operator!= (const aligned_allocator& other) const;

      void construct (T* const p, const T& t) const;

      void destroy (T* const p) const;

      // Returns true iff storage allocated from *this
      // can be deallocated from other, and vice versa.
      // Always returns true for stateless allocators.
      bool operator== (const aligned_allocator& other) const;

      // Default constructor, copy constructor, rebinding constructor, and destructor.
      // Empty for stateless allocators.
      aligned_allocator (void);
      aligned_allocator (const aligned_allocator&);
      template <typename U>
      aligned_allocator (const aligned_allocator<U, Alignment>&);

      ~aligned_allocator();

      inline bool is_pow2 (void) const;

      // The following will be different for each allocator.
      T* allocate (const size_t n) const;

      void deallocate (T * const p, const size_t n) const;

      // The following will be the same for all allocators that ignore hints.
      template <typename U>
      T * allocate (const size_t n, const U * /* const hint */) const;

      // Allocators are not required to be assignable, so
      // all allocators should have a private unimplemented
      // assignment operator. Note that this will trigger the
      // off-by-default (enabled under /Wall) warning C4626
      // "assignment operator could not be generated because a
      // base class assignment operator is inaccessible" within
      // the STL headers, but that warning is useless.
   //private:
   //   aligned_allocator& operator= (const aligned_allocator&);
   };

   //template <typename T>
   //struct mmx_aligned_allocator : public aligned_allocator<T, __Align_m128> {};

   //template <typename T>
   //struct page_aligned_allocator : public aligned_allocator<T, __Align_page> {};

   //template <typename T>
   //struct cache_aligned_allocator : public aligned_allocator<T, __Align_cache> {};

   template <typename T>
   struct numa_interleaved_allocator
         : public aligned_allocator<T, __Align_page>
   {
      //public:

      typedef aligned_allocator<T,__Align_page>	Parent;

      // The following must be the same for all allocators.
      template <typename U>
      struct rebind
      {
         typedef numa_interleaved_allocator<U> other;
      };

      void deallocate (T * const p, const size_t n) const;
      T* allocate (const size_t n) const;
      template <typename U>
      T * allocate (const size_t n, const U * /* const hint */) const;
   };

#if 0
   template <typename T>
   struct numa_local_allocator
         : public aligned_allocator<T, __Align_page>
   {
//   public:
      typedef aligned_allocator<T,__Align_page>	Parent;

      // The following must be the same for all allocators.
      template <typename U>
      struct rebind
      {
         typedef numa_local_allocator<U> other;
      };

      void deallocate (T * const p, const size_t n) const;
      T* allocate (const size_t n) const;

      template <typename U>
      T * allocate (const size_t n, const U * /* const hint */) const;
   };
#endif

} // namespace parallel

#endif
