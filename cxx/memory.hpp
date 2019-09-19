#ifndef __memory_h
#define __memory_h

#include <cstdio>
#include <cstdlib>

namespace HighOrderFEM
{

template <typename T>
T* allocate ( const size_t nelems )
{
   T *ptr = (T*) malloc( nelems * sizeof(T) );
   if (ptr == NULL)
   {
      fprintf(stderr,"Allocation error for %lu elements\n", nelems);
      return NULL;
   }
   else
      return ptr;
}

template <typename T>
int allocate ( const size_t nelems, T* &ptr )
{
   ptr = allocate<T>( nelems );
   if (ptr == NULL)
   {
      fprintf(stderr,"Allocation error for %lu elements at %d %s\n", nelems, __LINE__, __FILE__);
      exit(-1);
      return 1;
   }
   else
      return 0;
}

template <typename T>
void deallocate ( T* &ptr )
{
   if (ptr)
   {
      free(ptr);
      ptr = NULL;
   }
}

template <typename T>
int reallocate ( const size_t nelems, T* &ptr )
{
   T *new_ptr = (T*) realloc( (void *) ptr, nelems * sizeof(T) );
   if ( new_ptr == NULL )
   {
      fprintf(stderr,"Reallocation error for %lu elements at %d %s\n", nelems, __LINE__, __FILE__);
      deallocate( ptr );
      exit(-1);
      return 1;
   }
   else
   {
      ptr = new_ptr;
      return 0;
   }
}

} // namespace

#endif
