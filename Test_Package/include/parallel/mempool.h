#ifndef __parallel_mempool_h
#define __parallel_mempool_h

#ifdef USE_LIBNUMA
   #include <numa.h>
   #include <numaif.h>
   #if defined(LIBNUMA_API_VERSION) && (LIBNUMA_API_VERSION == 2)
      #warning 'including libnuma compat header'
      #include <numacompat1.h>
   #endif
//   #define USE_LIBNUMA (1)
#endif

#include <utils/config.h>

#include <unistd.h>
#include <stdint.h>
#include <errno.h>

#include <parallel/parallel.h>

namespace parallel
{
   class ThreadMemoryPool
   {
//      enum { _max_threads = 64 };
      enum { _max_threads = parallel::max_threads };

      typedef ThreadMemoryPool self;

      char	*memptr[_max_threads];
      int	onnode[_max_threads];
      size_t	memsize[_max_threads];
      bool	locked[_max_threads];

      int	max_threads;
      size_t	page_size;

      private:

      ThreadMemoryPool()
      {
         this->max_threads = sysconf( _SC_NPROCESSORS_ONLN );
         this->page_size = sysconf( _SC_PAGESIZE );

         assert(this->max_threads <= _max_threads);

         for (int i = 0; i < this->max_threads; ++i)
         {
            this->memptr[i] = NULL;
            this->memsize[i] = 0;
            this->onnode[i] = -1;
            this->locked[i] = false;
         }
         printf("initialized ThreadMemoryPool %d %lu\n", this->max_threads, this->page_size);
      }

      ThreadMemoryPool(const ThreadMemoryPool &);
      void operator=(const ThreadMemoryPool &);

      ~ThreadMemoryPool()
      {
         this->clear();
      }

      public:

      void clear (int i)
      {
         if (this->locked[i])
         {
            fprintf(stderr,"pool[%d] is locked\n", i);
            exit(-1);
         }
         if (this->memptr[i] && this->memsize[i])
         {
            free(this->memptr[i]);
            this->memptr[i] = NULL;
            this->memsize[i] = 0;
            printf("thread_memory_pool: cleared item %d\n", i);
         }
      }
      void clear (void)
      {
         for (int i = 0; i < max_threads; ++i)
            this->clear(i);
      }

      void* getPointer (int thread_id, size_t size, int node = -1)
      {
         // Is this a valid thread_id?
         if (thread_id < 0 || thread_id >= this->max_threads)
         {
            fprintf(stderr,"thread_id=%d is invalid %d\n", thread_id, this->max_threads);
            exit(-1);
         }

         // Is this locked?
         if (this->locked[thread_id])
         {
            fprintf(stderr,"thread memory pool %d is locked\n", thread_id);
            return NULL;
         }

         // Round the size request up to the maximum type size...
         if (size % sizeof(double))
         {
            size_t offset = sizeof(double) - (size % sizeof(double));
            //printf("size=%d, size %% sizeof(double) = %d, offset=%d, size+offset=%d\n", size, size % sizeof(double), offset, size + offset);
            size += offset;
         }

         // Is there enough space?
         bool realloc = (size > this->memsize[thread_id]) || (this->memptr[thread_id] == NULL);

         if (__DEBUG) printf("thread_id=%d, size=%lu, memsize=%lu, realloc=%d\n", thread_id, size, this->memsize[thread_id], realloc);

         if (realloc == true)
         {
            // Is this already allocated?
            this->clear(thread_id);

            int ierr = posix_memalign((void**)&this->memptr[thread_id], this->page_size, size);
            if (ierr)
            {
               if (ierr == EINVAL)
                  fprintf(stderr,"EINVAL error\n");
               else if (ierr == ENOMEM)
                  fprintf(stderr,"ENOMEM error\n");
               else
                  fprintf(stderr,"unknown error\n");
               exit(-1);
            }

            this->memsize[thread_id] = size;
         }

         this->locked[thread_id] = true;

         return this->memptr[thread_id];
      }

      void release (int thread_id)
      {
         this->locked[thread_id] = false;
         if (__DEBUG) printf("released %d\n", thread_id);
      }

      static self& getInstance(void)
      {
         static self instance;
         if (__DEBUG) printf("instance = %lx\n", &instance);
         return instance;
      }
   };

   //static thread_memory_pool_type thread_memory_pool;
   //thread_memory_pool_type thread_memory_pool;
   //Singleton<thread_memory_pool_type> thread_memory_pool_singleton;

#if 0
   template <typename _Tp>
   struct aligned_vector
   {
      typedef _Tp	value_type;

      value_type	*memptr; // pointer returned by malloc.
      value_type	*values; // pointer to the aligned address.
      size_t		num_elements; // length of the vector -- not the memory size!.
      size_t		max_elements; // maximum # of aligned elements -- not the memory size!.
      size_t		alignment;

      aligned_vector (size_t _alignment = sizeof(value_type))
         :
         memptr(NULL), values(NULL), num_elements(0), max_elements(0), alignment(_alignment)
      {}

      ~aligned_vector()
      {
         this->clear();
      }

      void clear (void)
      {
         if (this->memptr && this->max_elements)
         {
            free(this->memptr);
            if (__DEBUG > 0) printf("dealloc'd aligned_vector %lu\n", this->size());
            this->memptr = NULL;
            this->values = NULL;
            this->max_elements = 0;
            this->num_elements = 0;
         }
      }

      // Change the size of the internal memory.
      void reserve (const size_t new_max_elements)
      {
         // make sure we're a power-of-2.
         assert ((this->alignment & (this->alignment - 1)) == 0);

         if (new_max_elements <= this->max_elements)
            return;

         value_type *new_memptr = (value_type *) malloc(sizeof(value_type)*new_max_elements + (this->alignement-1));
         if (new_memptr == NULL)
         {
            fprintf(stderr,"Allocation error %s %d\n", __FILE__, __LINE__);
            exit(-1);
         }

         uintptr_t bitmask = ~((uintptr_t)(this->alignement - 1));

         value_type *new_values = (value_type *) (((uintptr_t)new_memptr + this->alignement - 1) & bitmask);

         if (__DEBUG > 0) printf("aligned vector %x %x %u\n", new_memptr, new_values, this->alignement);
         assert (((uintptr_t)new_values & ~bitmask) == 0);

         // Copy the old data ...
         size_t num_elements_to_copy = std::min(this->num_elements, new_max_elements);
         if (num_elements_to_copy > 0)
         {
            for (int i = 0; i < num_elements_to_copy; ++i)
               new_values[i] = this->values[i];
         }

         // Free the old data
         if (this->memptr && this->max_elements) free(this->memptr);

         // Reset the pointers.
         this->values = new_values;
         this->memptr = new_memptr;

         // Reset the limits.
         this->max_elements = new_max_elements;
      }

      // Change vector size and possibly realloc the memory.
      void resize (const size_t new_num_elements)
      {
         // If the old and new alignments match and the reserved space is sufficient, only update the vector length.
         if (new_num_elements <= this->max_elements)
         {
            //this->num_elements = new_num_elements;
         }
         else
         {
            // Otherwise, grow the storage.
            this->reserve(new_num_elements);
         }

         // Reset the vector length.
         this->num_elements = new_num_elements;

         return;
      }

      size_t size(void) const { return this->num_elements; }

      inline const value_type& operator[](const int i) const { return this->values[i]; }
      inline       value_type& operator[](const int i)       { return this->values[i]; }
   };
#endif

#if 0
   inline void partition (int &start, int &stop, int part_id, int num_parts, int lo, int hi)
   {
      // Total # of elements.
      int length = hi - lo;

      // Even chunck size.
      int chunk_size = (int)((float)length / num_parts);

      // Overshoot
      int remainder  = length - chunk_size*num_parts;
      //int remainder  = length % num_parts;

      start = lo + part_id * chunk_size;
      stop  = start + chunk_size;

      if (part_id < remainder) {
         start += part_id;
         stop  += (part_id + 1);
      }
      else {
         start += remainder;
         stop  += remainder;
      }
   }
   inline void partition (int parts[], int num_parts, int lo, int hi)
   {
      // Total # of elements.
      int length = hi - lo;

      // Even chunck size.
      int chunk_size = (int)((float)length / num_parts);

      // Overshoot
      int remainder  = length - chunk_size*num_parts;

      int sum = lo;
      for (int i = 0; i < num_parts; i++)
      {
         parts[i] = sum;
         sum += chunk_size;
         if (i < remainder) sum += 1;
         if (__DEBUG > 0) printf("parts[%d]=%d %d\n", i, parts[i], sum);
      }
      parts[num_parts] = hi; // is this necessary?

      if (__DEBUG > 0)
         for (int i = 0; i < num_parts; i++)
            printf("parts[%d]=%d %d\n", i, parts[i], parts[i+1]-parts[i]);
   }

   inline void* align_pointer (void *ptr0)
   {
      char *ptr = (char*)ptr0;

      // Align on a page boundary
      uintptr_t offset = ((uintptr_t)ptr) % numa_params.page_size;

      if (offset)
         if (offset < (numa_params.page_size / 2))
            ptr -= offset;
         else
            ptr += (numa_params.page_size - offset);

      if (__DEBUG > 0) printf("align_pointer: ptr=%x, old=%x, diff=%d\n", ptr, ptr0, (ptrdiff_t)(ptr - (char*)ptr0));

      return (void*)ptr;
   }

   inline int bind_to_node (void *ptr, size_t memsize, int to_node)
   {
      char *first = (char *)ptr;
      char *end = first + memsize;
      if (__DEBUG > 0) printf("first=%x, memsize=%lu, node=%d\n", first, memsize, to_node);

      if (1)
      {
         char *old = first;
#if 0
         // Align on a page boundary
         uintptr_t offset = ((uintptr_t)first) % numa_params.page_size;

         if (offset)
            if (offset < (numa_params.page_size / 2))
               first -= offset;
            else
               first += (numa_params.page_size - offset);
         if (__DEBUG > 0) printf("first=%x, offset=%lu, diff=%d\n", first, offset, (ptrdiff_t)(old - first));
#else
         first = (char*)align_pointer(first);
         //if (__DEBUG > 0) printf("first=%x, diff=%d\n", first, (ptrdiff_t)(old - first));
#endif
      }

      // Pin to the node.
      memsize = end - first;

      unsigned long mask = (1 << to_node);

      if (__DEBUG > 0) printf("memsize=%d, mask=%u\n", memsize, mask);

      if (memsize < 0)
      {
         fprintf(stderr,"page alignement made negative length -- too small\n");
         exit(-1);
      }

#if defined(USE_LIBNUMA)
#warning 'binding memory to numa nodes'
      errno = 0;
      if ((mbind (first, memsize, MPOL_BIND, &mask, 64, MPOL_MF_MOVE)) != 0)
      {
         perror("mbind...");
         exit(-1);
      }
#endif

      return 0;
   }

   inline int get_node_location (void *ptr)
   {
      ptr = align_pointer(ptr);

      int status[1]; // node location array
      status[0] = -1;

      void *pages[1]; // list of pages
      pages[1] = ptr;

#if defined(LIBNUMA_API_VERSION) && (LIBNUMA_API_VERSION > 1)
      errno = 0;
      int ierr = move_pages (0,		// this process
                             1,		// node count to move
                             pages,	// list of pages (pointers) to move/check
                             NULL,	// node list -- or NULL if query only
                             status,	// status of each page -- and node location for query
                             0);	// flag
      if (ierr || errno)
      {
         perror("error in move_pages: ");
         exit(-1);
      }

      if (status[0] < 0)
      {
         int err = -status[0];
              if (err == EACCES) fprintf(stderr,"move_pages: is mapped by multiple processes\n");
         else if (err == EBUSY ) fprintf(stderr,"move_pages: is busy\n");
         else if (err == EFAULT) fprintf(stderr,"move_pages: is zero\n");
         else if (err == EIO   ) fprintf(stderr,"move_pages: is EIO\n");
         else if (err == EINVAL) fprintf(stderr,"move_pages: is EINVAL\n");
         else if (err == ENOENT) fprintf(stderr,"move_pages: is not present\n");
         else if (err == ENOMEM) fprintf(stderr,"move_pages: is too large -- out of memory\n");
         else                    fprintf(stderr,"move_pages: unknown error\n");
      }
      else
         printf("move_pages: ptr=%x, node=%d\n", ptr, status[0]);
#else
      printf("move_pages not available in libnuma version < 2\n");
#endif

      return status[0];
   }

   template <typename Vector>
   void numa_interleave_vector (Vector &v)
   {
      typedef typename Vector::value_type value_type;
#if defined(USE_LIBNUMA)
#warning 'interleaving memory across numa nodes'
      value_type *values = v.values;
      size_t num_values  = v.size();

      errno = 0;
      numa_interleave_memory (values, num_values*sizeof(value_type), &numa_params.all_nodes_mask);
      if (errno)
      {
         perror("error in numa_interleave_memory: ");
         exit(-1);
      }
#endif
   }

   template <typename Vector, typename Vector2>    
   void numa_distribute_vector (Vector &v, const Vector2 &partition)
   {
      typedef typename Vector::value_type value_type;

      value_type *values = v.values;
      size_t num_values  = v.size();
      int num_partitions = partition.size()-1;

      for (int n = 0; n < num_partitions; n++)
      {
         size_t num_values_i = partition[n+1] - partition[n];
         value_type *values_start_i = &values[ partition[n] ];
//         value_type *values_stop_i  = values_start_i + num_values_i;

         if (__DEBUG > 0) printf("n=%d, num_values_i=%d, values_start_i=%x\n", n, num_values_i, values_start_i);
#if 0
         numa_tonode_memory (values_start_i, num_values_i * sizeof(value_type), n);
         if (errno)
         {
            perror("error in numa_interleave_memory: ");
            exit(-1);
         }
#else
         bind_to_node (values_start_i, num_values_i*sizeof(value_type), n);
#endif
      }
   }

#endif

} // splib

#endif
