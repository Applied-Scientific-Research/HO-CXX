#ifndef __has_parallel_numa_h
#define __has_parallel_numa_h

#if defined(USE_LIBNUMA)
   #ifndef _OPENMP
      #warning 'OpenMP is required for libnuma ... disabling libnuma'
      #undef USE_LIBNUMA
   #endif
#endif

#if defined(USE_LIBNUMA)
#warning 'Using libnuma ...'
   #include <numa.h>
   #include <numaif.h>
   #if defined(LIBNUMA_API_VERSION) && (LIBNUMA_API_VERSION == 2)
      #warning 'including libnuma compat header'
      #include <numacompat1.h>
   #endif
#endif

#include <parallel/parallel.h>
#include <splib/attributes.h>
#include <splib/base_matrix.h>

#include <vector>

namespace parallel
{
   struct ParallelPropertiesType
   {
      unsigned int attributes;

      // Optional row partitioning information
      int num_partitions;
      std::vector<int> partition;

      ParallelPropertiesType(void);
   };

   namespace numa_options
   {
      enum type
      {
         none        = 0,
         interleaved = (1 << 1),
         distributed = (1 << 2)
      };
   }

   struct numa_params_type
   {
      int lib_exists;
      int is_initialized;
      int num_procs;
      int num_nodes;
      int max_threads;
      int num_threads_per_node;
      size_t page_size;

      // Smallest vector we'll consider for numafication ...
      size_t min_pages;

#ifdef USE_LIBNUMA
      nodemask_t all_nodes_mask;
#endif

      std::vector<size_t> node_size;

      numa_params_type(void);

      int initialize(void);
   };

   //static numa_params_type numa_params;
   extern numa_params_type numa_params;

   int numa_lock_threads_to_nodes(void);
   int numa_unlock_threads_to_nodes(void);

   template <typename Vector, typename Vector2>
   int numa_distribute_vector (Vector &v, const Vector2 &partition);

   template <typename Vector>
   int numa_distribute_vector (Vector &v);

   template <typename Vector>    
   int numa_interleave_vector (Vector &v);

   template <typename Vector>
   int numafy_vector (Vector &v, const int _distribute = 0);

   template <typename Matrix>
   int numafy_matrix (Matrix &A);

   int numafy_base_matrix (splib::BaseMatrix *baseA);

} // namespace parallel

#endif
