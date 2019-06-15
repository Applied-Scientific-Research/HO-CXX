#ifndef __has_parallel_scan_h
#define __has_parallel_scan_h

#include <utils/timer.h>
#include <splib/attributes.h>
#include <splib/spblas.h>

#include <parallel/parallel.h>
#include <parallel/partitioner.h>

namespace parallel
{

namespace details
{

   // In-place parallel prefix scan (sum, etc.) ...
   // takes v of length n+1 and returns
   // v[1:n] = Sum_{k=1}^i v[k], with v[0] = 0
   template <typename Vector>
   void inplace_scan (Vector& v)
   {
      typedef typename Vector::value_type value_type;

      const value_type Zero(0);

      const int n = v.size();

#ifdef _OPENMP
      const int num_threads = parallel::get_num_threads();

      static value_type sums[parallel::max_threads];

      #pragma omp parallel
      {
         const int thread_id = omp_get_thread_num();

         // 1) Compute the segment summation.
         value_type my_sum(0);

         int i_start, i_stop;
         parallel::partitioner(i_start, i_stop, thread_id, num_threads, 0, n-1);

         for (int i = i_start; i < i_stop; ++i)
            my_sum += v[i];

         // 2) Write this segement sum to the shared sums[] ... and wait
         sums[thread_id] = my_sum;

         #pragma omp barrier

         value_type left_sum(0);

         // 3) Accumlate the sums from the left ...
         for (int k = 0; k < thread_id; ++k)
            left_sum += sums[k];

         // 4) Pass-2 updates the array ...
         for (int i = i_start; i < i_stop; ++i)
         {
            value_type vi = v[i];
            v[i] = left_sum;
            left_sum += vi;
         }

         if (thread_id == num_threads-1)
            v[n-1] = left_sum;
      }
#else
      value_type sum(0);
      for (int i = 0; i < n-1; ++i)
      {
         value_type vi = v[i];
         v[i] = sum;
         sum += vi;
      }
      v[n-1] = sum;
#endif
   }

} // details

template <typename Vector>
void
inplace_scan (Vector& v)
{
   utils::Timer timer;
   if (__PROFILE) timer.start();
   details::inplace_scan(v);
   if (__PROFILE) { timer.stop(); printf("inplace_scan: v.size() = %d, time = %f (ms)\n", v.size(), timer.time()); }
}

} // namespace parallel

#endif
