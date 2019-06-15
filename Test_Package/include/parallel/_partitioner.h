#ifndef __has_parallel_partitioner_h
#define __has_parallel_partitioner_h

namespace parallel
{

   inline void partitioner (int &start, int &stop, int part_id, int num_parts, int lo, int hi)
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
   template <typename T>
   inline void partitioner (T parts[], int num_parts, const T lo, const T hi)
   {
      // Total # of elements.
      T length = hi - lo;

      // Even chunck size.
      T chunk_size = (T)((float)length / num_parts);

      // Overshoot
      T remainder = length - chunk_size*num_parts;

      T sum = lo;
      for (int i = 0; i < num_parts; i++)
      {
         parts[i] = sum;
         sum += chunk_size;
         if (i < remainder) sum += 1;
         if (__DEBUG > 0) printf("parts[%d]=%d %d\n", i, parts[i], sum);
      }
      parts[num_parts] = hi; // is this necessary?

//      if (__DEBUG > 0)
//         for (int i = 0; i < num_parts; i++)
//            printf("parts[%d]=%d %d\n", i, parts[i], parts[i+1]-parts[i]);
   }

} // parallel

#endif
