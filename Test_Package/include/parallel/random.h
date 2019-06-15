#ifndef __has_parallel_random_h
#define __has_parallel_random_h

#include <stdlib.h>

#include <parallel/parallel.h>

namespace parallel
{
   struct RandomGeneratorType
   {
      drand48_data state_buffer;
      bool initialized;

      RandomGeneratorType(void) : initialized(false) {}

      void init (long int seed = 0)
      {
         if (this->initialized == false)
         {
            this->initialized = true;

            seed += 100;

            srand48_r (seed, &this->state_buffer);

            for (int i = 0; i < 10; i++)
            {
               double rval;
               drand48_r(&this->state_buffer, &rval);
               printf("i=%d, seed=%d, rval=%lf\n", i, seed, rval);
            }
         }
      }

      double operator()(void)
      {
         double rval;
         drand48_r(&this->state_buffer, &rval);
         return rval;
      }
   };

   static RandomGeneratorType RandomGenerator[parallel::max_threads];
}

#endif
