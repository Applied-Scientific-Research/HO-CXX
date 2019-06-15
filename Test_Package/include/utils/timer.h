#ifndef __has_utils_timer_h
#define __has_utils_timer_h

#include <utils/config.h>
#include <stdio.h>

namespace utils
{

// Generic timer ...
struct Timer
{
   typedef Timer self_type;

   double m_time;
   double m_time_start, m_time_stop;
   bool m_running;

   Timer (void);
   Timer (const double& time);
   Timer (const Timer& t);

   double operator-(const Timer& t) const;
   double operator+(const Timer& t) const;

   //self_type & operator= (const Timer& t);
   self_type & operator+=(const Timer& t);
   self_type & operator-=(const Timer& t);

   //const double & operator&(void) const;
   //      double & operator&(void)      ;

   void start (void);
   void stop  (void);
   void reset (void);
   double to_date (void) const;

   const double & time (void) const;
         double & time (void)      ;
};

// Specific singleton timers for cummulative operations ...
class GlobalTimers
{
public:
   typedef GlobalTimers T;

   // These are the defined timers ...

   Timer spmv, spmm, gemv;

   static T& getInstance()
   {
      static T instance; // Guaranteed to be destroyed.
                         // Instantiated on first use.
      return instance;
   }

   ~GlobalTimers(void)
   {
      if (__PROFILE)
      {
         printf("GlobalTimers::SpMV = %f (ms)\n", 1000.*this->spmv.time());
         printf("GlobalTimers::SpMM = %f (ms)\n", 1000.*this->spmm.time());
         printf("GlobalTimers::GEMV = %f (ms)\n", 1000.*this->gemv.time());
      }
   }

private:
   GlobalTimers() {};        // Constructor? (the {} brackets) are needed here.

                             // Dont forget to declare these too. You want to make sure they
                             // are unaccessable otherwise you may accidently get copies of
                             // your singleton appearing.
   GlobalTimers(T const&);   // Don't Implement (copy constructor)
   void operator=(T const&); // Don't implement (assignment operator)
};

} // utils

#endif
