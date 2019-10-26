#ifndef __wtimer_h
#define __wtimer_h

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#ifdef _OPENMP
# include <omp.h>
#endif

namespace HighOrderFEM
{

#ifdef _OPENMP
#warning 'Using OpenMP wall-clock timer'
typedef double TimerType;

TimerType getTimeStamp(void)
{
   return omp_get_wtime();
}
double getElapsedTime( const TimerType& a, const TimerType& b )
{
   return (b - a);
}
#else
typedef std::chrono::high_resolution_clock::time_point TimerType;

TimerType getTimeStamp(void)
{
   return std::chrono::high_resolution_clock::now();
}

double getElapsedTime( const TimerType& a, const TimerType& b )
{
   using namespace std::chrono;

   std::chrono::duration<double> t_span = duration_cast< duration<double> >(b - a);

   return t_span.count();
}
#endif

} // namespace

#endif
