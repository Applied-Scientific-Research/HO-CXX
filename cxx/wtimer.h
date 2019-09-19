#ifndef __wtimer_h
#define __wtimer_h

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

namespace HighOrderFEM
{

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

} // namespace

#endif
