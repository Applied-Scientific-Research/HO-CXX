#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "forall.hpp"

using namespace HighOrderFEM;

int main (int argc, char* argv[])
{
   int error = 0;

   {
      std::vector<int> v{0,1,2,3,4,5};
      std::vector<int> vout{1,2,3,4,5,6};

      forall( v.size(), [&](int i){ v[i] += 1; } );

      bool correct = true;
      for (int i = 0; i < v.size(); ++i)
         correct &= ( v[i] == vout[i] );

      if (not(correct))
      {
         fprintf(stderr,"Failed forall-1d test\n");
         fprintf(stderr,"Expected, Actual\n");
         for (int i = 0; i < v.size(); ++i)
            fprintf(stderr,"%d, %d\n", vout[i], v[i]);
      }

      error += 1;
   }

   {
      std::vector<double> v{0,1,2,3,4,5};
      //auto sum1 = forall_reduce( v.size(), [&](const int i) -> double { return v[i]*v[i]; } );
      auto func = [&](const int i) -> double { return v[i]*v[i]; };
      double sum1 = forall_reduce( v.size(), func );
      double expected1 = 55.0;
      printf("sum1= %e expected= %e error= %e\n", sum1, expected1, std::fabs(sum1-expected1)/expected1);

      Reduction::maxval< double> max;
      double max1 = forall_reduce( v.size(), func, max );
      double max2 = forall_reduce( v.size(), func, Reduction::maxval<double>() );
    //double max3 = forall_reduce( v.size(), func, maxval<int>() );
      double expected2 = 25.0;
      printf("max1= %e expected= %e error= %e\n", max2, expected2, std::fabs(max2-expected2)/expected2);
   }

   return 0;
}
