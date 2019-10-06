#ifndef __basis_hpp
#define __basis_hpp

#include "array.hpp"

namespace HighOrderFEM
{

namespace details
{
}

template <int _K, int _L, typename T = double>
struct BasisFunctions
{
   typedef T value_type;

   enum { K = _K };
   enum { L = _L };

   constexpr int getSolutionOrder(void) const { return (int)K; }
   constexpr int getGeometryOrder(void) const { return (int)L; }

   StaticArrayType< value_type[K] > sps;
   StaticArrayType< value_type[L] > gps;
   StaticArrayType< value_type[K] > wgt;

   StaticArrayType< value_type[2][K] > NodesRadau; // Left/Right basis.
   StaticArrayType< value_type[2][K] > NodesGradRadau;

   StaticArrayType< value_type[K][K] > SolnNodesLgrangeBasis, SolnNodesGradLgrangeBasis;
   StaticArrayType< value_type[2][K] > SolnBndryLgrangeBasis, SolnBndryGradLgrangeBasis;
   StaticArrayType< value_type[2][L] > GeomBndryLgrangeBasis, GeomBndryGradLgrangeBasis;
   StaticArrayType< value_type[K][L] > GeomNodesLgrangeBasis, GeomNodesGradLgrangeBasis;

   BasisFunctions(void)
   {
      static_assert( K >= 1 && K <= 5, "Testing 1 <= K <= 5 failed.");

      if (K == 1) {
         sps(0) = 0.0;                                     // this case doesn't work; I'd need to do special cases
         wgt(0) = 2.0;                                     // Gaussian weights
      }
      else if (K == 2) {
         sps(1) = 1.0/sqrt(3.0);                           // collocation nodes per elem are {-C1, C1}
         sps(0) = -sps(1);
         wgt(1) = 1.0;                                     // Gaussian weights
         wgt(0) = wgt(1);
      }
      else if (K == 3) {
         sps(2) = sqrt(0.6);                                // collocation nodes per elem are {-C1, 0, C1}
         sps(1) = 0.;
         sps(0) = -sps(2);                                  // sloppy but quick; we really don't need to save the -ve values
         wgt(2) = 5./9.;                                    // Gaussian weights
         wgt(1) = 8./9.;
         wgt(0) = wgt(2);
      }
      else if (K == 4) {
         sps(3) = sqrt((15.+2.*sqrt(30.))/35.);      // collocation nodes per elem are {-C2, -C1, C1, C2}
         sps(2) = sqrt((15.-2.*sqrt(30.))/35.);
         sps(1) = -sps(2);
         sps(0) = -sps(3);
         wgt(3) = (18.-sqrt(30.))/36.;                // Gaussian weights
         wgt(2) = (18.+sqrt(30.))/36.;
         wgt(1) = wgt(2);
         wgt(0) = wgt(3);
      }
      else if (K == 5) {
         sps(4) = sqrt(5.+2.*sqrt(70.)/7.)/3.;        // collocation nodes per elem are {-C2, -C1, 0, C1, C2}
         sps(3) = sqrt(5.-2.*sqrt(70.)/7.)/3.;
         sps(2) = 0.;
         sps(1) = -sps(3);
         sps(0) = -sps(4);
         wgt(4) = (322.-13.*sqrt(70.))/900.;          // Gaussian weights
         wgt(3) = (322.+13.*sqrt(70.))/900.;
         wgt(2) = 128./225.;
         wgt(1) = wgt(3);
         wgt(0) = wgt(4);
      }

      static_assert( L >= 2 && L <= 5, "Testing 1 <= L <= 5 failed.");

      if (L == 2) {
         gps(1) = 1.;
         gps(0) = -gps(1);
      }
      else if (L == 3) {
         gps(2) = 1.;
         gps(1) = 0.;
         gps(0) = -gps(2);
      }
      else if (L == 4) {
         gps(3) = 1.;
         gps(2) = 1./3.;
         gps(1) = -gps(2);
         gps(0) = -gps(3);
      }

      for (int i = 0; i < K; ++i)
         printf("sps[%d]= %e %e\n", i, sps(i), wgt(i));

      for (int i = 0; i < L; ++i)
         printf("gps[%d]= %e\n", i, gps(i));

      this->getRadauBasis();
   }

private:

   void getRadauBasis (void)
   {
      if (K == 1) {
         NodesRadau(0,0) = NodesRadau(0,1) = 0.5;
         NodesGradRadau(0,0) = -0.5;
         NodesGradRadau(0,1) =  0.5;
      }
      else {
         double coef  = ( K % 2 == 0 ) ? 0.5 : -0.5;
         double coefd = coef * K;

         for (int i = 0; i < K; ++i)
         {
           // R_k(x)|Right = R_k(-x)|Left
           NodesRadau(i,0) = coef * (std::legendre(K, sps(i)) - std::legendre(K-1, sps(i)));
           NodesRadau(i,1) = coef * (std::legendre(K,-sps(i)) - std::legendre(K-1,-sps(i)));

           // D[R_k(x),x]|Right = -D[R_k(-x),x]|Left
           NodesGradRadau(i,0) = coefd * (std::legendre(K, sps(i)) + std::legendre(K-1, sps(i))) / (1.0 + sps(i));
           NodesGradRadau(i,1) =-coefd * (std::legendre(K,-sps(i)) + std::legendre(K-1,-sps(i))) / (1.0 - sps(i));

         }
      }

      for (int i = 0; i < K; ++i)
         printf("NodesRadau(%d,:)= %e %e\n", i, NodesRadau(i,0), NodesRadau(i,1));
      for (int i = 0; i < K; ++i)
         printf("NodesGradRadau(%d,:)= %e %e\n", i, NodesGradRadau(i,0), NodesGradRadau(i,1));
   }


};

} // end namespace

#endif
