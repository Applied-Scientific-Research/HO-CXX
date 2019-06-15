#ifndef __has_precon_point_relax_h
#define __has_precon_point_relax_h

namespace precon
{
namespace point_relax
{
   namespace PointRelaxOptions
   {
      enum { None		= 0,
             ZeroGuess		= 1,
             Jacobi		= 2,
             GS_Forward		= 4,
             GS_Backward	= 8 };

      // Method used to estimate the largest eigenvalue of the Jacobi iteration matrix.
      enum { Arnoldi		= 1,
             Lanczos		= 2,
             Maxnorm		= 4 };

      bool is_zero_guess (int options);
      bool is_jacobi (int options);
      bool is_gs_forward (int options);
      bool is_gs_backward (int options);

      int enable_zero_guess (int options);
      int enable_jacobi (int options);
      int enable_gs_forward (int options);
      int enable_gs_backward (int options);

      int disable_zero_guess (int options);
      int disable_jacobi (int options);
      int disable_gs_forward (int options);
      int disable_gs_backward (int options);

   } // namespace PointRelaxOptions

   template <typename Matrix, typename Vector1, typename Vector2>
   void jacobi (const Matrix	&A,
                      Vector1	&x,
                const Vector2	&b,
                const typename Matrix::value_type omega = typename Matrix::value_type(1),
                const int	options = PointRelaxOptions::Jacobi);

   template <typename Matrix, typename Vector1, typename Vector2>
   void sgs (const Matrix	&A,
                   Vector1	&x,
             const Vector2	&b,
             const typename Matrix::value_type omega = typename Matrix::value_type(1),
             const int		options = PointRelaxOptions::None);

   template <typename Matrix, typename Vector1, typename Vector2>
   void gs (const Matrix	&A,
                  Vector1	&x,
            const Vector2	&b,
            const typename Matrix::value_type omega = typename Matrix::value_type(1),
            const int		options = PointRelaxOptions::None);

   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
   void gs_cf (const Matrix	&A,
                     Vector1	&x,
               const Vector2	&b,
               const Vector3	&c_points,
               const int	parallel_cc,
               const int	parallel_ff,
               const typename Matrix::value_type omega = typename Matrix::value_type(1),
               const int	options = PointRelaxOptions::None);

} // namespace point_relax
} // namespace precon

#endif
