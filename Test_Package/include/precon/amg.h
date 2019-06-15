#ifndef __has_precon_amg_h
#define __has_precon_amg_h

#include <precon/precons.h>
#include <precon/base_precon.h>
#include <precon/color.h>
#include <precon/point_relax.h>

#include <splib/dense_matrix.h>
#include <splib/vector.h>
#include <splib/csr_matrix.h>
#include <splib/base_matrix.h>
#include <splib/hyb_matrix.h>

#include <utils/timer.h>
#include <utils/type_traits.h>

#include <parallel/parallel.h>
#include <parallel/numa.h>
//#include <parallel/partitioned_objects.h>

#include <solver/base_solver.h>

#include <map>
#include <string>

namespace precon
{
   namespace SmootherTags
   {
      enum Tag { WJ = 1, SGS, GS };
   }

   std::string smoother_name (SmootherTags::Tag tag);

   struct SmootherOptions
   {
      SmootherTags::Tag	tag;
      int		num_pre_iters;
      int		num_post_iters;
      double		omega;
      double		rho;
      int		spectral_radius_method;
      int		coloring_threshold;
      int		max_colors;
      bool		color_matrix;
      double		coloring_truncation_rate;

      SmootherOptions(void);
      void print(FILE *fp);
   };

   namespace CoarseGridSolverTags
   {
      enum Tag { LU = 1, SVD, SLU, ITER };
   }

   std::string coarse_grid_solver_name (CoarseGridSolverTags::Tag tag);

   struct CoarseGridSolver
   {
      enum { max_size = 500 }; // max size for a dense solver

      typedef double	value_type;
      typedef int	index_type;

      typedef splib::DenseMatrix<value_type>	matrix_type;
      typedef splib::Vector<index_type>		pivot_type;
      typedef splib::Vector<value_type>		vector_type;

      CoarseGridSolverTags::Tag		solver_tag;
      matrix_type	LU;
      matrix_type	U, VT;
      pivot_type	p;
      vector_type	s, work;
      void		*slu_data; // internal storage handle for SuperLU.
      bool		slu_data_is_double;
      //float		*slu_data_f;
      //double		*slu_data_d;
      //solver::BaseSolver *iterSolver; // internal storage for a sparse iterative solver.
      void              *itsol_data; // internal storage for a sparse iterative solver.

      utils::Timer	timer;

      std::string name (void) const;

      CoarseGridSolver (void);
      ~CoarseGridSolver();

      template <typename Matrix>
      void setup (const Matrix &A);

      template <typename Vector1, typename Vector2>
      void solve (const Vector1 &b, Vector2 &x) const;

      template <typename Vector1, typename Vector2>
      void operator() (const Vector1 &b, Vector2 &x) const;
   };

   template <typename MatrixType>
   struct Level
   {
      typedef MatrixType			matrix_type;
      typedef typename MatrixType::value_type	value_type;
//      typedef typename MatrixType::index_type	index_type;
      typedef typename MatrixType::format_type	format_type;

      // Level # if in the hierarchy
      int			level;

      // (Pointer-to) the input matrix
      MatrixType		*A;
      splib::csr_matrix<float>	A_float;
      bool			keep_A;

      // Restriction matrix
      MatrixType		R;

      // Prolongation matrix
      MatrixType		P;

      // Strength-of-connection matrix
      MatrixType		S;

#define __enable_amg_baseMatrix
#if defined(__enable_amg_baseMatrix)
      splib::BaseMatrix		*baseA, *baseR, *baseP;
      bool			baseA_allocated, baseR_allocated, baseP_allocated;
#endif

      // C/F splitting vector
      splib::Vector<int>	cf_marker;
      splib::Vector<int>	cf_points;
      //size_t                    ncc, nff;

      // Multi-color matrix for MC-GS
      typedef point_relax::ColoredMatrix<value_type>	colored_matrix_type;
      colored_matrix_type	colored_matrix;

      // Scratch vectors
      typedef splib::Vector<value_type> vector_type;
      vector_type r, b_coarse, x_coarse;

      // Pre- and Post-smoothing data
      SmootherOptions		smoother_options;
      splib::Vector<value_type>	smoother_data_dinv;

#ifdef __EnablePartitioning
      typedef splib::Vector<value_type,parallel::numa_interleaved_allocator<value_type> > interleaved_vector_type;
      typedef parallel::partitioned_vector<vector_type> partitioned_vector_type;
      typedef parallel::partitioned_matrix<MatrixType> partitioned_matrix_type;

      interleaved_vector_type interleaved_x_coarse, interleaved_r;
      partitioned_vector_type partitioned_b_coarse;
      partitioned_matrix_type partitioned_A, partitioned_R, partitioned_P;
#endif

      // timing objects
      utils::Timer		smoother_timer;
      utils::Timer		timer;

      Level (int _level = -1);

      ~Level (void);
   };

   namespace StrengthTags
   {
      enum Tag { ClassicalStrength = 1,
                 SymmetricStrength,
                 NullStrength };
   }

   namespace CoarsenTags
   {
      enum Tag { ClassicCoarsen = 1,
                 PMISCoarsen };
   }

   namespace InterpTags
   {
      enum Tag { ClassicInterp = 1,
                 DirectInterp,
                 ExtendedInterp,
                 ExtendedPlusiInterp,
                 ExtendedPlusiiffInterp };
   }

template <typename Matrix>
struct AMG
   : public BasePrecon
{
   enum { PreconTag = PreconTags::AMG };
   enum { ValueTag = Matrix::ValueTag };
   enum { FormatTag = Matrix::FormatTag };

   typedef BasePrecon	Parent;

   typedef typename Matrix::value_type	value_type;
   typedef typename Matrix::index_type	index_type;

   typedef splib::csr_matrix<value_type> matrix_type;

   typedef Level<matrix_type> level_type;

   int		max_levels;
   int          max_coarse_grid_size;
   int		num_levels;
   value_type	theta;
   bool		mixed_precision_smoother;
   value_type	interp_trunc_factor;
   utils::Timer	timer;

   StrengthTags::Tag	strength_method;
   CoarsenTags::Tag	coarsen_method;
   InterpTags::Tag	interp_method;

   // List of levels ...
   level_type	*levels;

   // Pointer to a CSR-formatted version on the finest level -- may just be the input matrix ...
   matrix_type	*input_A;
   bool		A_is_allocated;
   bool		print_info;

   matrix_type  _S, _ST;

   CoarseGridSolver/*<matrix_type>*/	coarse_grid_solver;

   // Smoother options ... passed to the smoother build routines.
   typedef std::map<int,SmootherOptions> smoother_options_type;
   smoother_options_type smoother_options;

   //AMG (void);
   AMG (const BasePreconOptions& options);
   ~AMG();

   void numafy (void);
   void partition_data (void);

   template <typename OtherMatrixType>
   int build (const OtherMatrixType& A);

   int build (const splib::BaseMatrix *baseA);

   void print_hierarchy (void) const;

   int tag (void) const;
   int value_tag (void) const;
   int format_tag (void) const;
   std::string name (const bool full_name = true) const;

   // Solver Ax=b .... i.e., return x = A^(-1)b
   template <typename Vector1, typename Vector2>
   int solve (const Vector1& b, Vector2& x) /*const*/;

   //template <typename Vector1, typename Vector2>
   //int solve (const parallel::partitioned_vector<Vector1>& b, Vector2& x) /*const*/;

   //int solve (const splib::Vector<double>& b, splib::Vector<double>& x) /*const*/;
   //int solve (const splib::Vector<value_type>& b, splib::Vector<value_type>& x) /*const*/;
   //int solve (const splib::Vector<double>& b, splib::Vector<float >& x) /*const*/;
   //int solve (splib::BaseVector *baseb, splib::BaseVector *basex) /*const*/;

};

} // namespace precon

#endif
