#ifndef _solver_base_h
#define _solver_base_h

#include <solver/solvers.h>
#include <precon/precons.h>
#include <precon/base_precon.h>
#include <precon/identity.h>

#include <splib/vector.h>
#include <splib/base_matrix.h>
#include <splib/csr_matrix.h>
#include <splib/coo_matrix.h>
#include <splib/dia_matrix.h>

#include <string>

namespace solver
{

struct BaseSolver
{
   //int m_solver_tag;
   SolverTags::Tag solver_tag;
   //int m_precon_tag;

   //typedef splib::BaseMatrix	BaseMatrix;
   //typedef splib::BaseVector	BaseVector;
   //typedef precon::BasePrecon	BasePrecon;

   // Limits
   int max_iters;
   double max_resid;

   // Actual counters and final residual
   int num_iters;
   double resid;

   bool abs_convergence;
   int verbosity;

   int monitor_width;
   double min_rate;
   splib::Vector<double> resid_monitor;

   // Keep a copy of the Identity precon handy.
   precon::Identity _I;

   BaseSolver (SolverTags::Tag tag);
   //BaseSolver (SolverTags::Tag solver_tag, precon::PreconTags::Tag precon_tag);
   //BaseSolver (const int solver_tag, const int precon_tag);
   //BaseSolver (const int solver_tag, const int precon_tag);
   virtual ~BaseSolver ();

   virtual std::string name (void) const = 0;
   virtual int tag (void) const = 0;
   virtual int value_tag (void) const = 0;
   //virtual int format_tag (void) const = 0;
   //virtual int precon_tag (void) const = 0;

   //template <typename Matrix>
   //int build (const Matrix &A);
   //virtual int build (const splib::BaseMatrix &A);
   virtual int build (const splib::BaseMatrix *baseA);
   //virtual int setPrecon (const precon::BasePrecon *baseM);
   //virtual precon::BasePrecon * getPrecon (void) const = 0;

   //template <typename Matrix, typename Vector1, typename Vector2, typename Precon>
   //int solve (const Matrix& A, Vector1 &u, const Vector2& f, const Precon& M);
   //virtual int solve (splib::BaseVector *x, const splib::BaseVector *b) = 0;
   //virtual int solve (splib::BaseVector *x, const splib::BaseVector *b, const precon::BasePrecon *basePrecon) = 0;
   //virtual int solve (splib::Vector<double> &x, const splib::Vector<double> &b, /*const*/ precon::BasePrecon *basePrecon) = 0;

   //template <typename Matrix, typename Vector1, typename Vector2>
   //int solve (const Matrix& A, Vector1 &u, const Vector2& f);
   virtual int solve (splib::BaseMatrix *A,
                      splib::BaseVector  *x,
                      splib::BaseVector  *b,
                      precon::BasePrecon *M);// = 0;

   virtual int monitor (BaseSolver *solver, double resid, double resid_1);
};

//int solve (BaseSolver *solver, splib::BaseMatrix *A, splib::BaseVector *x, const splib::BaseVector *b, const precon::BasePrecon *M);

} // namespace solver

#endif
