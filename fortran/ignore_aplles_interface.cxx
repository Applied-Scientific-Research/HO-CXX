#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include <typeinfo>

#include <utils/config.h>
#include <utils/type_traits.h>
#include <utils/timer.h>

//#include <patch.h>

#ifdef _OPENMP
#include <omp.h>
#include <parallel/parallel.h>
#include <parallel/scoped_lock.h>
#include <parallel/numa.h>
#endif

#include <cstring>
using namespace std;

#include <splib/formats.h>

#include <splib/vector.h>
#include <splib/csr_matrix.h>
#include <splib/coo_matrix.h>
#include <splib/dia_matrix.h>
#include <splib/ell_matrix.h>
#include <splib/hyb_matrix.h>
#include <splib/copy.h>
#include <splib/filter.h>

#include <precon/precons.h>
#include <precon/identity.h>
#include <precon/relaxation.h>
#include <precon/amg.h>
#include <precon/krylov.h>

#include <solver/solvers.h>
#include <solver/pcg.h>
#include <solver/bicgstab.h>
#include <solver/gmres.h>

static int verbose = 1;
static int write_matrix = 0;

//static int	format_tag	= splib::FormatTags::CSR;

//typedef double Float;

#define NDIM 3

//typedef PatchVariableType<Float,NDIM> NodeVariable;
//typedef PatchIndexType<NDIM> Index;


template <typename T, typename ConstructorData>
T* _allocate (const ConstructorData &data)
{
   T *obj = new T(data);
   assert(obj != NULL);
   return obj;
}
template <typename T>
T* _allocate (void)
{
   T *obj = new T;
   assert(obj != NULL);
   return obj;
}

template <typename T>
T* _destroy (T *ptr)
{
   if (ptr == NULL)
      return NULL;

//   if (verbose) printf("APLLES_Destroy: ptr = %p %s %s\n", ptr, typeid(ptr).name(), typeid(*ptr).name());

   delete ptr;

   return NULL;
}

template <typename Matrix, typename ValueType>
precon::BasePrecon* create_precon_1 (char *PreconName, precon::BasePreconOptions &options, ValueType)
{
   using namespace precon;

   std::string precon_name_comp = PreconName;
   // Turn to lower case to allow for comparison w/ case-insensitive input
   for (unsigned int i = 0; i < precon_name_comp.size(); i++)
       precon_name_comp[i] = tolower(PreconName[i]);

   typedef typename Matrix::template rebind<ValueType>::other PreconMatrixType;
   if (verbose) printf("PreconMatrix = %s\n", typeid(PreconMatrixType).name());

   precon::BasePrecon *base_precon = NULL;

   if (precon_name_comp.compare("none") == 0 || precon_name_comp.compare("identity") == 0)
   {
      typedef precon::Identity precon_type;
      base_precon = _allocate<precon_type>(options);
   }
   else if (precon_name_comp.compare("jacobi") == 0)
   {
      typedef precon::Relaxation<PreconMatrixType,PreconTags::Jacobi> precon_type;
      base_precon = _allocate<precon_type>(options);
   }
   else if (precon_name_comp.compare("diagonal") == 0)
   {
      typedef precon::Relaxation<PreconMatrixType,PreconTags::Diagonal> precon_type;
      base_precon = _allocate<precon_type>(options);
   }
   else if (precon_name_comp.compare("sgs") == 0)
   {
      typedef precon::Relaxation<PreconMatrixType,PreconTags::SGS> precon_type;
      base_precon = _allocate<precon_type>(options);
   }
   else if (precon_name_comp.compare("krylov") == 0)
   {
      typedef precon::Krylov<PreconMatrixType> precon_type;
      base_precon = _allocate<precon_type>(options);
   }
   else if (precon_name_comp.compare("amg") == 0)
   {
      typedef precon::AMG<PreconMatrixType> precon_type;
      base_precon = _allocate<precon_type>(options);
   }
   else {
      fprintf(stderr,"ERROR: Unsupported Preconditioner: %s\n", PreconName);
      fprintf(stderr,"WARNING: Converting the Preconditioner to AMG\n");
      typedef precon::AMG<PreconMatrixType> precon_type;
      base_precon = _allocate<precon_type>(options);
   }

   if (verbose) {
     int precon_tag = base_precon->precon_tag;
     printf("precon_tag = %d %s\n", precon_tag, precon_name(precon_tag).c_str());
     printf("value_type = %s\n", typeid(ValueType).name());
     printf("Matrix     = %s\n", typeid(Matrix).name());
   }

   return base_precon;
}

// Build the preconditioner and return a pointer to a generic base ...
template <typename Matrix>
precon::BasePrecon* create_precon_0 (char *PreconName)
{
   using namespace precon;

   // Load the default options ... and read the XML file.
   precon::BasePreconOptions options;

   typedef typename Matrix::value_type value_type;

   if (options.mixed_precision and utils::is_same<value_type, float >::value)
   {
      fprintf(stderr,"mixed_precision disabled: both inputs are float()\n");
      options.mixed_precision = false;
   }

   if (options.mixed_precision)
      return create_precon_1<Matrix> (PreconName, options, float());
   else
      return create_precon_1<Matrix> (PreconName, options, value_type());
}

precon::BasePrecon* create_precon (char *PreconName, int format_tag)
{

   typedef double value_type;

   if      (format_tag == splib::FormatTags::CSR)
      return create_precon_0< splib::csr_matrix<value_type> >(PreconName);
   else if (format_tag == splib::FormatTags::COO)
      return create_precon_0< splib::coo_matrix<value_type> >(PreconName);
   else if (format_tag == splib::FormatTags::DIA)
      return create_precon_0< splib::dia_matrix<value_type> >(PreconName);
   else if (format_tag == splib::FormatTags::ELL)
      return create_precon_0< splib::ell_matrix<value_type> >(PreconName);
   else if (format_tag == splib::FormatTags::HYB)
      return create_precon_0< splib::hyb_matrix<value_type> >(PreconName);
   else
   {
      fprintf(stderr,"Error in create_precon: invalid matrix format tag %d\n", format_tag);
      return NULL;
   }
}

// Destroy the preconditioner ...
precon::BasePrecon*
destroy_precon (precon::BasePrecon *baseM)
{
   using namespace precon;

   if (baseM == NULL)
      return NULL;

   if (__DEBUG) printf("destroy_precon: %s %s\n", typeid(baseM).name(), typeid(*baseM).name());

   delete baseM;

   return NULL;
}

splib::BaseMatrix*
destroy_matrix (splib::BaseMatrix *baseA)
{
   using namespace splib;

   if (baseA == NULL)
      return NULL;

   if (__DEBUG) printf("destroy_matrix: %s %s\n", typeid(baseA).name(), typeid(*baseA).name());

   delete baseA;

   return NULL;
}

// Build the solver and return a pointer to a generic base ...
template <typename Matrix>
solver::BaseSolver* create_solver_0 (char *SolverName)
{
   std::string solver_name_comp = SolverName;
   // Turn to lower case to allow for comparison w/ case-insensitive input
   for (unsigned int i = 0; i < solver_name_comp.size(); i++)
       solver_name_comp[i] = tolower(SolverName[i]);

   solver::BaseSolver *base_solver = NULL;
   
   if (solver_name_comp.compare("pcg") == 0) {
      typedef solver::PCG<Matrix> solver_type;
      base_solver = _allocate<solver_type>();
   }
   else if (solver_name_comp.compare("bicgstab") == 0) {
      typedef solver::BiCGstab<Matrix> solver_type;
      base_solver = _allocate<solver_type>();
   }
   else if (solver_name_comp.compare("gmres") == 0) {
      typedef solver::GMRES<Matrix> solver_type;
      base_solver = _allocate<solver_type>();
      base_solver->solver_tag = solver::SolverTags::GMRES;
   }
   else if (solver_name_comp.compare("fgmres") == 0) {
      typedef solver::GMRES<Matrix> solver_type;
      base_solver = _allocate<solver_type>();
      base_solver->solver_tag = solver::SolverTags::FGMRES;
   }
   else if (solver_name_comp.compare("rgmres") == 0) {
      typedef solver::GMRES<Matrix> solver_type;
      base_solver = _allocate<solver_type>();
      base_solver->solver_tag = solver::SolverTags::RGMRES;
   }
   else {
      fprintf(stderr,"ERROR: Unsupported Solver: %s\n", SolverName);
      fprintf(stderr,"WARNING: Converting the Solver to FGMRES\n");
      typedef solver::GMRES<Matrix> solver_type;
      base_solver = _allocate<solver_type>();
      base_solver->solver_tag = solver::SolverTags::FGMRES;
   }

   if (verbose) {
     int solver_tag = base_solver->solver_tag;
     printf("solver_tag = %d %s\n", solver_tag, solver::solver_name(solver_tag).c_str());
   }

   return base_solver;
}

solver::BaseSolver* create_solver (char *SolverName, int format_tag)
{

   typedef double value_type;

   if      (format_tag == splib::FormatTags::CSR)
      return create_solver_0< splib::csr_matrix<value_type> >(SolverName);
   else if (format_tag == splib::FormatTags::COO)
      return create_solver_0< splib::coo_matrix<value_type> >(SolverName);
   else if (format_tag == splib::FormatTags::DIA)
      return create_solver_0< splib::dia_matrix<value_type> >(SolverName);
   else if (format_tag == splib::FormatTags::ELL)
      return create_solver_0< splib::ell_matrix<value_type> >(SolverName);
   else if (format_tag == splib::FormatTags::HYB)
      return create_solver_0< splib::hyb_matrix<value_type> >(SolverName);
   else
   {
      fprintf(stderr,"Error in create_solver: invalid matrix format tag %d\n", format_tag);
      return NULL;
   }
}

// Destroy the solver ...
solver::BaseSolver*
destroy_solver (solver::BaseSolver *base)
{
   if (base == NULL)
      return NULL;

   if (__DEBUG) printf("destroy_solver: %s %s\n", typeid(base).name(), typeid(*base).name());

   delete base;

   return NULL;
}

extern "C"
{

#include <aplles_interface.h>

int APLLES_Setup_Matrix_CSR (int nrows, int rowptr[], int colidx[], double values[], APLLES_MatrixHandle_t &A_handle)
{
   // Assign the matrix type; csr
   typedef double value_type;
   typedef splib::csr_matrix<value_type> Matrix;

   // Initialize Matrix A
   Matrix *A_ptr = new Matrix();

   // Assign passed on values to Matrix A
   int nnz = rowptr[nrows];
   A_ptr->set_view(nrows,nrows,nnz,rowptr,colidx,values);

   // Need the following casting for the Fortran interface
   A_handle = reinterpret_cast<APLLES_MatrixHandle_t>(A_ptr);

   if (verbose) printf("APLLES_Setup_Matrix_CSR: %p, %d %d\n", A_ptr, nrows, nnz);

   return 0;
}

int APLLES_Setup_Solver (APLLES_MatrixHandle_t &A_handle, char *Solver_Name, APLLES_SolverHandle_t & S_handle)
{
   // Cast A_handle back onto a Base Matrix ...
   typedef splib::BaseMatrix Matrix;
   Matrix *baseA = reinterpret_cast<Matrix *>(A_handle);

   // Create the solver given a Solver_Name type...
   solver::BaseSolver *baseSolver = create_solver(Solver_Name, baseA->format_tag());

   if (verbose) printf("baseSolver->name() = %s %d\n", baseSolver->name().c_str(), baseSolver->tag());

   // Build the solver
   int iret = baseSolver->build(baseA);

   // Cast baseSolver onto a S_handle for portability
   S_handle = reinterpret_cast<APLLES_SolverHandle_t>(baseSolver);

   if (verbose) printf("APLLES_Setup_Solver: baseS = %p, Solver_Name = %s\n", baseSolver, Solver_Name);

   if (iret) {
     printf("Error building the solver ...");
     return iret;
   }

   return 0;
}

int APLLES_Setup_Precon (APLLES_MatrixHandle_t &M_handle, char *Precon_Name, APLLES_PreconHandle_t &P_handle)
{
   // Cast M_handle back onto a Base Matrix ...
   typedef splib::BaseMatrix Matrix;
   Matrix *baseM = reinterpret_cast<Matrix *>(M_handle);

   // Create the precon given a Precon_Name type ...
   precon::BasePrecon *basePrecon = create_precon(Precon_Name, baseM->format_tag());

   if (verbose) printf("basePrecon->name() = %s %d\n", basePrecon->name().c_str(), basePrecon->tag());

   // Build the precon
   int iret = basePrecon->build(baseM);

   // Cast basePrecon onto a P_handle for portability
   P_handle = reinterpret_cast<APLLES_PreconHandle_t>(basePrecon);

   if (verbose) printf("APLLES_Setup_Precon: baseM = %p, Precon_Name = %s\n", basePrecon, Precon_Name);

   if (iret) {
     printf("Error building the preconditioner ...");
     return iret;
   }

   return 0;
}

int APLLES_Solve (APLLES_MatrixHandle_t &A_handle, double x_ptr[], double b_ptr[], APLLES_SolverHandle_t &S_handle, APLLES_PreconHandle_t &P_handle)
{
   splib ::BaseMatrix *baseA      = reinterpret_cast<splib ::BaseMatrix *>(A_handle);
   solver::BaseSolver *baseSolver = reinterpret_cast<solver::BaseSolver *>(S_handle);
   precon::BasePrecon *basePrecon = reinterpret_cast<precon::BasePrecon *>(P_handle);

   typedef double value_type;
   typedef splib::Vector<value_type> Vector;

   Vector b, x;
    
   // Initialize the unknown
   int nrows = baseA->num_rows;
   b.set_view(b_ptr,nrows);
   x.set_view(x_ptr,nrows);

   // Solve the problem
   int iret = baseSolver->solve (baseA, &x, &b, basePrecon);

   return iret;
}

int APLLES_Destroy_Matrix (APLLES_MatrixHandle_t &A_handle)
{
   if (A_handle == NULL)
      return 1;

   splib ::BaseMatrix *baseA = reinterpret_cast<splib::BaseMatrix *>(A_handle);

   if (verbose) printf("APLLES_Destroy_Matrix: baseA = %p\n", baseA);

   if (baseA == NULL)
      return 2;

   baseA = _destroy(baseA);

   A_handle = reinterpret_cast<APLLES_MatrixHandle_t>(baseA);

   return 0;
}
int APLLES_Destroy_Precon (APLLES_PreconHandle_t &M_handle)
{
   if (M_handle == NULL)
      return 1;

   precon::BasePrecon *baseM = reinterpret_cast<precon::BasePrecon *>(M_handle);

   if (verbose) printf("APLLES_Destroy_Precon: baseM = %p\n", baseM);

   if (baseM == NULL)
      return 2;

   baseM = _destroy(baseM);

   M_handle = reinterpret_cast<APLLES_PreconHandle_t>(baseM);

   return 0;
}
int APLLES_Destroy_Solver (APLLES_SolverHandle_t &S_handle)
{
   if (S_handle == NULL)
      return 1;

   solver::BaseSolver *baseS = reinterpret_cast<solver::BaseSolver *>(S_handle);

   if (verbose) printf("APLLES_Destroy_Solver: baseS = %p\n", baseS);

   if (baseS == NULL)
      return 2;

   baseS = _destroy(baseS);

   S_handle = reinterpret_cast<APLLES_SolverHandle_t>(baseS);

   return 0;
}

int APLLES_Initialize(void)
{
   parallel::initialize();

   return 0;
}

int APLLES_Terminate(void)
{
   return 0;
}

int APLLES_Matrix_Copy_CSR_To_Other (APLLES_MatrixHandle_t &CSRMatrix_handle, char *format_name, APLLES_MatrixHandle_t &OtherMatrix_handle)
{
   std::string format_name_comp = format_name;
   // Turn to lower case to allow for comparison w/ case-insensitive input
   for (unsigned int i = 0; i < format_name_comp.size(); i++)
       format_name_comp[i] = tolower(format_name[i]);

   typedef splib::csr_matrix<double> csr_matrix;
   csr_matrix *csr = reinterpret_cast<csr_matrix *>(CSRMatrix_handle);

   if (format_name_comp.compare("ell") == 0) {
      typedef splib::ell_matrix<double> Matrix;
      Matrix *othr = new Matrix();
      splib::matrix_copy(*csr, *othr);
      OtherMatrix_handle = reinterpret_cast<APLLES_MatrixHandle_t>(othr);
   }
   else if (format_name_comp.compare("coo") == 0) {
      typedef splib::coo_matrix<double> Matrix;
      Matrix *othr = new Matrix();
      splib::matrix_copy(*csr, *othr);
      OtherMatrix_handle = reinterpret_cast<APLLES_MatrixHandle_t>(othr);
   }
   else if (format_name_comp.compare("hyb") == 0) {
      typedef splib::hyb_matrix<double> Matrix;
      Matrix *othr = new Matrix();
      splib::matrix_copy(*csr, *othr);
      OtherMatrix_handle = reinterpret_cast<APLLES_MatrixHandle_t>(othr);
   }
   else if (format_name_comp.compare("dia") == 0) {
      typedef splib::dia_matrix<double> Matrix;
      Matrix *othr = new Matrix();
      splib::matrix_copy(*csr, *othr);
      OtherMatrix_handle = reinterpret_cast<APLLES_MatrixHandle_t>(othr);
   }
   else if (format_name_comp.compare("csr") == 0) {
      //fprintf(stderr,"WARNING: Unnecessary CSR matrix copy operation\n");
      //OtherMatrix_handle = CSRMatrix_handle;
      //return 1;
      //fprintf(stderr,"WARNING: Copy CSR -> CSR\n");
      typedef splib::csr_matrix<double> Matrix;
      Matrix *othr = new Matrix();
      splib::matrix_copy(*csr, *othr);
      OtherMatrix_handle = reinterpret_cast<APLLES_MatrixHandle_t>(othr);
   }
   else {
      fprintf(stderr,"ERROR: Invalid matrix format %s\n", format_name);
      fprintf(stderr,"WARNING: Returning CSR matrix format\n");
      OtherMatrix_handle = CSRMatrix_handle;
      return 1;
   }

   return 0;

}

} // extern "C"
