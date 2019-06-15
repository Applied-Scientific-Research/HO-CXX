
#ifdef __cplusplus
extern "C"
{
#endif

/* Prototypes for the BASIC Aplles interface */

struct __APLLES_MatrixHandle_t {};
struct __APLLES_SolverHandle_t {};
struct __APLLES_PreconHandle_t {};

typedef __APLLES_MatrixHandle_t* APLLES_MatrixHandle_t;
typedef __APLLES_SolverHandle_t* APLLES_SolverHandle_t;
typedef __APLLES_PreconHandle_t* APLLES_PreconHandle_t;

int APLLES_Setup_Matrix_CSR (int nrows, int rowptr[], int colidx[], double values[], APLLES_MatrixHandle_t &A_handle);
int APLLES_Setup_Solver (APLLES_MatrixHandle_t &A_handle, char *Solver_Name, APLLES_SolverHandle_t &S_handle);
int APLLES_Setup_Precon (APLLES_MatrixHandle_t &M_handle, char *Precon_Name, APLLES_PreconHandle_t &P_handle);
int APLLES_Solve (APLLES_MatrixHandle_t &A_handle, double x[], double b[], APLLES_SolverHandle_t &S_handle, APLLES_PreconHandle_t &P_handle);
int APLLES_Destroy_Matrix (APLLES_MatrixHandle_t &A_handle);
int APLLES_Destroy_Precon (APLLES_PreconHandle_t &A_handle);
int APLLES_Destroy_Solver (APLLES_SolverHandle_t & A_handle);
int APLLES_Initialize(void);
int APLLES_Terminate(void);
int APLLES_Matrix_Copy_CSR_To_Other (APLLES_MatrixHandle_t &CSRMatrix_handle, char *format_name, APLLES_MatrixHandle_t &OtherMatrix_handle);

#ifdef __cplusplus
}
#endif
