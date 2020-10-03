module APLLES_Solvers_Module

  use iso_c_binding

  implicit NONE

  public
  TYPE APLLES_MatrixHandleType
     TYPE(C_PTR) :: ptr
  END TYPE
  TYPE APLLES_SolverHandleType
     TYPE(C_PTR) :: ptr
  END TYPE
  TYPE APLLES_PreconHandleType
     TYPE(C_PTR) :: ptr
  END TYPE

! private
  public APLLES_Initialize, APLLES_Terminate
  public APLLES_Setup_Matrix_CSR, APLLES_Setup_Solver, APLLES_Setup_Precon, APLLES_Solve
  public APLLES_Destroy_Matrix, APLLES_Destroy_Solver, APLLES_Destroy_Precon
  public APLLES_Matrix_Copy_CSR_To_Other

  interface
    function APLLES_Initialize() bind(C, NAME='APLLES_Initialize')
      use iso_c_binding
      implicit NONE
      integer(C_INT):: APLLES_Initialize
    end function APLLES_Initialize

    function APLLES_Terminate() bind(C, NAME='APLLES_Terminate')
      use iso_c_binding
      implicit NONE
      integer(C_INT):: APLLES_Terminate
    end function APLLES_Terminate

    function APLLES_Setup_Matrix_CSR_C(nrows,rowptr,colidx,values,A_ptr) &
                                  bind(C, NAME='APLLES_Setup_Matrix_CSR')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Setup_Matrix_CSR_C
      integer(C_INT), value :: nrows
      integer(C_INT), dimension(*), intent(in):: rowptr, colidx
      real(C_DOUBLE), dimension(*), intent(in) :: values
      type(C_PTR), intent(inout) :: A_ptr
    end function APLLES_Setup_Matrix_CSR_C

    function APLLES_Destroy_Matrix_C(A_ptr) &
                                  bind(C, NAME='APLLES_Destroy_Matrix')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Destroy_Matrix_C
      type(C_PTR), intent(inout) :: A_ptr
    end function APLLES_Destroy_Matrix_C

    function APLLES_Destroy_Solver_C(S_ptr) &
                                  bind(C, NAME='APLLES_Destroy_Solver')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Destroy_Solver_C
      type(C_PTR), intent(inout) :: S_ptr
    end function APLLES_Destroy_Solver_C

    function APLLES_Destroy_Precon_C(P_ptr) &
                                  bind(C, NAME='APLLES_Destroy_Precon')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Destroy_Precon_C
      type(C_PTR), intent(inout) :: P_ptr
    end function APLLES_Destroy_Precon_C

    function APLLES_Setup_Solver_C(A_ptr,Solver_Name,S_ptr) &
                                bind(C, NAME='APLLES_Setup_Solver')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Setup_Solver_C
      type(C_PTR), intent(in) :: A_ptr
      type(C_PTR), intent(inout) :: S_ptr
      character(C_CHAR), dimension(*), intent(in) :: Solver_Name
    end function APLLES_Setup_Solver_C

    function APLLES_Setup_Precon_C (A_ptr, Precon_Name, P_ptr) &
                                bind(C, NAME='APLLES_Setup_Precon')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Setup_Precon_C
      type(C_PTR), intent(in) :: A_ptr
      type(C_PTR), intent(inout) :: P_ptr
      character(C_CHAR), dimension(*), intent(in) :: Precon_Name
    end function APLLES_Setup_Precon_C

    function APLLES_Solve_C (A_ptr, x, b, S_ptr, P_ptr) &
                       bind(C, NAME='APLLES_Solve')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Solve_C
      type(C_PTR), intent(in) :: A_ptr
      type(C_PTR), intent(in) :: S_ptr
      type(C_PTR), intent(in) :: P_ptr
      real(KIND=C_DOUBLE), dimension(*), intent(inout) :: x
      real(KIND=C_DOUBLE), dimension(*), intent(in) :: b
    end function APLLES_Solve_C

    function APLLES_Matrix_Copy_CSR_To_Other_C (CSR_ptr, format_name, Other_ptr) &
                                  bind(C, NAME='APLLES_Matrix_Copy_CSR_To_Other')
      use iso_c_binding
      implicit NONE
      integer(C_INT) :: APLLES_Matrix_Copy_CSR_To_Other_C
      type(C_PTR), intent(in) :: CSR_ptr
      type(C_PTR), intent(inout) :: Other_ptr
      character(C_CHAR), dimension(*), intent(in) :: format_name
    end function APLLES_Matrix_Copy_CSR_To_Other_C

  end interface

contains

  function APLLES_Setup_Matrix_CSR (nrows, rowptr, colidx, values, A_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    integer :: nrows
    integer, dimension(*), intent(in):: rowptr, colidx
    real(C_DOUBLE), dimension(*), intent(in) :: values
    type(APLLES_MatrixHandleType), intent(inout) :: A_handle

    integer :: ierr

    ierr = APLLES_Setup_Matrix_CSR_C (nrows, rowptr, colidx, values, A_handle%ptr)

  end function APLLES_Setup_Matrix_CSR

  function APLLES_Setup_Solver (A_handle, Solver_Name, S_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_MatrixHandleType), intent(in) :: A_handle
    type(APLLES_SolverHandleType), intent(inout) :: S_handle
    character(kind=C_CHAR,len=*), intent(in) :: Solver_Name

    integer :: ierr

    ierr = APLLES_Setup_Solver_C (A_handle%ptr, trim(Solver_Name)//C_NULL_CHAR, S_handle%ptr)

  end function APLLES_Setup_Solver

  function APLLES_Setup_Precon (M_handle, Precon_Name, P_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_MatrixHandleType), intent(in) :: M_handle
    type(APLLES_PreconHandleType), intent(inout) :: P_handle
    character(kind=C_CHAR,len=*), intent(in) :: Precon_Name

    integer :: ierr

    ierr = APLLES_Setup_Precon_C (M_handle%ptr, trim(Precon_Name)//C_NULL_CHAR, P_handle%ptr)

  end function APLLES_Setup_Precon

  function APLLES_Solve (A_handle, x, b, S_handle, P_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_MatrixHandleType), intent(in) :: A_handle
    real(KIND=C_DOUBLE), dimension(*), intent(inout) :: x
    real(KIND=C_DOUBLE), dimension(*), intent(in) :: b
    type(APLLES_SolverHandleType), intent(inout) :: S_handle
    type(APLLES_PreconHandleType), intent(inout) :: P_handle

    integer :: ierr

    ierr = APLLES_Solve_C (A_handle%ptr, x, b, S_handle%ptr, P_handle%ptr)

  end function APLLES_Solve

  function APLLES_Destroy_Precon (P_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_PreconHandleType), intent(inout) :: P_handle

    integer :: ierr

    ierr = APLLES_Destroy_Precon_C (P_handle%ptr)

  end function APLLES_Destroy_Precon

  function APLLES_Destroy_Solver (S_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_SolverHandleType), intent(inout) :: S_handle

    integer :: ierr

    ierr = APLLES_Destroy_Solver_C (S_handle%ptr)

  end function APLLES_Destroy_Solver

  function APLLES_Destroy_Matrix (A_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_MatrixHandleType), intent(inout) :: A_handle

    integer :: ierr

    ierr = APLLES_Destroy_Matrix_C (A_handle%ptr)

  end function APLLES_Destroy_Matrix

  function APLLES_Matrix_Copy_CSR_To_Other (csr_handle, format_name, A_handle) result(ierr)

    use iso_c_binding
    implicit NONE

! Input/output data    
    type(APLLES_MatrixHandleType), intent(in) :: csr_handle
    type(APLLES_MatrixHandleType), intent(inout) :: A_handle
    character(kind=C_CHAR,len=*), intent(in) :: format_name

    integer :: ierr

    ierr = APLLES_Matrix_Copy_CSR_To_other_C (csr_handle%ptr, trim(format_name)//C_NULL_CHAR, A_handle%ptr)

  end function APLLES_Matrix_Copy_CSR_To_Other

end module APLLES_Solvers_Module
