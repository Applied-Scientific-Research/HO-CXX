! High-order diffusion solver in 2D based on Methods 2 and 11 in Huynh's 2009 AIAA paper*

! * H.T. Huynh, AIAA2009-403
!       "A Reconstruction Approach to High-Order Schemes Including Discontinuous Galerkin for Diffusion"

! The solver can only accept quadrilateral meshes by design
! Solver orders 1 thru 5  (same order in x and y directions)
! Quad mesh orders 1 and 2  (same order in x and y directions)
! Explicit time integration orders 1, RK2, and RK4

! Benchmark problems implemented:
! (1) Square Domain: Dirichlet (u = 0) at East (E) & West (W); Neumann (n.du = 0) at North (N) & South (S)
! (2) Square Domain: u = 0 at E; n.du = 0 elsewhere
! (3) Square Domain: u = 0 at NEWS
! (4) Concentric Cylinders: R_i = 0.5; R_o = 1.0; BC: u = 0 at N & S

!
! By: Adrin Gharakhani, Sc.D.
!     Applied Scientific Research, Irvine, CA
!     Date: 10-15-17
!
! Project supported by NIH Grant 1 R01 EB022180-01A1

MODULE params
       implicit NONE
       real*8, parameter :: pi = 3.1415926535897931d0
       integer, parameter :: NDIM = 3
       integer, parameter :: FC = 2*NDIM - 1
       integer, parameter :: DefaultBC = 0, DirichletBC = 1, NeumannBC = 2

       integer Knod, Lnod, Nel, Nodes, NelB, NBndry
       ! i2f converts indexed notation (0-1 are left-right; 2-3 are south-north) to CCW FEM face numbers starting from the south
       ! face
!       integer, parameter, dimension(0:3) :: i2f = (/4, 2, 1, 3/), nbr = (/1, 0, 3, 2/)
       ! Most likely can just use the follwowing 3D for 2D as well
       integer, parameter, dimension(0:5) :: i2f = (/4, 2, 1, 3, 5, 6/), nbr = (/1, 0, 3, 2, 5, 4/)

       logical exact

END MODULE params


MODULE variables
       implicit NONE
                                              
       real*8, allocatable :: xcoord(:), ycoord(:), zcoord(:)  ! nodal positions of the geometry elements (structured elements for now)
       integer, allocatable :: nodeID(:,:)          ! quad-element node indices, starting from SW node and going CCW 1-2-3-4 (5-6-7-8)
       integer, allocatable :: elemID(:,:)          ! face connectivity; elemID(nf,el) points to a neighboring mesh/element that shares
                                                    ! face nf of mesh el
                                                    ! elemID(:,:) < 0 points to boundary meshes, where Abs[ elemID(:,:)] is the mesh number
       integer, allocatable :: faceID(:,:)          ! face connectivity; faceID(nf,el) points to the face number of elemID(nf,el) that shares
                                                    ! face nf of mesh el
       integer, allocatable :: BoundaryNodeID(:,:)  ! Boundary face node indices
                                                    ! For now it is the same as nodeID for those elements; i.e., it points to the same coords
                                                    ! BUT the general version would be completely independent, and linked to quads via elemID
                                                    ! WARNING: We MUST make sure the numbering order is consistent between nodeID and BoundaryNodeID
       integer, allocatable :: boundaryPointsElemID(:) ! index pointing to mesh element number containing that boundary element
       character(len=100), allocatable :: BndryTag(:)  ! Each boundary will be given a name/tag
       logical, allocatable :: NoSlip(:)            ! Assign NoSlip (= .true.), or Slip (.false.) for each tagged boundary patch
       integer, allocatable :: BndryNum(:)          ! Cumulative number of boundary elements; number of elements in each patch i is
                                                    ! BndryNum(i) - BndryNum(i-1),  where BndryNum(0)
       real*8, allocatable :: VelocJump(:,:,:)
       real*8, allocatable :: gps(:)                ! scaled/parametric coordinates of the geometric nodes (in 1-D)
       real*8, allocatable :: sps(:)                ! scaled/parametric coordinates of the solution points (in 1-D)
       real*8, allocatable :: wgt(:)                ! Gaussian weights at the solution points (in 1-D)

       real*8, allocatable, dimension(:,:,:,:) :: Vort0, Vort   ! initial Vort0(i,kx,ky) and convected+diffused Vort(i,kx,ky) values of the unknowns
                                                    ! per element i and at (kx,ky) nodes per elem
       real*8, allocatable, dimension(:,:,:,:) :: Uvel, Vvel, Wvel    ! u(i,kx,ky) and v(i,kx,ky) velocity components
                                                    ! per element i and at (kx,ky) nodes per elem
       real*8, allocatable :: Vol_Dx_iDxsi_j(:,:,:,:,:,:), Vol_Jac(:,:,:,:)
       real*8, allocatable, dimension(:,:,:,:) :: Face_Acoef, Face_Bcoef, Face_Ccoef, Face_Norm, Face_Jac
       ! These bases use the element nodes as collocation points
       real*8, allocatable :: GeomBndryLgrangeBasis(:,:)       ! Lagrange interpolation basis for Left and Right boundaries of elem
       real*8, allocatable :: GeomBndryGradLgrangeBasis(:,:)   ! Corresponding Basis for derivatives at L + R boundaries
       real*8, allocatable :: GeomNodesLgrangeBasis(:,:)       ! Basis for Lagrange interpol. at solution points
       real*8, allocatable :: GeomNodesGradLgrangeBasis(:,:)   ! Basis for derivatives of Lagrange interpol. at solution points
       ! These bases use the solution nodes as collocation points
       real*8, allocatable :: SolnBndryLgrangeBasis(:,:)       ! Lagrange interpolation basis for Left and Right boundaries of elem
       real*8, allocatable :: SolnBndryGradLgrangeBasis(:,:)   ! Corresponding Basis for derivatives at L + R boundaries
       real*8, allocatable :: SolnNodesGradLgrangeBasis(:,:)   ! Basis for derivatives of Lagrange interpol. at solution points
       real*8, allocatable :: NodesRadau(:,:)                  ! Left (0) and Right (1) Radau function at solution points
       real*8, allocatable :: NodesGradRadau(:,:)              ! Derivatives of Left (0) and Right (1) Radau function at solution points
       real*8 gLprime(0:1)
       ! Neumann or Dirichlet BC. The patches are numbered from the bottom/south and go counter-clockwise
       ! S=1; E=2; N=3; W=4
       ! For now, memory is over-allocated based on the original value of NelB (NelB gets modified properly in MeshSetup)
       real*8, allocatable :: BC_VelNorm(:,:,:), BC_VelParl(:,:,:), BC_Psi(:,:,:)
       real*8, allocatable ::  BC_Values(:,:,:)                  ! BC values for the Knod points of NelB boundary elements
       integer, allocatable ::  BC_Switch(:), BC_Switch_Lap(:)                  ! The BC switch for these patches is Dirichlet = 1; Neumann = 2

       real*8, allocatable :: BndrySrc(:,:,:,:)                  ! Poisson Solver's RHS term to be dotted by BC_Values

       character(len=30) :: aplles_solver_name, aplles_precon_name

END MODULE variables


MODULE CnvrtTensor2FemIndices

       ! We use tensor product of 1-D arrays, which is incosistent with FEM index notation
       ! This index conversion is the easiest way to handle this, especially when dealing
       ! w/ higher order elements.  BUT, we have to manually create the indices for each case

       ! index t2f(i,j,k) converts node (i,j,k) to traditional FEM nodes, or other FEM formats (to be included)

       USE params

       PRIVATE

       INTERFACE T2F
         MODULE PROCEDURE T2F1D, T2F2D, T2F3D
       END INTERFACE T2F
       PUBLIC T2F

CONTAINS

  FUNCTION T2F1D(i)
       implicit NONE

       integer t2f1d
       integer i

       IF (Lnod .eq. 2) THEN

         IF (i .eq. 1) t2f1d = 1
         IF (i .eq. 2) t2f1d = 2

       ELSEIF (Lnod .eq. 3) THEN

         IF (i .eq. 1) t2f1d = 1
         IF (i .eq. 2) t2f1d = 3
         IF (i .eq. 3) t2f1d = 2

       ELSEIF (Lnod .eq. 4) THEN

         IF (i .eq. 1) t2f1d = 1
         IF (i .eq. 2) t2f1d = 3
         IF (i .eq. 3) t2f1d = 4
         IF (i .eq. 4) t2f1d = 2

       ELSE
         print *,'Only up to cubic finite elements (Lnod = 4) are supported!'
         STOP
       ENDIF
  END FUNCTION T2F1D

  FUNCTION T2F2D(i,j)
       implicit NONE

       integer t2f2d
       integer i, j

       IF (Lnod .eq. 2) THEN

         IF (i .eq. 1 .and. j .eq. 1) t2f2d = 1
         IF (i .eq. 2 .and. j .eq. 1) t2f2d = 2
         IF (i .eq. 1 .and. j .eq. 2) t2f2d = 4
         IF (i .eq. 2 .and. j .eq. 2) t2f2d = 3

       ELSEIF (Lnod .eq. 3) THEN

         IF (i .eq. 1 .and. j .eq. 1) t2f2d = 1
         IF (i .eq. 2 .and. j .eq. 1) t2f2d = 5
         IF (i .eq. 3 .and. j .eq. 1) t2f2d = 2
         IF (i .eq. 1 .and. j .eq. 2) t2f2d = 8
         IF (i .eq. 2 .and. j .eq. 2) t2f2d = 9
         IF (i .eq. 3 .and. j .eq. 2) t2f2d = 6
         IF (i .eq. 1 .and. j .eq. 3) t2f2d = 4
         IF (i .eq. 2 .and. j .eq. 3) t2f2d = 7
         IF (i .eq. 3 .and. j .eq. 3) t2f2d = 3

       ELSEIF (Lnod .eq. 4) THEN

         IF (i .eq. 1 .and. j .eq. 1) t2f2d = 1
         IF (i .eq. 2 .and. j .eq. 1) t2f2d = 5
         IF (i .eq. 3 .and. j .eq. 1) t2f2d = 6
         IF (i .eq. 4 .and. j .eq. 1) t2f2d = 2
         IF (i .eq. 1 .and. j .eq. 2) t2f2d = 12
         IF (i .eq. 2 .and. j .eq. 2) t2f2d = 13
         IF (i .eq. 3 .and. j .eq. 2) t2f2d = 14
         IF (i .eq. 4 .and. j .eq. 2) t2f2d = 7
         IF (i .eq. 1 .and. j .eq. 3) t2f2d = 11
         IF (i .eq. 2 .and. j .eq. 3) t2f2d = 16
         IF (i .eq. 3 .and. j .eq. 3) t2f2d = 15
         IF (i .eq. 4 .and. j .eq. 3) t2f2d = 8
         IF (i .eq. 1 .and. j .eq. 4) t2f2d = 4
         IF (i .eq. 2 .and. j .eq. 4) t2f2d = 10
         IF (i .eq. 3 .and. j .eq. 4) t2f2d = 9
         IF (i .eq. 4 .and. j .eq. 4) t2f2d = 3

       ELSE
         print *,'Only up to cubic finite elements (Lnod = 4) are supported!'
         STOP
       ENDIF
  END FUNCTION T2F2D

  FUNCTION T2F3D(i,j,k)
       implicit NONE

       integer t2f3d
       integer i, j, k

       IF (Lnod .eq. 2) THEN

         IF (i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) t2f3d = 1
         IF (i .eq. 2 .and. j .eq. 1 .and. k .eq. 1) t2f3d = 2
         IF (i .eq. 1 .and. j .eq. 2 .and. k .eq. 1) t2f3d = 4
         IF (i .eq. 2 .and. j .eq. 2 .and. k .eq. 1) t2f3d = 3

         IF (i .eq. 1 .and. j .eq. 1 .and. k .eq. 2) t2f3d = 5
         IF (i .eq. 2 .and. j .eq. 1 .and. k .eq. 2) t2f3d = 6
         IF (i .eq. 1 .and. j .eq. 2 .and. k .eq. 2) t2f3d = 8
         IF (i .eq. 2 .and. j .eq. 2 .and. k .eq. 2) t2f3d = 7

       ELSE
         print *,'Only up to linear finite elements (Lnod = 2) are supported!'
         STOP
       ENDIF
  END FUNCTION T2F3D

END MODULE CnvrtTensor2FemIndices


PROGRAM ThreeD_Vorticity_Transport
       USE params
       USE variables
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer Nelx, Nely, Nelz, Lnod_in, prob_type, HuynhSolver_type, tIntegrator_type, numStep, dumpFreq, fast
       real*8 Reyn, fac, dxrat, dt
       integer i, el, idum, ic
       character(132) dum, mesh_filename

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle
 
       !!!! User Input Data
       open(unit=2, file='input.dat', status='old')
         read(2,'(a)') dum
         read(2,*) Nelx, Nely, Nelz     ! num. of meshes/elements in x, y, z dirs (structured grids; must generalize to unstructured)
         read(2,'(a)') dum
         read(2,*) Knod                 ! solution nodes per elem (must generalize to include x, y, z dirs)
         read(2,'(a)') dum
         read(2,*) Lnod_in              ! geometry parametric order per elem (must generalize to include x, y, z dirs)
         read(2,'(a)') dum
         read(2,'(a)') mesh_filename    ! SU-based mesh file; it's a hack; using prob_type >= 10 for now
         read(2,'(a)') dum
         read(2,'(L)') exact            ! Exact geometry for concentric cylinders
         read(2,'(a)') dum
         read(2,*) Reyn                 ! Reynolds number
         read(2,'(a)') dum
         read(2,*) dxrat                ! multiplier to expand the grid geometrically by factor dxrat   
         read(2,'(a)') dum
         read(2,*) fac                  ! multiplier to make grids randomly non-uniform; uniform: fac=0.
         read(2,'(a)') dum
         read(2,*) HuynhSolver_type     ! based on Huynh's scheme type in Table 6.1 of his diffusion paper
                                        ! types: 2=std DG; 11=higher order
!         HuynhSolver_type = 2

         read(2,'(a)') dum
         read(2,*) tIntegrator_type     ! time integration method; 1=Euler; 2=Mod. Euler; 4=RK4
         read(2,'(a)') dum
         read(2,*) prob_type            ! Solve different problems (types 1 and 2)
         IF (prob_type /= 1 .and. prob_type /= 3 .and. prob_type /= 10) THEN
           print *,'problem type unavailable ', prob_type
!           stop
         ENDIF
         IF (prob_type .lt. 7) THEN
          exact = .false.
         ENDIF

         read(2,'(a)') dum
         read(2,*) numStep              ! num. of timesteps
         read(2,'(a)') dum
         read(2,*) dt                   ! timestep size
         read(2,'(a)') dum
         read(2,*) dumpFreq             ! dump solution and stats at every dumpFreq steps; = 0 dumps at numStep
         IF (dumpFreq .eq. 0) dumpFreq = numStep

         read(2,'(a)') dum
         read(2,*) fast                 ! enter 0 to solve using original; anything else to solve using compact
         fast = 0

         read(2,'(a)') dum
         read(2,'(a)') aplles_solver_name
         read(2,'(a)') aplles_precon_name ! Select the solver/preconditioner at run-time       
         if (trim(Aplles_solver_name) .eq. "default") Aplles_solver_name = "fgmres"
         if (trim(Aplles_precon_name) .eq. "default") Aplles_precon_name = "amg"
         print *, 'Aplles solver/precon: ', trim(Aplles_solver_name), ' / ', trim(Aplles_precon_name)

       close(unit=2)

       Lnod = Lnod_in + 1
       ! Assuming a Cubic Domain (or Square in 2D)
       NBndry = 4
       IF (NDIM .eq. 2) THEN
         NBndry = 4
         Nelz = 1
       ENDIF
       Nel = Nelx * Nely * Nelz
       NelB = 2 * ( (Nelx + Nely) * Nelz + (NDIM - 2) * Nelx * Nely )
       Nodes = (Nelx*Lnod_in + 1) * (Nely*Lnod_in + 1) * ((NDIM - 2)*Nelz*Lnod_in + 1)

 !***** Do this for 3D *****
       ! Temp hack for generic mesh generation
       IF (prob_type .eq. 10) THEN
         OPEN (10, File = trim(mesh_filename), Status = 'unknown')

         NelB = 0
         read(10,'(a6,i10)') dum, idum
         read(10,'(a6,i10)') dum, Nel
         DO el = 1, Nel
           read(10,'(a)') dum
         ENDDO
         read(10,'(a6,i10)') dum, Nodes
         DO i = 1, Nodes
           read(10,'(a)') dum
         ENDDO
         ! *******
         ! See the relevant subroutine for explanation of the following
!         Nodes = Nodes - 1
         ! *******
         read(10,'(a6,i10)') dum, NBndry
         DO i = 1, NBndry
           read(10,'(a)') dum
           read(10,'(a14,i10)') dum, idum
           NelB = NelB + idum
           DO el = 1, idum
             read(10,'(a)') dum
           ENDDO
         ENDDO

         CLOSE (10)
         print *,'Nel= ',Nel,' Nodes= ',Nodes,' NelB= ',NelB,' NBndry= ',NBndry
       ENDIF

       ! Nothing to do here: Just allocate arrays
       CALL Allocate_Arrays(Nel, Nodes, NBndry, NelB, Knod, Lnod)

       IF (Knod .eq. 1) THEN
         sps(1) = 0.d0                                       ! this case doesn't work; I'd need to do special cases
         wgt(1) = 2.d0                                       ! Gaussian weights
       ELSEIF (Knod .eq. 2) THEN
         sps(2) = 1.d0/sqrt(3.0d0)                           ! collocation nodes per elem are {-C1, C1}
         sps(1) = -sps(2)
         wgt(2) = 1.d0                                       ! Gaussian weights
         wgt(1) = wgt(2)
       ELSEIF (Knod .eq. 3) THEN
         sps(3) = sqrt(0.6d0)                                ! collocation nodes per elem are {-C1, 0, C1}
         sps(2) = 0.d0
         sps(1) = -sps(3)                                    ! sloppy but quick; we really don't need to save the -ve values
         wgt(3) = 5.d0/9.d0                                  ! Gaussian weights
         wgt(2) = 8.d0/9.d0
         wgt(1) = wgt(3)
       ELSEIF (Knod .eq. 4) THEN
         sps(4) = sqrt((15.d0+2.d0*sqrt(30.d0))/35.d0)       ! collocation nodes per elem are {-C2, -C1, C1, C2}
         sps(3) = sqrt((15.d0-2.d0*sqrt(30.d0))/35.d0)
         sps(2) = -sps(3)
         sps(1) = -sps(4)
         wgt(4) = (18.d0-sqrt(30.d0))/36.d0                  ! Gaussian weights
         wgt(3) = (18.d0+sqrt(30.d0))/36.d0
         wgt(2) = wgt(3)
         wgt(1) = wgt(4)
       ELSEIF (Knod .eq. 5) THEN
         sps(5) = sqrt(5.d0+2.d0*sqrt(70.d0)/7.d0)/3.d0      ! collocation nodes per elem are {-C2, -C1, 0, C1, C2}
         sps(4) = sqrt(5.d0-2.d0*sqrt(70.d0)/7.d0)/3.d0
         sps(3) = 0.d0
         sps(2) = -sps(4)
         sps(1) = -sps(5)
         wgt(5) = (322.d0-13.d0*sqrt(70.d0))/900.d0          ! Gaussian weights
         wgt(4) = (322.d0+13.d0*sqrt(70.d0))/900.d0
         wgt(3) = 128.d0/225.d0
         wgt(2) = wgt(4)
         wgt(1) = wgt(5)
       ELSE
         print *,'Highest Knod allowed is 5!'
         STOP
       ENDIF

       IF (Lnod .eq. 2) THEN
         gps(2) = 1.d0
         gps(1) = -gps(2)
       ELSEIF (Lnod .eq. 3) THEN
         gps(3) = 1.d0
         gps(2) = 0.d0
         gps(1) = -gps(3)
       ELSEIF (Lnod .eq. 4) THEN
         gps(4) = 1.d0
         gps(3) = 1.d0/3.d0
         gps(2) = -gps(3)
         gps(1) = -gps(4)
       ELSE
         print *,'Only up to cubic elements (Lnod = 3) is allowd!'
         STOP
       ENDIF

       ! Nothing to do here: Set up bases for Lagrangian (and derivative) inter/extrapolation
       !                     at solution and boundary nodes using these nodes as collocation points
       CALL GetSolnLagrangianBases

       ! Nothing to do here: Set up bases for Lagrangian (and derivative) inter/extrapolation
       !                     at geometry nodes using these nodes as collocation points
       CALL GetGeomLagrangianBases

       ! Nothing to do here: Set up values and derivatives of Right Radau function at solution points
       CALL GetRadauBases

       ! Nothing to do here: Set up nodes and element meshes
       IF (prob_type .eq. 10) THEN
         CALL SetupMesh_SU(mesh_filename)
       ELSE
         CALL SetupMesh(Nelx, Nely, Nelz, fac, dxrat, prob_type)
       ENDIF

       ! Nothing to do here: UNLESS we want to add new problem types
       !                     Set up IC, BC, and Source for a given problem
       ! NOTE: Make sure to add the corresponding exact solution for the problem type in Function u_xct
       CALL SetupICBCandSrc(Nel, NelB, NBndry, Knod, Lnod, prob_type)

       ! Nothing to do here: Set up volume and surface metrics
       CALL GetMetrics

       ! Nothing to do here: Now solve the convection+diffusion problem
       CALL SolveConvectionDiffusion(HuynhSolver_type, tIntegrator_type, Reyn, dt, numStep, dumpFreq, prob_type)
       ! Nothing to do here: Now solve the diffusion problem
!       IF (fast .eq. 0) THEN
!       CALL SolveDiffusion(Nel, Knod, Lnod, HuynhSolver_type, tIntegrator_type, dt, numStep, dumpFreq, prob_type)
!       ELSE
!       CALL SolveDiffusionFast(Nel, Knod, Lnod, HuynhSolver_type, tIntegrator_type, dt, numStep, dumpFreq, prob_type)
!       ENDIF

       ! Nothing to do here: Just deallocate arrays
       CALL deAllocate_Arrays

END PROGRAM ThreeD_Vorticity_Transport


SUBROUTINE Allocate_Arrays(N, NN, NP, NB, K, L)
       USE params
       USE variables

       implicit NONE
       integer N, NN, NP, NB, K, L
       integer Lm1, Ksq, KZ, L2
       ! number of faces (4 in 2D, 6 in 3D)

       Lm1 = L - 1
       Ksq = K*K
       KZ = K**(NDIM-2)
       L2 = L**(NDIM-1)
 
       allocate (xcoord(NN), ycoord(NN), zcoord(NN))
       allocate (nodeID(L**NDIM,N))
       allocate (elemID(FC+1,N))
       allocate (BndryNum(0:NP), BndryTag(NP), NoSlip(NP))
       allocate (BoundaryNodeID(L2,NB), BC_Values(K,K,NB), BC_VelNorm(K,K,NB), BC_VelParl(K,K,NB), BC_Psi(K,K,NB))
       allocate (BC_Switch(NB), BC_Switch_Lap(NB))
       allocate (boundaryPointsElemID(NB))
       allocate (VelocJump(K,K,NB))
       allocate (BndrySrc(K,KZ,K**NDIM,NB))
       allocate (sps(K), wgt(K))
       allocate (gps(L))
       allocate (Vort0(K,KZ,KZ,N), Vort(K,KZ,KZ,N), Uvel(K,K,KZ,N), Vvel(K,K,KZ,N), Wvel(KZ,KZ,KZ,N))
       allocate (Vol_Dx_iDxsi_j(NDIM,NDIM,K,K,KZ,N), Vol_Jac(K,K,KZ,N))
       allocate (Face_Acoef(K,KZ,0:FC,N), Face_Bcoef(K,KZ,0:FC,N), Face_Ccoef(K,KZ,0:FC,N))
       allocate (Face_Norm(K,KZ,0:FC,N), Face_Jac(K,KZ,0:FC,N))
       allocate (SolnBndryLgrangeBasis(K,0:1), SolnBndryGradLgrangeBasis(K,0:1))
       allocate (NodesRadau(K,0:1), NodesGradRadau(K,0:1))
       allocate (SolnNodesGradLgrangeBasis(K,K))
       allocate (GeomBndryLgrangeBasis(L,0:1), GeomBndryGradLgrangeBasis(L,0:1))
       allocate (GeomNodesLgrangeBasis(L,K), GeomNodesGradLgrangeBasis(L,K))

END SUBROUTINE Allocate_Arrays


SUBROUTINE deAllocate_Arrays
       USE variables

       deallocate (xcoord, ycoord, zcoord)
       deallocate (nodeID)
       deallocate (elemID)
       deallocate (BndryNum, BndryTag, NoSlip)
       deallocate (BoundaryNodeID, BC_Values, BC_VelNorm, BC_VelParl, BC_Psi)
       deallocate (BC_Switch, BC_Switch_Lap)
       deallocate (boundaryPointsElemID)
       deallocate (VelocJump)
       deallocate (BndrySrc)
       deallocate (sps, wgt)
       deallocate (gps)
       deallocate (Vort0, Vort, Uvel, Vvel)
       deallocate (Vol_Dx_iDxsi_j, Vol_Jac)
       deallocate (Face_Acoef, Face_Bcoef, Face_Norm, Face_Jac)
       deallocate (SolnBndryLgrangeBasis, SolnBndryGradLgrangeBasis)
       deallocate (NodesRadau, NodesGradRadau)
       deallocate (SolnNodesGradLgrangeBasis)
       deallocate (GeomBndryLgrangeBasis, GeomBndryGradLgrangeBasis)
       deallocate (GeomNodesLgrangeBasis, GeomNodesGradLgrangeBasis)

END SUBROUTINE deAllocate_Arrays


SUBROUTINE GetSolnLagrangianBases
       USE params
       USE variables

       implicit NONE

       integer i, j, k
       real*8 Denom, GradNumer(0:1), NodesGradNumer(1:Knod)

       ! Get Lagrange extrapolation basis for internal nodes, the Left (-1.0) and Right (1.0) boundaries
       !    as well as derivatives at the boudaries and internal nodes
       SolnBndryLgrangeBasis = 0.d0
       SolnBndryGradLgrangeBasis = 0.d0
       SolnNodesGradLgrangeBasis = 0.d0
       DO k = 1, Knod
         SolnBndryLgrangeBasis(k,0:1) = 1.d0
         Denom = 1.d0
         DO j = 1, Knod
           IF (j .eq. k) CYCLE
           Denom = Denom * (sps(k) - sps(j))       ! Basis denominator is common to all evaluations
           ! Get the numerators for the extrapolations to L+R
           SolnBndryLgrangeBasis(k,1) = SolnBndryLgrangeBasis(k,1) * ( 1.d0 - sps(j))
           SolnBndryLgrangeBasis(k,0) = SolnBndryLgrangeBasis(k,0) * (-1.d0 - sps(j))
           GradNumer = 1.d0
           NodesGradNumer = 1.d0
           DO i = 1, Knod
             IF (i .eq. k .or. i .eq. j) CYCLE
             ! Get the numerators for derivatives of extrpolations to L+R
             GradNumer(1) = GradNumer(1) * ( 1.d0 - sps(i))
             GradNumer(0) = GradNumer(0) * (-1.d0 - sps(i))
             ! Get the numerators for derivatives of interpolation to interior nodes
             NodesGradNumer(1:Knod) = NodesGradNumer(1:Knod) * (sps(1:Knod) - sps(i))
           ENDDO
           SolnBndryGradLgrangeBasis(k,0:1) = SolnBndryGradLgrangeBasis(k,0:1) + GradNumer(0:1)
           SolnNodesGradLgrangeBasis(k,1:Knod) = SolnNodesGradLgrangeBasis(k,1:Knod) + NodesGradNumer(1:Knod)
         ENDDO
         SolnBndryLgrangeBasis(k,0:1) = SolnBndryLgrangeBasis(k,0:1) / Denom
         SolnBndryGradLgrangeBasis(k,0:1) = SolnBndryGradLgrangeBasis(k,0:1) / Denom
         ! Get grads DOT_PROD[SolnNodesGradLgrangeBasis(k,j),u(k)] at a node/point j using data values u(k) at nodes k
         SolnNodesGradLgrangeBasis(k,1:Knod) = SolnNodesGradLgrangeBasis(k,1:Knod) / Denom
       ENDDO

END SUBROUTINE GetSolnLagrangianBases


SUBROUTINE GetGeomLagrangianBases
       USE params
       USE variables

       implicit NONE

       integer i, j, k
       real*8 Denom, GradNumer(0:1), NodesGradNumer(1:Knod)

       ! Get Lagrange extrapolation basis for internal nodes, the Left (-1.0) and Right (1.0) boundaries
       !    as well as derivatives at the boudaries and internal nodes
       GeomBndryLgrangeBasis = 0.d0
       GeomBndryGradLgrangeBasis = 0.d0
       GeomNodesLgrangeBasis = 0.d0
       GeomNodesGradLgrangeBasis = 0.d0
       DO k = 1, Lnod
         GeomBndryLgrangeBasis(k,0:1) = 1.d0
         GeomNodesLgrangeBasis(k,1:Knod) = 1.d0
         Denom = 1.d0
         DO j = 1, Lnod
           IF (j .eq. k) CYCLE
           Denom = Denom * (gps(k) - gps(j))       ! Basis denominator is common to all evaluations
           ! Get the numerators for the extrapolations to L+R
           GeomBndryLgrangeBasis(k,1) = GeomBndryLgrangeBasis(k,1) * ( 1.d0 - gps(j))
           GeomBndryLgrangeBasis(k,0) = GeomBndryLgrangeBasis(k,0) * (-1.d0 - gps(j))
           GeomNodesLgrangeBasis(k,1:Knod) = GeomNodesLgrangeBasis(k,1:Knod) * (sps(1:Knod) - gps(j))
           GradNumer = 1.d0
           NodesGradNumer = 1.d0
           DO i = 1, Lnod
             IF (i .eq. k .or. i .eq. j) CYCLE
             ! Get the numerators for derivatives of extrpolations to L+R
             GradNumer(1) = GradNumer(1) * ( 1.d0 - gps(i))
             GradNumer(0) = GradNumer(0) * (-1.d0 - gps(i))
             ! Get the numerators for derivatives of interpolation to interior nodes
             NodesGradNumer(1:Knod) = NodesGradNumer(1:Knod) * (sps(1:Knod) - gps(i))
           ENDDO
           GeomBndryGradLgrangeBasis(k,0:1) = GeomBndryGradLgrangeBasis(k,0:1) + GradNumer(0:1)
           GeomNodesGradLgrangeBasis(k,1:Knod) = GeomNodesGradLgrangeBasis(k,1:Knod) + NodesGradNumer(1:Knod)
         ENDDO
         GeomBndryLgrangeBasis(k,0:1) = GeomBndryLgrangeBasis(k,0:1) / Denom
         GeomBndryGradLgrangeBasis(k,0:1) = GeomBndryGradLgrangeBasis(k,0:1) / Denom
         GeomNodesLgrangeBasis(k,1:Knod) = GeomNodesLgrangeBasis(k,1:Knod) / Denom
         GeomNodesGradLgrangeBasis(k,1:Knod) = GeomNodesGradLgrangeBasis(k,1:Knod) / Denom
       ENDDO
 
END SUBROUTINE GetGeomLagrangianBases


SUBROUTINE GetRadauBases
       USE params
       USE variables

       implicit NONE

       integer n
       real*8 coef, coefd

       ! Get values and derivatives of right Radau on the interior nodes
       ! R_k(x) = (-1)^k/2 * (P(x)_k - P(x)_(k-1)), where P(x)_k is Legendre poly of order k
       ! Above Radau was defined by Huynh, Eqn. (3.15)

       ! The derivative at the interior can be shown to be:
       ! D[R_k(x),x] = (k/2) * (-1)^k * (P(x)_k + P(x)_(k-1)) / (x+1)
       NodesRadau = 0.d0
       NodesGradRadau = 0.d0
       IF (Knod .eq. 0) THEN
         Print *,'ERROR: Order 0 Radau is Undefined!'
         STOP
       ELSEIF (Knod .eq. 1) THEN
         NodesRadau(1,0) =  0.5d0
         NodesRadau(1,1) =  0.5d0

         NodesGradRadau(1,0) = -0.5d0
         NodesGradRadau(1,1) =  0.5d0
       ELSE
         coef  = 0.5d0 * (-1)**Knod
         coefd = 0.5d0 * (Knod * (-1)**Knod)
         DO n = 1, Knod
           ! This is a lazy and inefficient approach because we're evaluating the recursions twice for each function
           ! But who cares at this point because it's a one-time evaluation at the start

           ! R_k(x)|Right = R_k(-x)|Left
           NodesRadau(n,0) = coef * (LegendreP(Knod, sps(n)) - LegendreP(Knod-1, sps(n)))
           NodesRadau(n,1) = coef * (LegendreP(Knod,-sps(n)) - LegendreP(Knod-1,-sps(n)))

           ! D[R_k(x),x]|Right = -D[R_k(-x),x]|Left
           NodesGradRadau(n,0) = coefd * (LegendreP(Knod, sps(n)) + LegendreP(Knod-1, sps(n))) / (1.d0 + sps(n))
           NodesGradRadau(n,1) =-coefd * (LegendreP(Knod,-sps(n)) + LegendreP(Knod-1,-sps(n))) / (1.d0 - sps(n))
         ENDDO
       ENDIF

CONTAINS

  FUNCTION LegendreP(n,x)
       implicit NONE
       integer n
       real*8 x, LegendreP

       integer i
       real*8 xn, xnm1, xnm2

       IF (n .eq. 0) THEN
         LegendreP = 1.d0
       ELSEIF (n .eq. 1) THEN
         LegendreP = x
       ELSE
         xnm2 = 1.d0
         xnm1 = x
         DO i = 2, n
           xn = ((2*i-1)*x*xnm1 - (i-1)*xnm2) / i
           xnm2 = xnm1
           xnm1 = xn
         ENDDO
         LegendreP = xn
       ENDIF
       
  END FUNCTION LegendreP

END SUBROUTINE GetRadauBases


SUBROUTINE GetMetrics
       USE params
       USE variables
       USE CnvrtTensor2FemIndices
       USE omp_lib

       implicit NONE

       integer el, jx, jy, jz, lx, ly, lz, i, ip, KZ
       real*8 dxsi, deta, dzeta
       real*8, dimension(NDIM, NDIM) :: Dx_iDxsi_j, Dxsi_iDx_j
       real*8, dimension(NDIM, Lnod, Lnod, Lnod) :: xyzloc

   INTERFACE

     PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

     END FUNCTION Cross_Product


     PURE FUNCTION NodalVolume(Dx_iDxsi_j)
       implicit NONE

       ! First index(i) represents dx, dy, dz; Second index(j) represents dxsi, deta, dzeta
       real*8, dimension(3,3), intent(in) :: Dx_iDxsi_j
       real*8 NodalVolume

     END FUNCTION NodalVolume

   END INTERFACE


       Vol_Dx_iDxsi_j(1:NDIM, 1:NDIM, 1:Knod, 1:Knod, 1:Knod, 1:Nel) = 0.d0
       Vol_Jac(1:Knod, 1:Knod, 1:Knod, 1:Nel) = 0.d0
       Face_Acoef(1:Knod, 1:Knod, 0:FC, 1:Nel) = 0.d0
       Face_Bcoef(1:Knod, 1:Knod, 0:FC, 1:Nel) = 0.d0
       Face_Ccoef(1:Knod, 1:Knod, 0:FC, 1:Nel) = 0.d0
       Face_Norm(1:Knod, 1:Knod, 0:FC, 1:Nel) = 0.d0
       Face_Jac(1:Knod, 1:Knod, 0:FC, 1:Nel) = 0.d0

       DO el = 1, Nel

         DO lz = 1, Lnod
           DO ly = 1, Lnod
             DO lx = 1, Lnod
               xyzloc(1,lx,ly,lz) = xcoord(nodeID(t2f(lx,ly,lz),el))
               xyzloc(2,lx,ly,lz) = ycoord(nodeID(t2f(lx,ly,lz),el))
               xyzloc(3,lx,ly,lz) = zcoord(nodeID(t2f(lx,ly,lz),el))
             ENDDO
           ENDDO
         ENDDO

         ! Geometric metric stuff at all solution nodes (per element)
         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod

                             ! grad in x-dir                    no grad in y- and z- dir
                     dxsi  = GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * &
                                                                GeomNodesLgrangeBasis(lz,jz)
                             ! grad in y-dir                    no grad in z- and x- dir
                     deta  = GeomNodesGradLgrangeBasis(ly,jy) * GeomNodesLgrangeBasis(lz,jz) * &
                                                                GeomNodesLgrangeBasis(lx,jx)
                             ! grad in z-dir                    no grad in x- and y- dir
                     dzeta = GeomNodesGradLgrangeBasis(lz,jz) * GeomNodesLgrangeBasis(lx,jx) * &
                                                                GeomNodesLgrangeBasis(ly,jy)

                     ! Dx_Dxsi, Dy_Dxsi, Dz_Dxsi
                     Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el) = Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el) + dxsi  * xyzloc(:,lx,ly,lz)

                     ! Dx_Deta, Dy_Deta, Dz_Deta
                     Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el) = Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el) + deta  * xyzloc(:,lx,ly,lz)

                     ! Dx_Dzeta, Dy_Dzeta, Dz_Dzeta
                     Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el) = Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el) + dzeta * xyzloc(:,lx,ly,lz)

                   ENDDO
                 ENDDO
               ENDDO

               Vol_Jac(jx,jy,jz,el) = NodalVolume(Vol_Dx_iDxsi_j(1:3,1:3,jx,jy,jz,el))
             ENDDO
           ENDDO
         ENDDO

         ! Geometric metric stuff at all faces (per element)
         DO i = 0, 1
           ! Left (0) and Right (1) Face metric
           DO jz = 1, Knod
             DO jy = 1, Knod
               Dx_iDxsi_j = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod

                             ! grad in x-dir                    no grad in y- and z- dir
                     dxsi  = GeomBndryGradLgrangeBasis(lx,i)  * GeomNodesLgrangeBasis(ly,jy) * &
                                                                GeomNodesLgrangeBasis(lz,jz)
                             ! grad in y-dir                    no grad in z- and x- dir
                     deta  = GeomNodesGradLgrangeBasis(ly,jy) * GeomNodesLgrangeBasis(lz,jz) * &
                                                                GeomBndryLgrangeBasis(lx,i)
                             ! grad in z-dir                    no grad in x- and y- dir
                     dzeta = GeomNodesGradLgrangeBasis(lz,jz) * GeomBndryLgrangeBasis(lx,i)  * &
                                                                GeomNodesLgrangeBasis(ly,jy)

                     ! dx_dxsi, dy_dxsi, dz_dxsi
                     Dx_iDxsi_j(:,1) = Dx_iDxsi_j(:,1) + dxsi  * xyzloc(:,lx,ly,lz)

                     ! dx_deta, dy_deta, dz_deta
                     Dx_iDxsi_j(:,2) = Dx_iDxsi_j(:,2) + deta  * xyzloc(:,lx,ly,lz)

                     ! dx_dzeta, dy_dzeta, dz_dzeta
                     Dx_iDxsi_j(:,3) = Dx_iDxsi_j(:,3) + dzeta * xyzloc(:,lx,ly,lz)

                   ENDDO
                 ENDDO
               ENDDO
               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Dx_iDxsi_j(:,2),Dx_iDxsi_j(:,3))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Dx_iDxsi_j(:,3),Dx_iDxsi_j(:,1))
               ! Jac * dzta / dx_j
               Dxsi_iDx_j(:,3) = Cross_Product(Dx_iDxsi_j(:,1),Dx_iDxsi_j(:,2))

               Face_Jac(jy,jz,i,el) = NodalVolume(Dx_iDxsi_j)
               Face_Acoef(jy,jz,i,el) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,1))
               Face_Bcoef(jy,jz,i,el) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,2))
               Face_Ccoef(jy,jz,i,el) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,3))
               Face_Norm(jy,jz,i,el)  = Sqrt(Face_Acoef(jy,jz,i,el))
             ENDDO
           ENDDO
         ENDDO

         DO i = 0, 1
           ip = i + 2
           ! South (2) and North (3) Face metric
           DO jx = 1, Knod
             DO jz = 1, Knod
               Dx_iDxsi_j = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod

                             ! grad in y-dir                    no grad in z- and x- dir
                     deta  = GeomBndryGradLgrangeBasis(ly,i)  * GeomNodesLgrangeBasis(lz,jz) * &
                                                                GeomNodesLgrangeBasis(lx,jx)
                             ! grad in z-dir                    no grad in x- and y- dir
                     dzeta = GeomNodesGradLgrangeBasis(lz,jz) * GeomNodesLgrangeBasis(lx,jx) * &
                                                                GeomBndryLgrangeBasis(ly,i)
                             ! grad in x-dir                    no grad in y- and z- dir
                     dxsi  = GeomNodesGradLgrangeBasis(lx,jx) * GeomBndryLgrangeBasis(ly,i)  * &
                                                                GeomNodesLgrangeBasis(lz,jz)

                     ! dx_dxsi, dy_dxsi, dz_dxsi
                     Dx_iDxsi_j(:,1) = Dx_iDxsi_j(:,1) + dxsi  * xyzloc(:,lx,ly,lz)

                     ! dx_deta, dy_deta, dz_deta
                     Dx_iDxsi_j(:,2) = Dx_iDxsi_j(:,2) + deta  * xyzloc(:,lx,ly,lz)

                     ! dx_dzeta, dy_dzeta, dz_dzeta
                     Dx_iDxsi_j(:,3) = Dx_iDxsi_j(:,3) + dzeta * xyzloc(:,lx,ly,lz)

                   ENDDO
                 ENDDO
               ENDDO
               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Dx_iDxsi_j(:,2),Dx_iDxsi_j(:,3))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Dx_iDxsi_j(:,3),Dx_iDxsi_j(:,1))
               ! Jac * dzta / dx_j
               Dxsi_iDx_j(:,3) = Cross_Product(Dx_iDxsi_j(:,1),Dx_iDxsi_j(:,2))

               ! We don't use cyclic indexing; always use the smaller index first (allows for easier implementation)
               Face_Jac(jx,jz,ip,el)   = NodalVolume(Dx_iDxsi_j)
               Face_Acoef(jx,jz,ip,el) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,2))
               Face_Bcoef(jx,jz,ip,el) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,3))
               Face_Ccoef(jx,jz,ip,el) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,1))
               Face_Norm(jx,jz,ip,el)  = Sqrt(Face_Acoef(jx,jz,ip,el))
             ENDDO
           ENDDO
         ENDDO

         DO i = 0, 1
           ip = i + 4
           ! Bottom (4) and Top (5) Face metric
           DO jy = 1, Knod
             DO jx = 1, Knod
               Dx_iDxsi_j = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod

                             ! grad in z-dir                    no grad in x- and y- dir
                     dzeta = GeomBndryGradLgrangeBasis(lz,i)  * GeomNodesLgrangeBasis(lx,jx) * &
                                                                GeomNodesLgrangeBasis(ly,jy)
                             ! grad in x-dir                    no grad in y- and z- dir
                     dxsi  = GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * &
                                                                GeomBndryLgrangeBasis(lz,i)

                             ! grad in y-dir                    no grad in z- and x- dir
                     deta  = GeomNodesGradLgrangeBasis(ly,jy) * GeomBndryLgrangeBasis(lz,i)  * &
                                                                GeomNodesLgrangeBasis(lx,jx)

                     ! dx_dxsi, dy_dxsi, dz_dxsi
                     Dx_iDxsi_j(:,1) = Dx_iDxsi_j(:,1) + dxsi  * xyzloc(:,lx,ly,lz)

                     ! dx_deta, dy_deta, dz_deta
                     Dx_iDxsi_j(:,2) = Dx_iDxsi_j(:,2) + deta  * xyzloc(:,lx,ly,lz)

                     ! dx_dzeta, dy_dzeta, dz_dzeta
                     Dx_iDxsi_j(:,3) = Dx_iDxsi_j(:,3) + dzeta * xyzloc(:,lx,ly,lz)

                   ENDDO
                 ENDDO
               ENDDO
               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Dx_iDxsi_j(:,2),Dx_iDxsi_j(:,3))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Dx_iDxsi_j(:,3),Dx_iDxsi_j(:,1))
               ! Jac * dzta / dx_j
               Dxsi_iDx_j(:,3) = Cross_Product(Dx_iDxsi_j(:,1),Dx_iDxsi_j(:,2))

               Face_Jac(jx,jy,ip,el)   = NodalVolume(Dx_iDxsi_j)
               Face_Acoef(jx,jy,ip,el) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,3))
               Face_Bcoef(jx,jy,ip,el) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,1))
               Face_Ccoef(jx,jy,ip,el) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,2))
               Face_Norm(jx,jy,ip,el)  = Sqrt(Face_Acoef(jx,jy,ip,el))
             ENDDO
           ENDDO
         ENDDO

       ENDDO

END SUBROUTINE GetMetrics


PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

       Cross_Product(1) = x(2) * y(3) - y(2) * x(3)
       Cross_Product(2) = x(3) * y(1) - y(3) * x(1)
       Cross_Product(3) = x(1) * y(2) - y(1) * x(2)

END FUNCTION Cross_Product


PURE FUNCTION NodalVolume(Dx_iDxsi_j)
       implicit NONE

       ! First index(i) represents dx, dy, dz; Second index(j) represents dxsi, deta, dzeta
       real*8, dimension(3,3), intent(in) :: Dx_iDxsi_j
       real*8 NodalVolume

       real*8 cofactor(3)

   INTERFACE

     PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

     END FUNCTION Cross_Product

   END INTERFACE

       cofactor(1:3) = Cross_Product(Dx_iDxsi_j(1:3,1), Dx_iDxsi_j(1:3,2))

       NodalVolume = Dot_Product(cofactor(1:3), Dx_iDxsi_j(1:3,3))

END FUNCTION NodalVolume


SUBROUTINE SetupMesh_SU(mesh_filename)
       USE params
       USE variables

       implicit NONE

       integer L1, L2, FC1, idum, el, i, j, k, nb, ip1, jm1, eln, ndum
       character(20) cdum
       character(132) mesh_filename
       integer, allocatable, dimension(:,:) :: FaceNodes

       L1 = Lnod**(NDIM-1)
       L2 = Lnod**NDIM
       FC1 = 2*NDIM

       ! FaceNodes(i,j) points to the ith node of face j in both 2D and 3D (in CCW direction as per standard FEM)
       !                This is created to simplify searching for neighboring elements and assigning elemID
       Allocate (FaceNodes(2*(NDIM-1),2*NDIM))
       DO i = 1, 4
         ip1 = mod(i,4) + 1
         FaceNodes(1,i) = i
         FaceNodes(2,i) = ip1
         IF (NDIM .eq. 3) THEN
           FaceNodes(3,i) = ip1 + 4
           FaceNodes(4,i) = i + 4
         ENDIF
       ENDDO
       IF (NDIM .eq. 3) THEN
         DO i = 1, 4
           FaceNodes(i,5) = i
           FaceNodes(i,6) = i + 4 
         ENDDO
       ENDIF


       OPEN (10, File = trim(mesh_filename), Status = 'unknown')

       ! Dimension Number
       read(10,'(a6,i10)') cdum, ndum

       IF (ndum .ne. NDIM) THEN
         print *,'Problem Dimension = ',ndum,' NDIM = ',NDIM
         STOP
       ENDIF

       ! Nelem
       read(10,'(a6,i10)') cdum, Nel
       print *,cdum, Nel
       ! Connectivity data; Node numbers for each mesh (standard FEM notation?)
       DO el = 1, Nel
         read(10,*) idum, nodeID(1:L2,el), idum
       ENDDO

       ! nodeID is 0-based in SU2; convert to 1-based (this is NOT required for the circle/ellipse case below)
       nodeID = nodeID + 1

       ! Number of nodes plus one!!!!
       ! For now, this is a hack because we know there's an extra point inside the circle/ellipse
       ! Used by Gmsh for geometry generation. So, we have to take that out. But, this may not even
       ! be important because we mainly care about mesh nodes pointing to the correct node
       
       ! NOTE: For now, we're assuming mesh numbers are consecutive and always start from one
       !       We may want to generalize this later (IF necessary)
       read(10,'(a6,i10)') cdum, Nodes

       ! *********
       ! For the Sudden Expansion problem we don't have an extra point to take out
!       Nodes = Nodes - 1
       print *,cdum, Nodes
       ! *********
!       read(10,'(a)') cdum
       ! *********

       ! (x,y) coords of the mesh collocation points
       DO i = 1, Nodes
         read(10,*) xcoord(i), ycoord(i), zcoord(i)
       ENDDO

       ! Number of (solid and fluid) boundaries
       read(10,'(a6,i10)') cdum, NBndry
       print *,cdum, NBndry
       nb = 0
       BndryNum(0) = 0
       DO i = 1, NBndry
        ! Name of boundaries
         read(10,'(a12,a)') cdum, BndryTag(i)
         print *,cdum, BndryTag(i)
        ! Number of elements for each boundary
         read(10,'(a14,i10)') cdum, BndryNum(i)
         print *,cdum, BndryNum(i)
         DO j = 1, BndryNum(i)
           nb = nb + 1
           ! Node number for each boundary; this is the same as the numbers for the mesh that contains this boundary element
           read(10,*) idum, BoundaryNodeID(1:L1,nb)
         ENDDO
         ! We change the definition here, so BndryNum gives the cumulative value at each boundary/patch
         ! Therefore, the number of elements for each boundary i is BndryNum(i) - BndryNum(i-1)
         BndryNum(i) = BndryNum(i) + BndryNum(i-1)
       ENDDO

       ! BoundaryNodeID is 0-based in SU2; convert to 1-based (this is NOT required for the circle/ellipse case above)
       BoundaryNodeID = BoundaryNodeID + 1

       CLOSE (10)

       ! Poor man's approach to finding neighboring elements in a structured grid environment
       ! It also links the mesh faces to the boundary elements
       ! Seems to work!!

       elemID = 0
       boundaryPointsElemID = 0
       DO el = 1, Nel
         ! Source Element Face No.
         DO i = 1, FC1

           IF (elemID(i,el) .eq. 0) THEN
  Inner1:    DO eln = 1, Nel
               IF (eln .eq. el) cycle
               ! Target (neighbor) Element Face No.
               IF (NDIM .eq. 2) THEN
                 DO j = 1, FC1
                   IF (nodeID(FaceNodes(1,i),el) .eq. nodeID(FaceNodes(2,j),eln) .and.  &
                       nodeID(FaceNodes(2,i),el) .eq. nodeID(FaceNodes(1,j),eln) ) THEN
                     elemID(i,el) = eln
                     elemID(j,eln) = el
                     exit Inner1
                   ENDIF
                 ENDDO
               ENDIF
               IF (NDIM .eq. 3) THEN
                 DO j = 1, FC1
                   IF (nodeID(FaceNodes(1,i),el) .eq. nodeID(FaceNodes(2,j),eln) .and.  &
                       nodeID(FaceNodes(2,i),el) .eq. nodeID(FaceNodes(1,j),eln) .and.  &
                       nodeID(FaceNodes(3,i),el) .eq. nodeID(FaceNodes(4,j),eln) .and.  &
                       nodeID(FaceNodes(4,i),el) .eq. nodeID(FaceNodes(3,j),eln) ) THEN
                     elemID(i,el) = eln
                     elemID(j,eln) = el
                     exit Inner1
                   ENDIF
                 ENDDO
               ENDIF
             ENDDO Inner1
           ENDIF

           IF (elemID(i,el) .eq. 0) THEN
  Inner2:    DO eln = 1, NelB
               ! Target (neighbor) Element Face No.
               IF (NDIM .eq. 2) THEN
                 IF (nodeID(FaceNodes(1,i),el) .eq. BoundaryNodeID(2,eln) .and.  &
                     nodeID(FaceNodes(2,i),el) .eq. BoundaryNodeID(1,eln) ) THEN
                   elemID(i,el) = -eln
                   boundaryPointsElemID(eln) = el
                   exit Inner2
                 ENDIF
               ENDIF
               IF (NDIM .eq. 3) THEN
                 IF (nodeID(FaceNodes(1,i),el) .eq. BoundaryNodeID(1,eln) .and.  &
                     nodeID(FaceNodes(2,i),el) .eq. BoundaryNodeID(2,eln) .and.  &
                     nodeID(FaceNodes(3,i),el) .eq. BoundaryNodeID(3,eln) .and.  &
                     nodeID(FaceNodes(4,i),el) .eq. BoundaryNodeID(4,eln) ) THEN
                   elemID(i,el) = -eln
! IF (i.eq.5 .or. i.eq.6) elemID(i,el) = el
                   boundaryPointsElemID(eln) = el
                   exit Inner2
                 ENDIF
               ENDIF
             ENDDO Inner2
           ENDIF

         ENDDO
       ENDDO
if (.false.) then
print *,'nel ',nel
do i = 1,nel
print *,i,elemID(:,i)
enddo
print *,' '
print  *,'nodes ',nodes
do i = 1,nel
print *,i,nodeID(:,i)
enddo
print *,' '
print *,'nelb ',nelb
do i = 1,nelb
print *,i,boundaryPointsElemID(i)
enddo
print *,'nbndry ',nbndry
do i = 1,nbndry
print *,i,trim(BndryTag(i)),BndryNum(i-1)+1,BndryNum(i)
print *,'boundarynodeID '
do eln=BndryNum(i-1)+1,BndryNum(i)
print *,eln,BoundaryNodeID(:,eln)
enddo
enddo
print *,' '
do i = 1, nodes
print *,i, xcoord(i), ycoord(i), zcoord(i)
enddo
endif

       ! Paraview Geometry can handle only up to Lnod = 3
       ! Quick Paraview viz of geometry (useful when checking effect of fac randomization of grids)
       OPEN(unit = 8, file = 'Geometry.vtk', status = 'unknown')
       write(8,'(a)') '# vtk DataFile Version 3.0'
       write(8,'(a)') '2D Unstructured Grid of Quads'
       write(8,'(a)') 'ASCII'
       write(8,'(a)') ' '
       write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'

       write(8,'(a,i8,a)') 'POINTS ', Nodes, ' float'
       DO i = 1, Nodes
         write(8,*) real(xcoord(i)),real(ycoord(i)),real(zcoord(i))
       ENDDO

       write(8,'(a)') ' '
       write(8,'(a,i8,1x,i8)') 'CELLS ', Nel, (L2+1)*Nel
       DO i = 1, Nel
         write(8,*) L2, nodeID(1:L2,i)-1
       ENDDO

       write(8,'(a)') ' '
       write(8,'(a,i8)') 'CELL_TYPES ',Nel
       DO i = 1, Nel
         ! This works ONLY for Quad and Hex elements
         write(8,*) 3 * (NDIM + 1)
       ENDDO
       CLOSE(unit = 8)

       deAllocate (FaceNodes)

END SUBROUTINE SetupMesh_SU


SUBROUTINE SetupMesh(Nelx, Nely, Nelz, fac, dxrat, prob_type)
       USE params
       USE variables

       implicit NONE
       integer Nelx, Nely, Nelz, prob_type
       real*8 fac, dxrat

       ! A few grid types
       ! square box w/ random perturbation = 1
       ! square box w/ sinusoidal perturbation = 2
       ! a pair of concentric circles = 3
       ! square boxes 1 and 2 for periodic BC
       integer grid_type

       integer i, j, k, Lm1, nby, nl, nr, ns, nn, nb, nt, nc, nd, nd_sav, nadd
       real*8 xL, yL, zL, dx, dy, dz, rnd_num((Lnod**NDIM)*NDIM*Nel)
       real*8 R_in, R_out, theta_in, theta_out, dr, dtheta, R_j, theta_i, percent, tmp
       real*8 phi_in, phi_out, dphi, phi_k
       logical periodic
       logical left_boundary, right_boundary, south_boundary, north_boundary, bot_boundary, top_boundary

       periodic = .FALSE.
       CALL RANDOM_NUMBER(rnd_num)

       xL = 2.249d0*2*pi
       yL = 3.d0
       xL = 2.25d0
       yL = 2*2.25d0

       xL = 1.d0
       yL = 1.d0
       zL = 1.d0

       IF (prob_type .eq. 1) THEN
         xL = 2.d-1
         yL = 1.d0
         zL = 1.d1
       ENDIF

       grid_type = 2
       IF (Abs(dxrat - 1.d0) .gt. 1.d-4) grid_type = 4
       IF (prob_type .le. 2) periodic = .TRUE.
       IF (prob_type .gt. 6) grid_type = 3

       periodic = .TRUE.

       ! Quick Paraview viz of geometry (useful when checking effect of fac randomization of grids)
       OPEN(unit = 8, file = 'Geometry.vtk', status = 'unknown')
       write(8,'(a)') '# vtk DataFile Version 3.0'
       write(8,'(a)') '3D Unstructured Grid of Hexas'
       write(8,'(a)') 'ASCII'
       write(8,'(a)') ' '
       write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'

       ! poor man's attempt at creating unstructured mesh using simple structured data; just for testing
       ! Set geometry nodes (not the same as solution nodes)
       IF (grid_type .eq. 1 .or. grid_type .eq. 2) THEN

         IF (grid_type .eq. 1 .and. Lnod .gt. 2) THEN
           print *,'Only bilinear elements allowed for grid_type 1! '
           STOP
         ENDIF

         Lm1 = Lnod - 1
         dx = xL / (Nelx * Lm1)
         dy = yL / (Nely * Lm1)
         dz = zL / (Nelz * Lm1)

         nd = 0
         DO k = 1, Lm1*Nelz
           DO j = 1, Lm1*Nely
             DO i = 1, Lm1*Nelx
               nd = nd + 1
               xcoord(nd) = (i-1)*dx
               ycoord(nd) = (j-1)*dy
               zcoord(nd) = (k-1)*dz
               IF (fac .gt. 0.d0) THEN
                 IF (grid_type .eq. 1) THEN
                   xcoord(nd) = xcoord(nd) + fac*(0.5-rnd_num(nd))*dx
                   ycoord(nd) = ycoord(nd) + fac*(0.5-rnd_num(Nel+nd))*dy
                   zcoord(nd) = zcoord(nd) + fac*(0.5-rnd_num(2*Nel+nd))*dz
                 ELSE
                   tmp = sin((i-1)*2*pi*dx) * sin((j-1)*2*pi*dy) * sin((k-1)*2*pi*dz) * fac
!                   tmp = sin((i-1)*2*pi*dx) * sin((j-1)*2*pi*dy) * fac
                   xcoord(nd) = xcoord(nd) + tmp
                   ycoord(nd) = ycoord(nd) + tmp
                   zcoord(nd) = zcoord(nd) + tmp
                 ENDIF
               ENDIF
             ENDDO
             nd = nd + 1
             xcoord(nd) = xL
             ycoord(nd) = (j-1)*dy
             zcoord(nd) = (k-1)*dz
           ENDDO

           DO i = 1, Lm1*Nelx
             nd = nd + 1
             xcoord(nd) = (i-1)*dx
             ycoord(nd) = yL
             zcoord(nd) = (k-1)*dz
           ENDDO
           nd = nd + 1
           xcoord(nd) = xL
           ycoord(nd) = yL
           zcoord(nd) = (k-1)*dz
         ENDDO

         DO j = 1, Lm1*Nely
           DO i = 1, Lm1*Nelx
             nd = nd + 1
!                   tmp = sin((i-1)*2*pi*dx) * sin((j-1)*2*pi*dy) * fac
             xcoord(nd) = (i-1)*dx !+ tmp
             ycoord(nd) = (j-1)*dy !+ tmp
             zcoord(nd) = zL
           ENDDO
           nd = nd + 1
           xcoord(nd) = xL
           ycoord(nd) = (j-1)*dy
           zcoord(nd) = zL
         ENDDO

         DO i = 1, Lm1*Nelx
           nd = nd + 1
           xcoord(nd) = (i-1)*dx
           ycoord(nd) = yL
           zcoord(nd) = zL
         ENDDO
         nd = nd + 1
         xcoord(nd) = xL
         ycoord(nd) = yL
         zcoord(nd) = zL

         ! This spherical coordinate-based unit was used to check the correctness of the metrics computation
         IF (.false.) THEN
           R_in = 0.5d0
           R_out = 1.d0
           theta_in = pi/6 !+ percent*2.d0*pi
!         theta_in = 0
           theta_out = pi/3 !*(1.0d0 + percent)
           phi_in = pi/4
!         phi_in = 0
           phi_out = pi/4+pi/6
           dr = (R_out - R_in) / (Nelx * Lm1)
           dtheta = (theta_out - theta_in) / (Nely * Lm1)
           dphi = (phi_out - phi_in) / (Nelz * Lm1)

           nd = 0
           DO k = 1, Lm1*Nelz+1
             DO j = 1, Lm1*Nely+1
               DO i = 1, Lm1*Nelx+1
                 nd = nd + 1
                 theta_i = theta_in + (j-1)*dtheta
                 R_j = R_in + (i-1)*dr
                 phi_k = phi_in + (k-1)*dphi
                 if (i .eq. Lm1*Nelx+1) R_j = R_out
                 xcoord(nd) = R_j * cos(theta_i) * cos(phi_k)
                 ycoord(nd) = R_j * sin(theta_i) * cos(phi_k)
                 zcoord(nd) = R_j * sin(phi_k)
               ENDDO
             ENDDO
           ENDDO
         ENDIF

       ELSEIF (grid_type .eq. 3) THEN

         Lm1 = Lnod - 1
         R_in = 0.5d0
         R_out = 1.d0
!         zL = 1.d0

         dr = (R_out - R_in) / (Nely * Lm1)
         dz = zL / (Nelz * Lm1)
!         print *,'enter percent rotation '
!         read *, percent
         theta_in = 0.d0 !+ percent*2.d0*pi
         theta_out = 2.d0*pi !*(1.0d0 + percent)
!         print *,theta_out - theta_in,2.d0*pi
         dtheta = (theta_out - theta_in) / (Nelx * Lm1)

         nd = 0
         DO k = 1, Lm1*Nelz+1
           DO j = 1, Lm1*Nely+1
             DO i = 1, Lm1*Nelx
               nd = nd + 1
               theta_i = theta_in + (i-1)*dtheta
               R_j = R_in + (j-1)*dr
               if (j .eq. Lm1*Nely+1) R_j = R_out
               xcoord(nd) = R_j * cos(-theta_i)
               ycoord(nd) = R_j * sin(-theta_i)
               zcoord(nd) = (k-1)*dz
               if (k .eq. Lm1*Nelz+1) zcoord(nd) = zL
             ENDDO
           ENDDO
         ENDDO

       ELSEIF (grid_type .eq. 4) THEN

         Lm1 = Lnod - 1
         ! only for Lm1 = 1 and even Nelx/Nely
         IF (Lm1 .ne. 1) THEN
           print *,'only linear elements allowed'
           stop
         ENDIF
         IF (Mod(Nelx,2) .eq. 1 .or. Mod(Nely,2) .eq. 1 .or. Mod(Nelz,2) .eq. 1) THEN
           print *,'only even element counts allowed'
           stop
         ENDIF

         dx = 0.5d0 * xL * (1.d0 - dxrat) / (1.d0 - dxrat**(Nelx/2))
         dy = 0.5d0 * yL * (1.d0 - dxrat) / (1.d0 - dxrat**(Nely/2))
         dz = 0.5d0 * zL * (1.d0 - dxrat) / (1.d0 - dxrat**(Nelz/2))
         print *,'Smallest and largest dx ', dx, dx * (1.d0 - dxrat**(Nelx/2-1)) / (1.d0 - dxrat)
         nd = 0
         DO k = 1, Nely + 1
           DO j = 1, Nely + 1
             DO i = 1, Nelx + 1
               nd = nd + 1
               IF (i .le. Nelx/2) THEN
                 xcoord(nd) = dx * (1.d0 - dxrat**(i-1)) / (1.d0 - dxrat)
               ELSEIF (i .eq. Nelx/2 + 1) THEN
                 xcoord(nd) = 0.5d0 * xL
               ELSEIF (i .le. Nelx) THEN
                 xcoord(nd) = xL - dx * (1.d0 - dxrat**(Nelx + 1 - i)) / (1.d0 - dxrat)
               ELSE
                 xcoord(nd) = xL
               ENDIF
               IF (j .le. Nely/2) THEN
                 ycoord(nd) = dy * (1.d0 - dxrat**(j-1)) / (1.d0 - dxrat)
               ELSEIF (j .eq. Nely/2 + 1) THEN
                 ycoord(nd) = 0.5d0 * yL
               ELSEIF (j .le. Nely) THEN
                 ycoord(nd) = yL - dy * (1.d0 - dxrat**(Nely + 1 - j)) / (1.d0 - dxrat)
               ELSE
                 ycoord(nd) = yL
               ENDIF
               IF (k .le. Nelz/2) THEN
                 zcoord(nd) = dz * (1.d0 - dxrat**(k-1)) / (1.d0 - dxrat)
               ELSEIF (k .eq. Nelz/2 + 1) THEN
                 zcoord(nd) = 0.5d0 * zL
               ELSEIF (k .le. Nelz) THEN
                 zcoord(nd) = zL - dz * (1.d0 - dxrat**(Nelz + 1 - k)) / (1.d0 - dxrat)
               ELSE
                 zcoord(nd) = zL
               ENDIF
             ENDDO
           ENDDO
         ENDDO

       ELSE
         print *,'selected grid_type is not supported! '
         STOP
       ENDIF

       ! Paraview viz
       write(8,'(a,i8,a)') 'POINTS ', nd, ' float'
       DO i = 1, nd
         write(8,*) real(xcoord(i)),real(ycoord(i)),real(zcoord(i))
       ENDDO
       nd_sav = nd

       ! Set mesh connectivity to nodes (nodeID); also get connectivity to neighboring meshes (elemID)
       ! nodeID starts from SW corners and goes CCW; elemID starts with S face and goes CCW
       ! NOTE: This is of course a poor man's (structured grid) approach to the general problem
       !       Also, the -ve numbers refer to boundary faces, which are assigned in the BC_IC setup routine
       nby = 0
       nc = 0
       nd = 0
       DO k = 1, Nelz
         DO j = 1, Nely
           DO i = 1, Nelx
             nc = nc + 1
             nd = nd + 1

             left_boundary = .false.
             right_boundary = .false.
             south_boundary = .false.
             north_boundary = .false.
             bot_boundary = .false.
             top_boundary = .false.

             IF (k .eq. 1) THEN
               IF (periodic) THEN
                 elemID(5,nc) = nc + Nelx*Nely*(Nelz-1)
               ELSE
               bot_boundary = .true.
               nby = nby + 1
               nb = nby
               elemID(5,nc) = -nby
               boundaryPointsElemID(nby) = nc
               ENDIF
             ELSE
               elemID(5,nc) = nc - Nelx*Nely
             ENDIF
             IF (k .eq. Nelz) THEN
               IF (periodic) THEN
                 elemID(6,nc) = nc - Nelx*Nely*(Nelz-1)
               ELSE
               top_boundary = .true.
               nby = nby + 1
               nt = nby
               elemID(6,nc) = -nby
               boundaryPointsElemID(nby) = nc
             ENDIF
             ELSE
               elemID(6,nc) = nc + Nelx*Nely
             ENDIF

             IF (j .eq. 1) THEN
!             IF (periodic) THEN
!               elemID(1,nc) = nc + Nelx * (Nely - 1)
!             ELSE
               south_boundary = .true.
               nby = nby + 1
               ns = nby
               elemID(1,nc) = -nby
               boundaryPointsElemID(nby) = nc
!             ENDIF
             ELSE
               elemID(1,nc) = nc - Nelx
             ENDIF
             IF (j .eq. Nely) THEN
!             IF (periodic) THEN
!               elemID(3,nc) = nc - Nelx * (Nely - 1)
!             ELSE
               north_boundary = .true.
               nby = nby + 1
               nn = nby
               elemID(3,nc) = -nby
               boundaryPointsElemID(nby) = nc
!             ENDIF
             ELSE
               elemID(3,nc) = nc + Nelx
             ENDIF

             ! for concentric circles there is no boundary in the i direction
             ! also true of periodic BC
             IF (i .eq. 1) THEN
!               IF (grid_type .eq. 3 .or. periodic) THEN
             IF (grid_type .eq. 3) THEN
                 elemID(4,nc) = nc + Nelx - 1
               ELSE
                 left_boundary = .true.
                 nby = nby + 1
                 nl = nby
                 elemID(4,nc) = -nby
                 boundaryPointsElemID(nby) = nc
               ENDIF
             ELSE             
               elemID(4,nc) = nc - 1
             ENDIF
             IF (i .eq. Nelx) THEN
!               IF (grid_type .eq. 3 .or. periodic) THEN
             IF (grid_type .eq. 3) THEN
                 elemID(2,nc) = nc + 1 - Nelx
               ELSE
                 right_boundary = .true.
                 nby = nby + 1
                 nr = nby
                 elemID(2,nc) = -nby
                 boundaryPointsElemID(nby) = nc
               ENDIF
             ELSE
               elemID(2,nc) = nc + 1
             ENDIF

             nodeID(1,nc) = nd
             nodeID(2,nc) = nd + Lm1           
             nodeID(3,nc) = nd + Lm1 * (Lm1*Nelx + 2)
             nodeID(4,nc) = nd + Lm1 * (Lm1*Nelx + 1)

             IF (grid_type .eq. 3) THEN
               nodeID(3,nc) = nodeID(3,nc) - Lm1
               nodeID(4,nc) = nodeID(4,nc) - Lm1
               IF (i .eq. Nelx) THEN
                 nodeID(2,nc) = nodeID(2,nc) - Lm1*Nelx
                 nodeID(3,nc) = nodeID(3,nc) - Lm1*Nelx
               ENDIF
             ENDIF

             nadd = (Nelx + 1) * (Nely + 1)
             IF (grid_type .eq. 3) nadd = Nelx * (Nely + 1)
             nodeID(5,nc) = nodeID(1,nc) + nadd
             nodeID(6,nc) = nodeID(2,nc) + nadd
             nodeID(7,nc) = nodeID(3,nc) + nadd
             nodeID(8,nc) = nodeID(4,nc) + nadd

             IF (left_boundary) THEN
             ! CCW rotation for meshing
!             BoundaryNodeID(2,nl) = nodeID(1,nc)
!             BoundaryNodeID(1,nl) = nodeID(4,nc)
             ! For now, I'm assuming all sides are in positive x/y direction 
             ! NOTE: need to account for direction later
               BoundaryNodeID(2,nl) = nodeID(4,nc)
               BoundaryNodeID(1,nl) = nodeID(1,nc)
!  Should this be 3 and 4??
               BoundaryNodeID(4,nl) = nodeID(8,nc)
               BoundaryNodeID(3,nl) = nodeID(5,nc)
             ENDIF

             IF (right_boundary) THEN
               BoundaryNodeID(1,nr) = nodeID(2,nc)
               BoundaryNodeID(2,nr) = nodeID(3,nc)

               BoundaryNodeID(3,nr) = nodeID(6,nc)
               BoundaryNodeID(4,nr) = nodeID(7,nc)
             ENDIF

             IF (south_boundary) THEN
               BoundaryNodeID(1,ns) = nodeID(1,nc)
               BoundaryNodeID(2,ns) = nodeID(2,nc)

               BoundaryNodeID(3,ns) = nodeID(5,nc)
               BoundaryNodeID(4,ns) = nodeID(6,nc)
             ENDIF

             IF (north_boundary) THEN
             ! CCW rotation for meshing
!             BoundaryNodeID(1,nn) = nodeID(3,nc)
!             BoundaryNodeID(2,nn) = nodeID(4,nc)
             ! For now, I'm assuming all sides are in positive x/y direction 
             ! NOTE: need to account for direction later
               BoundaryNodeID(1,nn) = nodeID(4,nc)
               BoundaryNodeID(2,nn) = nodeID(3,nc)

               BoundaryNodeID(3,nn) = nodeID(8,nc)
               BoundaryNodeID(4,nn) = nodeID(7,nc)
             ENDIF

             IF (bot_boundary) THEN
               BoundaryNodeID(1,nb) = nodeID(1,nc)
               BoundaryNodeID(2,nb) = nodeID(2,nc)

               BoundaryNodeID(4,nb) = nodeID(3,nc)
               BoundaryNodeID(3,nb) = nodeID(4,nc)
             ENDIF

             IF (top_boundary) THEN
               BoundaryNodeID(1,nt) = nodeID(5,nc)
               BoundaryNodeID(2,nt) = nodeID(6,nc)

               BoundaryNodeID(4,nt) = nodeID(7,nc)
               BoundaryNodeID(3,nt) = nodeID(8,nc)
             ENDIF

           IF (.false.) THEN

             IF (Lm1 .gt. 1) THEN
               IF (Lm1 .eq. 2) THEN
               nd = nd + 1
               nodeID(5,nc) = nodeID(1,nc) + 1
               nodeID(7,nc) = nodeID(4,nc) + 1
               nodeID(8,nc) = nodeID(1,nc) + Lm1*Nelx + 1
               IF (grid_type .eq. 3) nodeID(8,nc) = nodeID(8,nc) - 1
               nodeID(6,nc) = nodeID(8,nc) + Lm1
               IF (i .eq. Nelx .and. grid_type .eq. 3) nodeID(6,nc) = nodeID(6,nc) - Lm1*Nelx
               nodeID(9,nc) = nodeID(8,nc) + 1

               IF (left_boundary)  BoundaryNodeID(3,nl) = nodeID(8,nc)
               IF (right_boundary) BoundaryNodeID(3,nr) = nodeID(6,nc)
               IF (south_boundary) BoundaryNodeID(3,ns) = nodeID(5,nc)
               IF (north_boundary) BoundaryNodeID(3,nn) = nodeID(7,nc)

               ELSEIF (Lm1 .eq. 3) THEN
               nd = nd + 2
               nodeID(5,nc) = nodeID(1,nc) + 1
               nodeID(6,nc) = nodeID(1,nc) + 2
               nodeID(9,nc) = nodeID(4,nc) + 2
               nodeID(10,nc) = nodeID(4,nc) + 1
               nodeID(12,nc) = nodeID(1,nc) + Lm1*Nelx + 1
               nodeID(11,nc) = nodeID(1,nc) + 2 * (Lm1*Nelx + 1)
               IF (grid_type .eq. 3) THEN
                 nodeID(12,nc) = nodeID(12,nc) - 1
                 nodeID(11,nc) = nodeID(11,nc) - 2
               ENDIF
               nodeID(7,nc) = nodeID(12,nc) + Lm1
               nodeID(8,nc) = nodeID(11,nc) + Lm1
               IF (i .eq. Nelx .and. grid_type .eq. 3) THEN
                 nodeID(7,nc) = nodeID(7,nc) - Lm1*Nelx
                 nodeID(8,nc) = nodeID(8,nc) - Lm1*Nelx
               ENDIF
               nodeID(13,nc) = nodeID(12,nc) + 1
               nodeID(14,nc) = nodeID(12,nc) + 2
               nodeID(15,nc) = nodeID(11,nc) + 2
               nodeID(16,nc) = nodeID(11,nc) + 1

               IF (left_boundary) THEN
                 BoundaryNodeID(3,nl) = nodeID(12,nc)
                 BoundaryNodeID(4,nl) = nodeID(11,nc)
               ENDIF
               IF (right_boundary) THEN
                 BoundaryNodeID(3,nr) = nodeID(7,nc)
                 BoundaryNodeID(4,nr) = nodeID(8,nc)
               ENDIF
               IF (south_boundary) THEN
                 BoundaryNodeID(3,ns) = nodeID(5,nc)
                 BoundaryNodeID(4,ns) = nodeID(6,nc)
               ENDIF
               IF (north_boundary) THEN
                 BoundaryNodeID(3,nn) = nodeID(10,nc)
                 BoundaryNodeID(4,nn) = nodeID(9,nc)
               ENDIF

               ELSE
               print *,'ERROR: Up to bicubic elements supported! Set Lnod = 3, 2 or 1 '
               STOP
               ENDIF
             ENDIF
           ENDIF
           ENDDO

           ! accounting for the node Nelx+1, except for grid_type 3
           IF (grid_type .ne. 3) nd = nd + 1
           ! doing same for higher order meshes
           IF (.false.) THEN
           IF (Lm1 .eq. 2) THEN
             nd = nd + Lm1*Nelx + 1
             IF (grid_type .eq. 3) nd = nd - 1
           ENDIF
           IF (Lm1 .eq. 3) THEN
             nd = nd + 2 * (Lm1*Nelx + 1)
             IF (grid_type .eq. 3) nd = nd - 2
           ENDIF
           ENDIF
         ENDDO

         nd = nd + Nelx
         IF (grid_type .ne. 3) nd = nd + 1
       ENDDO


!       DO i = 1, nc 
!        print *,i,nodeID(1:(Lm1+1)**2,i)
!        print *,i,nodeID(1:8,i)
!        print *,i,elemID(1:6,i)
!       ENDDO

       ! Need this for the periodic case
       NelB = nby

       IF (Lnod > 3) RETURN

       ! Paraview viz
       write(8,'(a)') ' '
       write(8,'(a,i8,1x,i8)') 'CELLS ', nc, (Lnod**NDIM+1)*nc
       DO i = 1, nc
         write(8,*) Lnod**NDIM, nodeID(1:Lnod**NDIM,i)-1
       ENDDO
       write(8,'(a)') ' '
       write(8,'(a,i8)') 'CELL_TYPES ',nc
       DO i = 1, nc
!         write(8,*) Lnod**3 + 1
         write(8,*) 12
       ENDDO
       CLOSE(unit = 8)

END SUBROUTINE SetupMesh


SUBROUTINE SetupICBCandSrc(N, NB, NP, K, L, prob_type)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer N, NB, NP, K, L, prob_type

       integer i, j, kx, ky, kz, lx, ly, lz, jx, jy, el
       real*8 x, y, z, xy, eps, ys, xm, ym, HalfHeight, y2, y1, tmp

       real*8 y_in_min, y_in_max, HalfHeight_in, y_out_min, y_out_max, HalfHeight_out

       real*8 theta2,theta1,thetam,thetad,rad2,rad1,radm,radd
       real*8, dimension(Lnod,Lnod) :: xloc, yloc

       eps = 1.d-4

!       OPEN(UNIT = 3, file = 'Solution_Nodes.csv', status = 'unknown')

       ! For now, we are assuming a constant value of BC across a single patch at the boundary
       ! That is, the geometry is topologically a square, each side of which has a constant
       ! Neumann or Dirichlet BC (or periodic). The patches are numbered from the bottom and go counter-clockwise
       ! S=1; E=2; N=3; W=4.  The BC switch for these patches is Dirichlet=1, or Neumann=2

       BC_Switch(1:NB) = DefaultBC
       BC_Values(1:K,1:K,1:NB) = 0.d0   ! This is legacy and should be removed later

       BC_VelNorm(1:K,1:K,1:NB) = 0.d0   ! Wall normal velocity
       BC_VelParl(1:K,1:K,1:NB) = 0.d0   ! Wall parallel velocity
       BC_Psi(1:K,1:K,1:NB) = 0.d0       ! Wall normal velocity (in terms of psi)

       NoSlip(1:NP) = .true.         ! Set NoSlip to all walls (as default)

       IF (prob_type .eq. 1) THEN

         ! Periodic in x, Neumann in y, moving wall in z; bottom velocity is 1; top is 0

!          BC_Switch(1:NB) = DirichletBC
          DO i = 1, NB
            IF (ABS(zcoord(BoundaryNodeID(1,i))) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(2,i))) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(3,i))) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(4,i))) < 1.d-6) THEN
                BC_Switch(i) = DirichletBC
                BC_VelParl(1:K,1:K,i) = 1.d0   ! Wall parallel velocity
            ENDIF
            IF ((ABS(ycoord(BoundaryNodeID(1,i))) < 1.d-6 .and. &
                 ABS(ycoord(BoundaryNodeID(2,i))) < 1.d-6 .and. &
                 ABS(ycoord(BoundaryNodeID(3,i))) < 1.d-6 .and. &
                 ABS(ycoord(BoundaryNodeID(4,i))) < 1.d-6) .or. &
                (ABS(ycoord(BoundaryNodeID(1,i))) < 1.d-6 .and. &
                 ABS(ycoord(BoundaryNodeID(2,i))-1) < 1.d-6 .and. &
                 ABS(ycoord(BoundaryNodeID(3,i))-1) < 1.d-6 .and. &
                 ABS(ycoord(BoundaryNodeID(4,i))-1) < 1.d-6)) THEN
                 BC_Switch(i) = NeumannBC
                 BC_VelParl(1:K,1:K,i) = 0.d0   ! Wall parallel velocity
            ENDIF
          ENDDO

         Vort0 = 0.d0

       ELSEIF (prob_type .eq. 3) THEN

          ! Cubic Cavity Problem; top velocity is 1

          BC_Switch(1:NB) = DirichletBC
          DO i = 1, NB
            IF (ABS(zcoord(BoundaryNodeID(1,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(2,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(3,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(4,i)) - 1.d0) < 1.d-6) THEN
              BC_VelParl(1:K,1:K,i) = 1.d0   ! Wall parallel velocity
            ENDIF
          ENDDO

       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    along z-dir
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO
!               vort0(kx,ky,kz,el) = exp(-(x**2+y**2+z**2))/pi
               vort0(kx,ky,kz,el) = 0
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       ELSEIF (prob_type .eq. 5) THEN

          ! Cubic Cavity Problem; top velocity is 1

          BC_Switch(1:NB) = NeumannBC

       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    along z-dir
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO
               Vort0(kx,ky,kz,el) = (1.5d0*Cos(pi*x) - 0.5d0*Cos(3.d0*pi*x)) * &
                                    (1.5d0*Cos(pi*y) - 0.5d0*Cos(3.d0*pi*y)) * &
                                    (1.5d0*Cos(pi*z) - 0.5d0*Cos(3.d0*pi*z)) ! IC
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       ELSEIF (prob_type .eq. 55) THEN

          ! Cubic Cavity Problem; top velocity is 1

          BC_Switch(1:NB) = DirichletBC
          DO i = 1, NB
            IF (ABS(zcoord(BoundaryNodeID(1,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(2,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(3,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(4,i)) - 1.d0) < 1.d-6) THEN
              BC_VelParl(1:K,1:K,i) = 0.d0   ! Wall parallel velocity
            ENDIF
          ENDDO

       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    along z-dir
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO
!               vort0(kx,ky,kz,el) = exp(-(x**2+y**2+z**2))/pi
               vort0(kx,ky,kz,el) = Sin(pi*x)*Sin(pi*y)*Sin(pi*z) * 3 * pi**2
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       ELSEIF (prob_type .eq. 6) THEN

          ! Cubic Cavity Problem; top velocity is 1

          BC_Switch(1:NB) = DirichletBC
          BC_VelParl(1:K,1:K,NB) = 0.d0   ! Wall parallel velocity
          BC_Values(1:K,1:K,NB) = 0.d0   ! Wall parallel velocity
if (.false.) then
          DO i = 1, NB
            IF (ABS(zcoord(BoundaryNodeID(1,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(2,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(3,i)) - 1.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(4,i)) - 1.d0) < 1.d-6 .or.  &
                ABS(zcoord(BoundaryNodeID(1,i)) - 0.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(2,i)) - 0.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(3,i)) - 0.d0) < 1.d-6 .and. &
                ABS(zcoord(BoundaryNodeID(4,i)) - 0.d0) < 1.d-6) THEN
!                  BC_Switch(i) = NeumannBC
                  BC_Values(1:K,1:K,i) = 0.d0   ! Wall parallel velocity
            ENDIF
          ENDDO
endif
          DO i = 1, NB
            IF (ABS(ycoord(BoundaryNodeID(1,i)) - 1.d0) < 1.d-6 .and. &
                ABS(ycoord(BoundaryNodeID(2,i)) - 1.d0) < 1.d-6 .and. &
                ABS(ycoord(BoundaryNodeID(3,i)) - 1.d0) < 1.d-6 .and. &
                ABS(ycoord(BoundaryNodeID(4,i)) - 1.d0) < 1.d-6) THEN
              BC_VelParl(1:K,1:K,i) = 1.d0   ! Wall parallel velocity
            ENDIF
          ENDDO

       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    along z-dir
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO
!               vort0(kx,ky,kz,el) = exp(-(x**2+y**2+z**2))/pi
!               vort0(kx,ky,kz,el) = Sin(pi*x)*Sin(pi*y) * 2 * pi**2
!               vort0(kx,ky,kz,el) = Sin(pi*x)*Sin(pi*y)*Cos(pi*z) * 3 * pi**2
             ENDDO
           ENDDO
         ENDDO
       ENDDO

 !!!!! THIS OPTION NEEDS TO BE CORRECTED FOR 3D
       ELSEIF (prob_type .eq. 10) THEN

         ! Need to make BC_Psi generic to account for width of inlet automatically instead of the current manual setup
         ! We want a uniform inlet velocity of 1, so Psi = y / (2*HalfHeight)
         ! But we want the values at the wall solution points, so y = 0.5 * ( (y_2 + y_1) + sps * (y_2 - y_1) )
         ! where y_2 and y_1 are the coords of the end points of a boundary element and sps is the solution point position in local
         ! coords
         BC_Psi = 0.d0
         Vort0 = 0.d0
         BC_Switch(1:NB) = DirichletBC
         NoSlip = .FALSE.

         ! This is hack
         y_in_max = -1000.
         y_in_min =  1000.
         y_out_max = -1000.
         y_out_min =  1000.
         DO i = 1, NBndry
           IF (trim(BndryTag(i)) .eq. 'Inlet') THEN
             DO j = BndryNum(i-1)+1,BndryNum(i)
               y_in_max = Max(y_in_max, MaxVal(ycoord(BoundaryNodeID(:,j))))
               y_in_min = Min(y_in_min, MinVal(ycoord(BoundaryNodeID(:,j))))
             ENDDO
           ENDIF
           IF (trim(BndryTag(i)) .eq. 'Outlet') THEN
             DO j = BndryNum(i-1)+1,BndryNum(i)
               y_out_max = Max(y_out_max, MaxVal(ycoord(BoundaryNodeID(:,j))))
               y_out_min = Min(y_out_min, MinVal(ycoord(BoundaryNodeID(:,j))))
             ENDDO
           ENDIF
         ENDDO

         HalfHeight_in = 0.5d0 * (y_in_max - y_in_min)
         HalfHeight_out = 0.5d0 * (y_out_max - y_out_min)

         DO i = 1, NBndry

           IF (trim(BndryTag(i)) .eq. 'Top') THEN
             ! Set slip BC to not allow diffusion
             NoSlip(i) = .true.
             DO j = BndryNum(i-1)+1, BndryNum(i)
               BC_Psi(1:K,1:K,j) = 0.5d0
             ENDDO
           ENDIF

           IF (trim(BndryTag(i)) .eq. 'Bottom_In' .or. &
               trim(BndryTag(i)) .eq. 'Bottom_Out' .or. &
               trim(BndryTag(i)) .eq. 'Step' .or. &
               trim(BndryTag(i)) .eq. 'Bottom' ) THEN
             ! Set slip BC to not allow diffusion
             NoSlip(i) = .true.
             DO j = BndryNum(i-1)+1, BndryNum(i)
               BC_Psi(1:K,1:K,j) = 0.d0
             ENDDO
           ENDIF

           IF (trim(BndryTag(i)) .eq. 'Periodic_Front' .or. &
               trim(BndryTag(i)) .eq. 'Periodic_Back' ) THEN
             ! Set slip BC to not allow diffusion
             NoSlip(i) = .false.
             DO j = BndryNum(i-1)+1, BndryNum(i)
               BC_Psi(1:K,1:K,j) = 0.d0
               BC_Switch(j) = NeumannBC
               BC_Psi(1:K,1:K,j) = 0.d0
               BC_Values(1:K,1:K,j) = 0.d0
             ENDDO
           ENDIF

           IF (trim(BndryTag(i)) .eq. 'Inlet') THEN
             ! Set slip BC to not allow diffusion
             NoSlip(i) = .false.
             DO j = BndryNum(i-1)+1, BndryNum(i)
               DO kz = 1, Knod
                 DO ky = 1, Knod
                   y = 0.d0
                   DO lz = 1, Lnod
                     DO ly = 1, Lnod
                             ! along y-dir                    along z-dir
                       tmp = GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                     ! xyz-coord of geom
                       y = y + tmp * ycoord(BoundaryNodeID(t2f(ly,lz),j))
                     ENDDO
                   ENDDO
                   BC_Psi(ky,kz,j) = (y - y_in_min)
                   BC_Psi(ky,kz,j) = 0.
                   BC_Switch(j) = NeumannBC
!                   BC_VelNorm(ky,kz,j) = 1.d0
                 ENDDO
               ENDDO
!               BC_Psi(1:K,j) = (0.25d0/HalfHeight) * (ycoord(BoundaryNodeID(2,j)) + ycoord(BoundaryNodeID(1,j)) + &
!                                          sps(1:K) * (ycoord(BoundaryNodeID(2,j)) - ycoord(BoundaryNodeID(1,j)) ) )
             ENDDO
           ENDIF

           IF (trim(BndryTag(i)) .eq. 'Outlet') THEN
             ! Set slip BC to not allow diffusion
             NoSlip(i) = .false.
             DO j = BndryNum(i-1)+1, BndryNum(i)
               DO kz = 1, Knod
                 DO ky = 1, Knod
                   y = 0.d0
                   DO lz = 1, Lnod
                     DO ly = 1, Lnod
                             ! along y-dir                    along z-dir
                       tmp = GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                     ! xyz-coord of geom
                       y = y + tmp * ycoord(BoundaryNodeID(t2f(ly,lz),j))
                     ENDDO
                   ENDDO
!                   BC_Psi(ky,kz,j) = 0.5*y
                   BC_Psi(ky,kz,j) = 0.
                   BC_Switch(j) = NeumannBC
                   BC_VelNorm(ky,kz,j) = 0.5d0
                 ENDDO
               ENDDO
!               BC_Psi(1:K,j) = (0.25d0/HalfHeight) * (ycoord(BoundaryNodeID(2,j)) + ycoord(BoundaryNodeID(1,j)) + &
!                                          sps(1:K) * (ycoord(BoundaryNodeID(2,j)) - ycoord(BoundaryNodeID(1,j)) ) )
             ENDDO
           ENDIF

         ENDDO

         BC_Switch_Lap = BC_Switch

       ELSE
         print *,'only prob_type = 3 (square cavity) is allowed '
         stop
       ENDIF

!       CLOSE(UNIT = 3)

END SUBROUTINE SetupICBCandSrc


SUBROUTINE Paraview_Output(nt)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer nt
       real*8 Reyn, fac, dxrat, dt

       integer el, kx, ky, kz, lx, ly, lz
       real*8 x, y, z, tmp

   INTERFACE

     FUNCTION itoa(i) RESULT(res)

       character(11) res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
     END FUNCTION itoa

   END INTERFACE


       OPEN(unit = 8, file = 'Vorticity'//trim(itoa(nt))//'.vtk', status = 'unknown')
       write(8,'(a)') '# vtk DataFile Version 3.0'
       write(8,'(a)') '3D Unstructured Grid of Quads'
       write(8,'(a)') 'ASCII'
       write(8,'(a)') ' '
!       write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'
       write(8,'(a)') 'DATASET POLYDATA'

       write(8,'(a,i8,a)') 'POINTS ', Knod**NDIM * Nel, ' float'
       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    along z-dir
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO
               write(8,*) real(x),real(y),real(z)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       write(8,'(a)') ' '
       write(8,'(a,i8,a)') 'POINT_DATA ', Knod**NDIM * Nel
       write(8,'(a)') 'SCALARS vorticity float 1'
       write(8,'(a)') 'LOOKUP_TABLE default'
       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               write(8,*) real(Vort(kx,ky,kz,el))
             ENDDO
           ENDDO
         ENDDO
       ENDDO

if(.true.) then
       write(8,'(a)') 'VECTORS velocity float'

       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               write(8,*) real(Uvel(kx,ky,kz,el)),real(Vvel(kx,ky,kz,el)),real(Wvel(kx,ky,kz,el))
             ENDDO
           ENDDO
         ENDDO
       ENDDO
endif
       CLOSE(unit = 8)

END SUBROUTINE Paraview_Output


SUBROUTINE dumpResult(numStep, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type, tIntegrator_type, numStep, prob_type
       real*8 Reyn, dt

       integer el, lx, ly, lz, kx, ky, kz, ntime
       real*8 x, y, z, tmp, time, timep, Vort_ex, Uvel_ex, Linf_Norm, L1_Norm, L2_Norm
       character(132) filename

       real*8 theta2,theta1,thetam,thetad,rad2,rad1,radm,radd
       real*8, dimension(Lnod,Lnod) :: xloc, yloc
       real*8, dimension(Lnod) :: xloc1D, yloc1D

   INTERFACE

     FUNCTION itoa(i) RESULT(res)

       character(11) res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
     END FUNCTION itoa

!     FUNCTION Vort_xct(x, y, z, t, Re, prob_type)
!       implicit NONE
!       real*8 x, y, z, t, Re, Vort_xct
!       integer prob_type
!     END FUNCTION Vort_xct

   END INTERFACE

       time = numStep*dt
       ntime = (time + 0.5d0*dt)
       IF (prob_type .le. 2) THEN
         timep = time - ntime
       ELSE
         timep = time
       ENDIF


       filename =  &
               'Problem_'//trim(itoa(prob_type))//'_Method_'//trim(itoa(HuynhSolver_type))//'_MeshOrder_'//trim(itoa(Knod))//      &
               '_TimeOrder_'//trim(itoa(tIntegrator_type))//'_NumMesh_'//trim(itoa(Nel))//'_dt_'//trim(itoa(INT(dt*1000000)))//    &
               '_nstep_'//trim(itoa(numStep))

       open(unit=9, file=trim(filename)//'.csv',status='unknown')

       Linf_Norm = -1.d0
       L1_Norm = 0.d0
       L2_Norm = 0.d0
       DO el = 1, Nel

         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod

               ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    x-coord of geom
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO

!               Vort_ex = Vort_xct(x, y, z, timep, Reyn, prob_type)
               Linf_Norm = Max(Linf_Norm, Abs(Vort_ex-Vort(kx,ky,kz,el)))
               L1_Norm = L1_Norm + Abs(Vort_ex-Vort(kx,ky,kz,el))  *wgt(kx)*wgt(ky)*wgt(kz)*Vol_Jac(kx,ky,kz,el)
               L2_Norm = L2_Norm + (Vort_ex-Vort(kx,ky,kz,el))**2  *wgt(kx)*wgt(ky)*wgt(kz)*Vol_Jac(kx,ky,kz,el)
!             write(*,*) x, y, Vort(i,kx,ky), Vort_ex, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
               write(9,*) x, y, Vort(kx,ky,kz,el), Vort_ex, abs(Vort_ex-Vort(kx,ky,kz,el)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
!             write(9,*) x, y, Vort(kx,ky,i), Uvel(kx,ky,i), Vvel(kx,ky,i)
             ENDDO
             write(9,*) ' '
           ENDDO
           write(9,*) ' '
         ENDDO

       ENDDO

       close(unit=9)

       open(unit=9, file=trim(filename)//'.stat',status='unknown')

       IF (.TRUE.) THEN

!       L1_Norm = L1_Norm / (Nel*Knod**3)
!       L2_Norm = Sqrt(L2_Norm / (Nel*Knod**3))
       L2_Norm = Sqrt(L2_Norm)
       write(9,*) trim('dt, 1.d0/Nel, Nel, Nel*Knod^3, Linf_Norm, L1_Norm, L2_Norm')
       write(9,*) dt, 1.d0/Nel, Nel, Nel*Knod**3, Linf_Norm, L1_Norm, L2_Norm
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod**3, Linf_Norm, L1_Norm, L2_Norm
       close(unit=9)

       ENDIF

       IF (.TRUE.) THEN
       open(unit=9, file=trim(filename)//'_du.csv',status='unknown')

       Linf_Norm = -1.d0
       L1_Norm = 0.d0
       L2_Norm = 0.d0
       DO el = 1, Nel

         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod

               ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    x-coord of geom
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO

!               Uvel_ex = Uvel_xct(x, y, z, timep, Reyn, prob_type)
               Linf_Norm = Max(Linf_Norm, Abs(Uvel_ex-Uvel(kx,ky,kz,el)))
               L1_Norm = L1_Norm + Abs(Uvel_ex-Uvel(kx,ky,kz,el))  *wgt(kx)*wgt(ky)*wgt(kz)*Vol_Jac(kx,ky,kz,el)
               L2_Norm = L2_Norm + (Uvel_ex-Uvel(kx,ky,kz,el))**2  *wgt(kx)*wgt(ky)*wgt(kz)*Vol_Jac(kx,ky,kz,el)
               write(9,*) x, y, Uvel(kx,ky,kz,el), Uvel_ex, abs(Uvel_ex-Uvel(kx,ky,kz,el)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
             ENDDO
             write(9,*) ' '
           ENDDO
           write(9,*) ' '
         ENDDO

       ENDDO

       close(unit=9)

       open(unit=9, file=trim(filename)//'_du.stat',status='unknown')

!       L1_Norm = L1_Norm / (Nel*Knod**2)
!       L2_Norm = Sqrt(L2_Norm / (Nel*Knod**2))
       L2_Norm = Sqrt(L2_Norm)
       write(9,*) trim('dt, 1.d0/Nel, Nel, Nel*Knod^2, Linf_Norm, L1_Norm, L2_Norm')
       write(9,*) dt, 1.d0/Nel, Nel, Nel*Knod**3, Linf_Norm, L1_Norm, L2_Norm
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod**3, Linf_Norm, L1_Norm, L2_Norm
       close(unit=9)

       ENDIF

  RETURN

       open(unit=9, file=trim(filename)//'_boundary.csv',status='unknown')

       Linf_Norm = -1.d0
       L1_Norm = 0.d0
       L2_Norm = 0.d0
       DO el = 1, Nelb/2

         DO lx = 1, Lnod
           xloc1D(lx) = xcoord(boundarynodeID(t2f(lx),el))
           yloc1D(lx) = ycoord(boundarynodeID(t2f(lx),el))
         ENDDO

         DO ky = 1, Knod
         DO kx = 1, Knod
           ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
           x = 0.d0
           y = 0.d0
           DO lx = 1, Lnod
                     ! along x-dir                  x-coord of geom
             x = x + GeomNodesLgrangeBasis(lx,kx) * xloc1D(lx)
             y = y + GeomNodesLgrangeBasis(lx,kx) * yloc1D(lx)
           ENDDO

!             Vort_ex = Vort_xct(x, y, z, timep, Reyn, prob_type)
             Linf_Norm = Max(Linf_Norm, Abs(Vort_ex-VelocJump(kx,ky,el)))
             L1_Norm = L1_Norm + Abs(Vort_ex-VelocJump(kx,ky,el))
             L2_Norm = L2_Norm + (Vort_ex-VelocJump(kx,ky,el))**2
!             write(*,*) x, y, Vort(i,kx,ky), Vort_ex, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
             write(9,*) x, y, VelocJump(kx,ky,el), Vort_ex, abs(Vort_ex-VelocJump(kx,ky,el)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
!             write(9,*) x, y, Vort(kx,ky,i), Uvel(kx,ky,i), Vvel(kx,ky,i)
           ENDDO
           ENDDO
           write(9,*) ' '
         ENDDO
       L1_Norm = L1_Norm / (Nel*Knod)
       L2_Norm = Sqrt(L2_Norm / (Nel*Knod))
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod, Linf_Norm, L1_Norm, L2_Norm

       close(unit=9)

       open(unit=9, file=trim(filename)//'_boundary_du.csv',status='unknown')

       Linf_Norm = -1.d0
       L1_Norm = 0.d0
       L2_Norm = 0.d0
       DO el = 1, Nelb

         DO lx = 1, Lnod
           xloc1D(lx) = xcoord(boundarynodeID(t2f(lx),el))
           yloc1D(lx) = ycoord(boundarynodeID(t2f(lx),el))
         ENDDO

         DO ky = 1, Knod
         DO kx = 1, Knod
           ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
           x = 0.d0
           y = 0.d0
           DO lx = 1, Lnod
                     ! along x-dir                  x-coord of geom
             x = x + GeomNodesLgrangeBasis(lx,kx) * xloc1D(lx)
             y = y + GeomNodesLgrangeBasis(lx,kx) * yloc1D(lx)
           ENDDO

!             Uvel_ex = Uvel_xct(x, y, z, timep, Reyn, prob_type)
             Linf_Norm = Max(Linf_Norm, Abs(Uvel_ex-VelocJump(kx,ky,el)))
             L1_Norm = L1_Norm + Abs(Uvel_ex-VelocJump(kx,ky,el))
             L2_Norm = L2_Norm + (Uvel_ex-VelocJump(kx,ky,el))**2
             write(9,*) x, y, VelocJump(kx,ky,el), Uvel_ex, abs(Uvel_ex-VelocJump(kx,ky,el))
           ENDDO
           write(9,*) ' '
         ENDDO
         ENDDO
       L1_Norm = L1_Norm / (Nel*Knod)
       L2_Norm = Sqrt(L2_Norm / (Nel*Knod))
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod, Linf_Norm, L1_Norm, L2_Norm

       close(unit=9)

END SUBROUTINE dumpResult


SUBROUTINE Paraview_Dump(HuynhSolver_type, nt, t)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer HuynhSolver_type, nt
       real*8 t

       integer el, kx, ky, kz, lx, ly, lz
       real*8 x, y, z, tmp
       real*8, dimension(Knod, Knod, Knod, Nel) :: psi

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

   INTERFACE

     FUNCTION itoa(i) RESULT(res)

       character(11) res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
     END FUNCTION itoa

   END INTERFACE


       OPEN(unit = 8, file = 'Vorticity'//trim(itoa(nt))//'.vtk', status = 'unknown')

       write(8,'(a)') '# vtk DataFile Version 3.0'
       write(8,'(a)') '3D Unstructured Grid of Quads'
       write(8,'(a)') 'ASCII'
       write(8,'(a)') ' '
!       write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'
       write(8,'(a)') 'DATASET POLYDATA'

       write(8,'(a,i8,a)') 'POINTS ', Knod**NDIM*Nel, ' float'
       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
               x = 0.d0
               y = 0.d0
               z = 0.d0
               DO lz = 1, Lnod
                 DO ly = 1, Lnod
                   DO lx = 1, Lnod
                           ! along x-dir                  along y-dir                    along z-dir
                     tmp = GeomNodesLgrangeBasis(lx,kx) * GeomNodesLgrangeBasis(ly,ky) * GeomNodesLgrangeBasis(lz,kz)
                                   ! xyz-coord of geom
                     x = x + tmp * xcoord(nodeID(t2f(lx,ly,lz),el))
                     y = y + tmp * ycoord(nodeID(t2f(lx,ly,lz),el))
                     z = z + tmp * zcoord(nodeID(t2f(lx,ly,lz),el))
                   ENDDO
                 ENDDO
               ENDDO
               write(8,*) real(x),real(y),real(z)
   Vort0(kx,ky,kz,el) = (1.5d0*Cos(pi*x)*exp(-pi**2*t) - 0.5d0*Cos(3.d0*pi*x)*exp(-9.d0*pi**2*t)) * &
                        (1.5d0*Cos(pi*y)*exp(-pi**2*t) - 0.5d0*Cos(3.d0*pi*y)*exp(-9.d0*pi**2*t)) * &
                        (1.5d0*Cos(pi*z)*exp(-pi**2*t) - 0.5d0*Cos(3.d0*pi*z)*exp(-9.d0*pi**2*t))
             ENDDO
           ENDDO
         ENDDO
       ENDDO

 if(.false.) then
       write(8,'(a)') ' '
       write(8,'(a,i8,1x,i8)') 'CELLS ', Knod**NDIM*Nel, 2*Knod**NDIM*Nel
       DO el = 1, Knod**NDIM*Nel
!!****  We need to change "1" to something else for 3D
         write(8,*) "1",el-1
       ENDDO
       write(8,'(a)') ' '
       write(8,'(a,i8)') 'CELL_TYPES ', Knod**NDIM*Nel
       DO el = 1, Knod**NDIM*Nel
!!****  We need to change "1" to something else for 3D
         write(8,*)  "1"
       ENDDO
endif

       write(8,'(a)') ' '
       write(8,'(a,i8,a)') 'POINT_DATA ', Knod**NDIM * Nel
       write(8,'(a)') 'SCALARS exact float'
       write(8,'(a)') 'LOOKUP_TABLE default'
       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               write(8,*) real(Vort0(kx,ky,kz,el))
             ENDDO
           ENDDO
         ENDDO
       ENDDO
!       write(8,'(a,i8,a)') 'POINT_DATA ', Knod**NDIM*Nel
       write(8,'(a)') 'SCALARS vorticity float'
       write(8,'(a)') 'LOOKUP_TABLE default'
       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               write(8,*) real(Vort(kx,ky,kz,el))
             ENDDO
           ENDDO
         ENDDO
       ENDDO

if(.false.) then
       BC_Switch = DirichletBC
       BC_Values = BC_Psi
!       CALL GetLaplacian(HuynhSolver_type, Vort, psi, A_handle, P_handle, S_handle)                                ! Stage 1
!       CALL GetLaplacGrads(HuynhSolver_type, psi, 1)

!endif

       write(8,'(a)') 'VECTORS velocity float'

       DO el = 1, Nel
         DO kz = 1, Knod
           DO ky = 1, Knod
             DO kx = 1, Knod
               write(8,*) real(Uvel(kx,ky,kz,el)),real(Vvel(kx,ky,kz,el)),real(Wvel(kx,ky,kz,el))
             ENDDO
           ENDDO
         ENDDO
       ENDDO
endif

       CLOSE(unit = 8)


END SUBROUTINE Paraview_Dump


FUNCTION itoa(i) RESULT(res)

       character(11) res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
       write(tmp,'(i0)') i
!       res = trim(tmp)
       res = tmp

END FUNCTION itoa


SUBROUTINE SetLaplacian(HuynhSolver_type, A_handle, P_handle, S_handle)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer HuynhSolver_type

       integer jx, jy, jz, lx, ly, lz, el, eln, bel, ijP, ijPm, idir, idir2, idir3, ibnd, ibndm
       integer i, j, k, l, m, p, q, ij, jk, jl, jm, lk, mk, mp, pl, pq
       integer jm1, km1, lm1, mm1, qm1, km1s, mm1s, pm1s, jk0, jm0, mk0
       integer ijk, ijm, imk, mjk, ilm, lmk, lmp, ljm, mpq
       integer Knm1, Knsq, KnsqmK, Kmxsq, ierr, Kmax
       integer colctr, colnbr(0:FC), nnz, nrows

       real*8 tmp, tmpB, tmpC, tmpx, tmpy, tmpz, Ovr_Jac
       real*8, dimension(NDIM, NDIM) :: Dxsi_iDx_j
       real*8, dimension(Knod, Knod) :: Acoef, Bcoef, Ccoef
       real*8, dimension(Knod, Knod, 0:FC) :: FaceA, FaceB, FaceC, NormA
       real*8, dimension(Knod, Knod, 0:1, 0:1) :: SBLB_i_NGR_j, SBGLB_SBLB_i_NGR_j
       real*8, dimension(Knod, Knod) :: SNGLB_SBLBdotNGR, SNGLB_2xSBLBdotNGR
       real*8, dimension(Knod, Knod, Knod, NDIM) :: SolnA, SolnB, SolnC
       ! small dense matrix to obtain comx (per element) for a given Neumann BC
       real*8, dimension(Knod**(NDIM-1), Knod**(NDIM-1)) :: NeuMatrix_Orig, NeuMatrix

       real*8, allocatable, dimension(:,:,:) :: LaplacianCenter
       real*8, allocatable, dimension(:,:,:,:) :: LaplacianNeigbr

       integer, dimension(:), allocatable :: rowptr, colidx
       real*8, dimension(:), allocatable :: values

       real*8, dimension(:), allocatable :: b, x

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

       integer, dimension(:), allocatable :: rowptr_filtered, colidx_filtered
       real*8, dimension(:), allocatable :: values_filtered

   INTERFACE

     PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

     END FUNCTION Cross_Product

     FUNCTION inv(A) RESULT(Ainv)
       integer dp
       parameter (dp=8)
       real(dp), dimension(:,:), intent(in) :: A
       real(dp), dimension(size(A,1),size(A,2)) :: Ainv

!       real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
!       integer, dimension(size(A,1)) :: ipiv   ! pivot indices
!       integer :: n, info
     END FUNCTION inv

   END INTERFACE

       Knm1 = Knod - 1
       Knsq = Knod * Knod
       KnsqmK = Knsq - Knod
       Kmax = Knod**NDIM
       Kmxsq = Kmax * Kmax
       allocate (LaplacianCenter(Kmax,Kmax,Nel), LaplacianNeigbr(Kmax,Kmax,0:FC,Nel))

!         gLprime(0) = -0.25d0*Knod*Knod
!         gLprime(1) =  0.25d0*Knod*knod
       gLprime(0) = -0.25d0*Knod*(Knod+1)
       gLprime(1) =  0.25d0*Knod*(Knod+1)

       BndrySrc = 0.d0
       LaplacianCenter = 0.d0
       LaplacianNeigbr = 0.d0

       DO i = 1, Knod
         DO j = 1, Knod
           ! Note: SBLB_i0_NGR_j0(i , j) + SBLB_i1_NGR_j1(Knod+1-i , Knod+1-j) = 0
           !       SBLB_i0_NGR_j1(i , j) + SBLB_i1_NGR_j0(Knod+1-i , Knod+1-j) = 0
           !       But we don't take advantage of this, if there is any

           ! E00(i,j), E01(i,j), E11(i,j), E10(i,j) in notes
           DO l = 0, 1
             DO m = 0, 1
               ! SBLB_i_NGR_j(i,j,l,m) is Elm(i,j) in notes (E00(i,j), E01(i,j), E11(i,j), E10(i,j))
               SBLB_i_NGR_j(i,j,l,m) = 0.5d0 * SolnBndryLgrangeBasis(i,l) * NodesGradRadau(j,m)

               ! SBGLB_SBLB_i_NGR_j(i,j,l,m) is Flm(i,j) in notes (F00(i,j), F01(i,j), F11(i,j), F10(i,j))
               SBGLB_SBLB_i_NGR_j(i,j,l,m) = 0.5d0 * (SolnBndryGradLgrangeBasis(i,l) - gLprime(l)*SolnBndryLgrangeBasis(i,l)) * &
                                                      NodesGradRadau(j,m)
             ENDDO
           ENDDO
           ! C(i,j) in notes;  C(i,i) = 0  (should use this during matrix setup)
           SNGLB_SBLBdotNGR(i,j) = SolnNodesGradLgrangeBasis(i,j) - (SBLB_i_NGR_j(i,j,0,0) + SBLB_i_NGR_j(i,j,1,1))
           ! D(i,j) in notes;  D(Knod/2+1 , Knod/2+1) = 0 for Knod = 1, 3, 5  (should use this during matrix setup)
           SNGLB_2xSBLBdotNGR(i,j) = SolnNodesGradLgrangeBasis(i,j) - 2.d0*(SBLB_i_NGR_j(i,j,0,0) + SBLB_i_NGR_j(i,j,1,1))
         ENDDO
       ENDDO

       DO el = 1, Nel

         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
               Ovr_Jac = 1.d0 / Vol_Jac(jx,jy,jz,el)

               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el))
               ! Jac * dzta / d_j
               Dxsi_iDx_j(:,3) = Cross_Product(Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el))

               ! Axx, Axy, Axz in notes
               SolnA(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnB(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,2)) * Ovr_Jac
               SolnC(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,3)) * Ovr_Jac

               ! Ayy in notes
               ! Saving in cyclic order makes it possible to lump operations in x, y, and z directions in one uniform call
! *** May actually have to do jy,jx,jz for better performance (looks that way)
               SolnA(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,2)) * Ovr_Jac
               SolnB(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnC(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,3)) * Ovr_Jac

               ! Azz in notes
               SolnA(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,3)) * Ovr_Jac
               SolnB(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnC(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,2)) * Ovr_Jac
             ENDDO

             DO ijP = 0, FC
               FaceA(jy,jz,ijP) = Face_Acoef(jy,jz,ijP,el) / Face_Jac(jy,jz,ijP,el)
               ! FaceB = FaceC = 0 for orthogonal (rectangular?) grids - should take advantage of it during matrix setup
               FaceB(jy,jz,ijP) = Face_Bcoef(jy,jz,ijP,el) / Face_Jac(jy,jz,ijP,el)
               FaceC(jy,jz,ijP) = Face_Ccoef(jy,jz,ijP,el) / Face_Jac(jy,jz,ijP,el)
               NormA(jy,jz,ijP) = Face_Norm(jy,jz,ijP,el)
             ENDDO
           ENDDO
         ENDDO

         DO k = 1, Knod
           km1s = (k-1) * Knsq
           DO j = 1, Knod
             jm1 = (j-1) * Knod
             jk = jm1 + km1s
             DO i = 1, Knod
                ij = i + jm1
               ijk = i + jk
               DO m = 1, Knod
                 mm1  = (m-1) * Knod
                 mm1s = mm1 * Knod
                 mjk  = m + jk
                  mk  = mm1 + km1s
                 imk = i + mk
                  jm = jm1 + mm1s
                 ijm = i + jm
                 tmpx = 0.d0
                 tmpy = 0.d0
                 tmpz = 0.d0
                 DO l = 1, Knod
                   ! Axx * D * C
                   tmpx = tmpx + SolnA(l,j,k,1) * SNGLB_2xSBLBdotNGR(l,i) * SNGLB_SBLBdotNGR(m,l)
                   tmpy = tmpy + SolnA(l,i,k,2) * SNGLB_2xSBLBdotNGR(l,j) * SNGLB_SBLBdotNGR(m,l)
                   tmpz = tmpz + SolnA(l,i,j,3) * SNGLB_2xSBLBdotNGR(l,k) * SNGLB_SBLBdotNGR(m,l)
                 ENDDO
                 LaplacianCenter(mjk,ijk,el) = LaplacianCenter(mjk,ijk,el) + tmpx
                 LaplacianCenter(imk,ijk,el) = LaplacianCenter(imk,ijk,el) + tmpy
                 LaplacianCenter(ijm,ijk,el) = LaplacianCenter(ijm,ijk,el) + tmpz

 if(.true.) then
                 DO l = 1, Knod
                   lmk = l + mk
                   ljm = l + jm
                   ilm = i + (l-1) * Knod + mm1s
                   LaplacianCenter(lmk,ijk,el) = LaplacianCenter(lmk,ijk,el) + ( &
                                                 SolnB(l,j,k,1) * SNGLB_2xSBLBdotNGR(l,i) * SNGLB_SBLBdotNGR(m,j) + &
                                                 SolnB(m,i,k,2) * SNGLB_2xSBLBdotNGR(m,j) * SNGLB_SBLBdotNGR(l,i) )
                   LaplacianCenter(ljm,ijk,el) = LaplacianCenter(ljm,ijk,el) + ( &
                                                 SolnC(l,j,k,1) * SNGLB_2xSBLBdotNGR(l,i) * SNGLB_SBLBdotNGR(m,k) + &
                                                 SolnB(m,i,j,3) * SNGLB_2xSBLBdotNGR(m,k) * SNGLB_SBLBdotNGR(l,i) )
                   LaplacianCenter(ilm,ijk,el) = LaplacianCenter(ilm,ijk,el) + ( &
                                                 SolnC(l,i,k,2) * SNGLB_2xSBLBdotNGR(l,j) * SNGLB_SBLBdotNGR(m,k) + &
                                                 SolnC(m,i,j,3) * SNGLB_2xSBLBdotNGR(m,k) * SNGLB_SBLBdotNGR(l,j) )
                 ENDDO
 endif
               ENDDO
             ENDDO
           ENDDO
         ENDDO

         ijP = 0
         ! Extrapolation operations in x (=1), y (=2), and z (=3) directions
         DO idir = 1, NDIM
           idir2 = Min(idir-1,1) * Knm1
           idir3 = Max(idir-2,0) * KnsqmK

           ! Operations at the Left/South/Bottom (=0) and Right/North/Top (=1) boundaries
           DO ibnd = 0, 1
             ibndm = 1 - ibnd
             ! mesh to the left/south/bottom of left/south/bottom face or right/north/top of right/north/top face
             eln = elemID(i2f(ijP),el)

             IF (eln .lt. 0) THEN

               bel = -eln
               IF (BC_Switch(bel) .eq. DirichletBC) THEN

                 DO k = 1, Knod
                   km1s = (k-1) * Knsq
                   DO j = 1, Knod
                     jm1 = (j-1) * Knod
                     jk = jm1 + km1s
!                    tmp = 2.d0 * gLprime(ibnd) * FaceA(j,ijP) * BC_Values(j,bel)
                     tmp = 2.d0 * gLprime(ibnd) * FaceA(j,k,ijP)
                     DO i = 1, Knod
                       ! idir2*(i-j) basically saves the info in the transposed position of the matrix for idir = 2
                       ! idir2*(i-j) + idir3*(i-k) does this for idir = 3
                       !    This allows lumping all operations for x, y and z directions into one generic loop
                       ijk = i + jk + idir2*(i-j) + idir3*(i-k)
                       BndrySrc(j,k,ijk,bel) = BndrySrc(j,k,ijk,bel) + tmp * NodesGradRadau(i,ibnd)
                       DO m = 1, Knod
                         mm1 = (m-1) * Knod
                         mk = mm1 + km1s
                         mjk = m + jk + idir2*(m-j) + idir3*(m-k)
                         jm = jm1 + mm1*Knod
!                      BndrySrc(ij,el) = BndrySrc(ij,el) + &
!                                        SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,ibnd) &
!                                                                                 * BC_Values(j,bel) - &
!                                       NodesGradRadau(i,ibnd) * BC_Values(m,bel) &
!                                                           * ( SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) + &
!                                                               FaceB(j,el,ijP) * SolnNodesGradLgrangeBasis(m,j) )
                         BndrySrc(j,k,ijk,bel) = BndrySrc(j,k,ijk,bel) + &
                                                          SolnA(m,j,k,idir) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,ibnd)
 if(.true.) then
! check out SolnB; for idir=1 (and i,m,k) we really want idir=2 (and m,i,k) ... this will have to wait till later
                         BndrySrc(m,k,ijk,bel) = BndrySrc(m,k,ijk,bel) + NodesGradRadau(i,ibnd) &
                                                            * ( SolnB(i,m,k,idir) * SNGLB_2xSBLBdotNGR(m,j) + &
                                                                FaceB(j,k,ijP) * SolnNodesGradLgrangeBasis(m,j) )

! check out SolnC; for idir=1 (and i,j,m) we really want idir=3 (and m,i,j) ... this will have to wait till later
                         BndrySrc(j,m,ijk,bel) = BndrySrc(j,m,ijk,bel) + NodesGradRadau(i,ibnd) &
                                                            * ( SolnC(i,j,m,idir) * SNGLB_2xSBLBdotNGR(m,k) + &
                                                                FaceC(j,k,ijP) * SolnNodesGradLgrangeBasis(m,k) )
 endif
                         tmpx = 0.d0
                         DO l = 1, Knod
                           tmpx = tmpx + SolnA(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,l,ibnd,ibnd)
                         ENDDO
                         LaplacianCenter(mjk,ijk,el) = LaplacianCenter(mjk,ijk,el) - tmpx + 2.d0 * FaceA(j,k,ijP) * &
                                                     ( SBGLB_SBLB_i_NGR_j(m,i,ibnd,ibnd) - &
                                                       gLprime(ibnd) * SBLB_i_NGR_j(m,i,ibnd,ibnd) )

 if(.true.) then
                         tmpy = SolnB(i,m,k,idir) * SNGLB_2xSBLBdotNGR(m,j)
                         DO l = 1, Knod
                           lmk = l + mk + idir2*(l-m) + idir3*(l-k)
                           LaplacianCenter(lmk,ijk,el) = LaplacianCenter(lmk,ijk,el) - tmpy * SBLB_i_NGR_j(l,i,ibnd,ibnd)
                         ENDDO

                         tmpz = SolnC(i,j,m,idir) * SNGLB_2xSBLBdotNGR(m,k)
                         DO l = 1, Knod
                           ljm = l + jm  + idir2*(l-j) + idir3*(l-m)
                           LaplacianCenter(ljm,ijk,el) = LaplacianCenter(ljm,ijk,el) - tmpz * SBLB_i_NGR_j(l,i,ibnd,ibnd)
                         ENDDO
 endif

                       ENDDO
                     ENDDO
                   ENDDO
                 ENDDO

               ELSEIF (BC_Switch(bel) .eq. NeumannBC) THEN

                 NeuMatrix_Orig = 0.d0
                 DO l = 1, Knod
                   lm1 = (l-1) * Knod
                   DO k = 1, Knod
                     km1 = (k-1) * Knod
                     lk = l + km1
                     DO j = 1, Knod
                       jl = j + lm1
                       jk = j + km1
                       NeuMatrix_Orig(jk,lk) = NeuMatrix_Orig(jk,lk) + FaceB(j,k,ijP) * SolnNodesGradLgrangeBasis(l,j)
                       NeuMatrix_Orig(jk,jl) = NeuMatrix_Orig(jk,jl) + FaceC(j,k,ijP) * SolnNodesGradLgrangeBasis(l,k)
                     ENDDO
                     NeuMatrix_Orig(lk,lk) = NeuMatrix_Orig(lk,lk) + 2.d0 * gLprime(ibnd) * FaceA(l,k,ijP)
                   ENDDO
                 ENDDO

                 NeuMatrix = inv(NeuMatrix_Orig)

                 DO k = 1, Knod
                   km1  = (k-1) * Knod
                   km1s = (k-1) * Knsq
                   DO j = 1, Knod
                     jm1 = (j-1) * Knod
                     jk = jm1 + km1s
                     DO i = 1, Knod
                       ijk = i + jk + idir2*(i-j) + idir3*(i-k)
!                    BndrySrc(ij,el) = BndrySrc(ij,el) + NodesGradRadau(i,ibnd) * NormA(j,ijP) * BC_Values(j,bel)
                       BndrySrc(j,k,ijk,bel) = BndrySrc(j,k,ijk,bel) + NodesGradRadau(i,ibnd) * NormA(j,k,ijP)
                       DO m = 1, Knod
                         mm1 = (m-1) * Knod
                         mk = mm1 + km1s
                         mjk = m + jk + idir2*(m-j) + idir3*(m-k)
                         jm = jm1 + mm1*Knod

!                     tmpx = 0.d0
!                     tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(i,ibnd)
!                     DO l = 1, Knod
!                       tmpx = tmpx + NeuMatrix(j,l) * NormA(l,ijP) * BC_Values(l,bel)
!                       BndrySrc(ij,el) = BndrySrc(ij,el) + tmpy * NeuMatrix(m,l) * NormA(l,ijP) * BC_Values(l,bel)
!                     ENDDO
!                     BndrySrc(ij,el) = BndrySrc(ij,el) - tmpx * SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,ibnd)
                         tmpx = SolnA(m,j,k,idir) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,ibnd)
                         tmpy = SolnB(i,m,k,idir) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(i,ibnd)
                         tmpz = SolnC(i,j,m,idir) * SNGLB_2xSBLBdotNGR(m,k) * NodesGradRadau(i,ibnd)
                         jk0 = j + km1
                         jm0 = j + mm1
                         mk0 = m + km1
                         DO l = 1, Knod
                           lm1 = (l-1) * Knod
                           DO p = 1, Knod
                             pl = p + lm1
                             BndrySrc(p,l,ijk,bel) = BndrySrc(p,l,ijk,bel) + NormA(p,l,ijP) * &
                                                     (tmpx*NeuMatrix(jk0,pl) + tmpy*NeuMatrix(mk0,pl) + tmpz*NeuMatrix(jm0,pl))

                           ENDDO
                         ENDDO

                         tmpx = 0.d0
                         DO l = 1, Knod
                           tmpx = tmpx + SolnA(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,l,ibnd,ibnd)
                         ENDDO
                         LaplacianCenter(mjk,ijk,el) = LaplacianCenter(mjk,ijk,el) - tmpx

                         tmpB = SolnB(i,m,k,idir) * SNGLB_2xSBLBdotNGR(m,j)
                         tmpC = SolnC(i,j,m,idir) * SNGLB_2xSBLBdotNGR(m,k)
                         DO l = 1, Knod
                           lmk = l + mk + idir2*(l-m) + idir3*(l-k)
                           LaplacianCenter(lmk,ijk,el) = LaplacianCenter(lmk,ijk,el) - tmpB * SBLB_i_NGR_j(l,i,ibnd,ibnd)
                           ljm = l + jm + idir2*(l-j) + idir3*(l-m)
                           LaplacianCenter(ljm,ijk,el) = LaplacianCenter(ljm,ijk,el) - tmpC * SBLB_i_NGR_j(l,i,ibnd,ibnd)
                         ENDDO

                         tmpx = 0.d0
                         DO l = 1, Knod
                           tmpx = tmpx + SolnA(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * &
                                         (SBGLB_SBLB_i_NGR_j(m,l,ibnd,ibnd) - gLprime(ibnd) * SBLB_i_NGR_j(m,l,ibnd,ibnd))
                         ENDDO
                         tmpx = -2.d0 * tmpx

                         tmp = -2.d0 * (SBGLB_SBLB_i_NGR_j(m,i,ibnd,ibnd) - gLprime(ibnd) * SBLB_i_NGR_j(m,i,ibnd,ibnd))

!!
!                         tmpy =0.d0
!                         DO l = 1, Knod
!                           tmpy = tmpy + SolnB(i,l,k,idir) * SNGLB_2xSBLBdotNGR(l,j) * NeuMatrix(l,m)
!                         ENDDO
!                         tmpy = 2.d0 * tmpy * FaceA(m,ijP)
!!

                         DO q = 1, Knod
                           qm1 = (q - 1) * Knod
                           DO p = 1, Knod
                             pq = p + qm1
                             tmpy = 0.d0
                             tmpz = 0.d0
                             DO l = 1, Knod
                               lk = l + km1
                               jl = j + (l - 1) * Knod
                               tmpy = tmpy + SolnB(i,l,k,idir) * SNGLB_2xSBLBdotNGR(l,j) * NeuMatrix(lk,pq)
                               tmpz = tmpz + SolnC(i,j,l,idir) * SNGLB_2xSBLBdotNGR(l,k) * NeuMatrix(jl,pq)
                             ENDDO

!                             lmk = l + mk + idir2*(l-m) + idir3*(l-k)
                             mpq = m + (pq - 1) * Knod + idir2*(m-p) + idir3*(m-q)
                             LaplacianCenter(mpq,ijk,el) = LaplacianCenter(mpq,ijk,el) + FaceA(p,q,ijP) * &
                                                           (tmpx * NeuMatrix(jk0,pq) + tmp * (tmpy + tmpz))
!                                                           (tmpx * NeuMatrix(pq,jk0) + tmp * (tmpy + tmpz))

!                             tmpx = tmpx + SolnA(p,j,idir) * SNGLB_2xSBLBdotNGR(p,i) * &
!                                                          (SBGLB_SBLB_i_NGR_j(l,p,ibnd,ibnd) - &
!                                                           gLprime(ibnd) * SBLB_i_NGR_j(l,p,ibnd,ibnd))
                           ENDDO
!                           lmk = l + mk + idir2*(l-m) + idir3*(l-k)
!                           LaplacianCenter(lmk,ijk,el) = LaplacianCenter(lmk,ijk,el)  &
!                                                                               + 2.d0 * tmpx * FaceA(m,ijP) * NeuMatrix(j,m) &
!                                                                               - tmpy * (SBGLB_SBLB_i_NGR_j(l,i,ibnd,ibnd)   &
!                                                                               - gLprime(ibnd) * SBLB_i_NGR_j(l,i,ibnd,ibnd))
                         ENDDO
                       ENDDO
                     ENDDO
                   ENDDO
                 ENDDO

               ENDIF

             ELSE

               DO k = 1, Knod
                 DO j = 1, Knod
                   ijPm = nbr(ijP)
                   Acoef(j,k) = Face_Acoef(j,k,ijPm,eln) / Face_Jac(j,k,ijPm,eln)
                   Bcoef(j,k) = 0.5d0 * ( FaceB(j,k,ijP) + Face_Bcoef(j,k,ijPm,eln) / Face_Jac(j,k,ijPm,eln) )
                   Ccoef(j,k) = 0.5d0 * ( FaceC(j,k,ijP) + Face_Ccoef(j,k,ijPm,eln) / Face_Jac(j,k,ijPm,eln) )
                 ENDDO
               ENDDO

               DO k = 1, Knod
                 km1s = (k-1) * Knsq
                 DO j = 1, Knod
                   jm1 = (j-1) * Knod
                   jk = jm1 + km1s
                   DO i = 1, Knod
                     ijk = i + jk + idir2*(i-j) + idir3*(i-k)
                     DO m = 1, Knod
                       mm1 = (m-1) * Knod
                       mk = mm1 + km1s
                       mjk = m + jk + idir2*(m-j) + idir3*(m-k)
                       jm = jm1 + mm1*Knod
                       tmpx = 0.d0
                       DO l = 1, Knod
                                       ! Axx * D * E01
                         tmpx = tmpx + SolnA(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,l,ibndm,ibnd)
                       ENDDO
                       LaplacianNeigbr(mjk,ijk,ijP,el) = LaplacianNeigbr(mjk,ijk,ijP,el) + tmpx + &
                                                       Acoef(j,k) * SBGLB_SBLB_i_NGR_j(m,i,ibndm,ibnd) + &
                                                       gLprime(ibnd) * FaceA(j,k,ijP) * SBLB_i_NGR_j(m,i,ibndm,ibnd)
                       ! Ae(j) * F11(m,i) + gL(w) * A(j) * E11(m,i)
                       LaplacianCenter(mjk,ijk,el) = LaplacianCenter(mjk,ijk,el) + ( &
                                                   FaceA(j,k,ijP) * SBGLB_SBLB_i_NGR_j(m,i,ibnd,ibnd) + &
                                                   gLprime(ibndm) * Acoef(j,k) * SBLB_i_NGR_j(m,i,ibnd,ibnd) )

 if(.true.) then
                       tmpB = Bcoef(j,k) * SolnNodesGradLgrangeBasis(m,j)
                       tmpC = Ccoef(j,k) * SolnNodesGradLgrangeBasis(m,k)
                       DO l = 1, Knod
                         lmk = l + mk + idir2*(l-m) + idir3*(l-k)
                         LaplacianCenter(lmk,ijk,el) = LaplacianCenter(lmk,ijk,el) + tmpB * SBLB_i_NGR_j(l,i,ibnd,ibnd)
                         LaplacianNeigbr(lmk,ijk,ijP,el) = LaplacianNeigbr(lmk,ijk,ijP,el) + tmpB * SBLB_i_NGR_j(l,i,ibndm,ibnd)

                         ljm = l + jm  + idir2*(l-j) + idir3*(l-m)
                         LaplacianCenter(ljm,ijk,el) = LaplacianCenter(ljm,ijk,el) + tmpC * SBLB_i_NGR_j(l,i,ibnd,ibnd)
                         LaplacianNeigbr(ljm,ijk,ijP,el) = LaplacianNeigbr(ljm,ijk,ijP,el) + tmpC * SBLB_i_NGR_j(l,i,ibndm,ibnd) 
                       ENDDO

                       ! Temporary quick solution; problem is for y-dir a we normally don't use the cylic method but in this case 
                       !                                   a cyclic method is required for the general method to work
                       !                           Need to look carefully for a general cyclic method for everything
                       IF (idir .eq. 1) THEN
                         DO l = 1, Knod
                           lmk = l + mk !+ idir2*(l-m) + idir3*(l-k)
                           LaplacianNeigbr(lmk,ijk,2+ibnd,el) = LaplacianNeigbr(lmk,ijk,2+ibnd,el) + &
                                        SolnB(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,j,ibndm,ibnd) 

                           ljm = l + jm !+ idir2*(l-j) + idir3*(l-m)
                           LaplacianNeigbr(ljm,ijk,4+ibnd,el) = LaplacianNeigbr(ljm,ijk,4+ibnd,el) + &
                                        SolnC(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,k,ibndm,ibnd) 
                         ENDDO
                       ELSEIF (idir .eq. 2) THEN
                         DO l = 1, Knod
                           lmk = l + mk + idir2*(l-m) !+ idir3*(l-k)
                           LaplacianNeigbr(lmk,ijk,ibnd,el) = LaplacianNeigbr(lmk,ijk,ibnd,el) + &
                                        SolnB(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,j,ibndm,ibnd) 

                           ljm = l + jm + idir2*(l-j) !+ idir3*(l-m)
                           LaplacianNeigbr(ljm,ijk,4+ibnd,el) = LaplacianNeigbr(ljm,ijk,4+ibnd,el) + &
                                        SolnC(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,k,ibndm,ibnd) 
                         ENDDO
                       ELSEIF (idir .eq. 3) THEN
                         DO l = 1, Knod
                           lmk = l + mk + idir2*(l-m) + idir3*(l-k)
                           LaplacianNeigbr(lmk,ijk,ibnd,el) = LaplacianNeigbr(lmk,ijk,ibnd,el) + &
                                        SolnB(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,j,ibndm,ibnd) 

                           ljm = l + jm + idir2*(l-j) + idir3*(l-m)
                           LaplacianNeigbr(ljm,ijk,2+ibnd,el) = LaplacianNeigbr(ljm,ijk,2+ibnd,el) + &
                                        SolnC(l,j,k,idir) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i_NGR_j(m,k,ibndm,ibnd)  
                         ENDDO
                       ENDIF
 endif
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO

             ENDIF
             ijP = ijP + 1
           ENDDO
         ENDDO

       ENDDO

!      Everything above can be used in conjunction with another Laplacian to get fast diffusion; for now, below is to solve
!      the Poisson equation, but it needs to be modified to replace Laplacian (for diffusion) variables into matrix A
       nrows = Kmax * Nel
       nnz = Kmxsq * Nel
       DO el = 1, Nel
         DO ijP = 0, FC
           IF (elemID(i2f(ijP),el) > 0) nnz = nnz + Kmxsq
         ENDDO
       ENDDO

       Allocate (rowptr(nrows+1), colidx(nnz), values(nnz))

       rowptr(nrows+1) = nnz + 1

       nrows = 0
       nnz = 0
       DO el = 1, Nel
         colctr = (el - 1) * Kmax
         DO ijP = 0, FC
           colnbr(ijP) = (elemID(i2f(ijP),el) - 1) * Kmax
         ENDDO
         DO k = 1, Knod
           km1s = (k-1) * Knsq
           DO j = 1, Knod
             jm1 = (j-1) * Knod
              jk = jm1 + km1s
             DO i = 1, Knod
               ijk = i + jk
               nrows = nrows + 1
               rowptr(nrows) = nnz + 1
               DO p = 1, Knod
                 pm1s = (p-1) * Knsq
                 DO m = 1, Knod
                   mm1 = (m-1) * Knod
                    mp = mm1 + pm1s
                   DO l = 1, Knod
                     lmp = l + mp
                     nnz = nnz + 1
                     values(nnz) = LaplacianCenter(lmp,ijk,el)
                     colidx(nnz) = colctr + lmp
                   ENDDO
                 ENDDO
               ENDDO
               DO ijP = 0, FC
                 IF (elemID(i2f(ijP),el) > 0) THEN
                   DO p = 1, Knod
                     pm1s = (p-1) * Knsq
                     DO m = 1, Knod
                       mm1 = (m-1) * Knod
                        mp = mm1 + pm1s
                       DO l = 1, Knod
                         lmp = l + mp
                         nnz = nnz + 1
                         values(nnz) = LaplacianNeigbr(lmp,ijk,ijP,el)
                         colidx(nnz) = colnbr(ijP) + lmp
                       ENDDO
                     ENDDO
                   ENDDO
                 ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO
           
       ! Filter out zero's
!       print *, 'Zero Count: ', count( values(1:nnz) == 0.0d0 ), nnz, rowptr(nrows+1)
!       print *, 'Zero Count: ', count( ABS(values(1:nnz)) < 1.0d-14 ), nnz, rowptr(nrows+1)
       print *, 'Zero Count: ', nnz-count( ABS(values(1:nnz)) < 1.0d-14 ), rowptr(nrows+1)-1
!       nnz = nnz - count( values(1:nnz) == 0.0d0 )
       nnz = nnz - count( ABS(values(1:nnz)) < 1.0d-14 )

print *,'MATRIX MAX VALUE ',maxval(abs(values(1:nnz)))

       allocate( values_filtered(nnz), colidx_filtered(nnz), rowptr_filtered(nrows+1) )
       nnz = 0
       do i = 1, nrows
          rowptr_filtered(i) = nnz+1
          do l = rowptr(i), rowptr(i+1)-1
!             if (values(l) .ne. 0.0d0 ) then
             if (ABS(values(l)) .ge. 1.0d-14 ) then
                nnz = nnz + 1
                values_filtered(nnz) = values(l)
                colidx_filtered(nnz) = colidx(l)
             endif
          enddo
       enddo
       rowptr_filtered(nrows+1) = nnz+1
       rowptr_filtered = rowptr_filtered - 1
       colidx_filtered = colidx_filtered - 1

! Convert Fortran 1-index to C/C++ 0-index
       rowptr = rowptr - 1
       colidx = colidx - 1

       ierr = APLLES_Initialize()

       ! Use the filtered or raw matrix?
       if (.false.) then
         ierr = APLLES_Setup_Matrix_CSR(nrows,rowptr,colidx,values,A_handle)
       else
         ierr = APLLES_Setup_Matrix_CSR(nrows,rowptr_filtered,colidx_filtered,values_filtered,A_handle)
       endif

       ierr = APLLES_Matrix_Copy_CSR_To_Other(A_handle, "CSR", A_handle)

       ierr = APLLES_Setup_Solver(A_handle,trim(Aplles_solver_name),S_handle)

       ierr = APLLES_Setup_Precon(A_handle,trim(Aplles_precon_name),P_handle)

       deAllocate (rowptr, colidx, values)
       deAllocate (LaplacianCenter,LaplacianNeigbr)

       deAllocate (rowptr_filtered, colidx_filtered, values_filtered)

END SUBROUTINE SetLaplacian


SUBROUTINE GetLaplacian(HuynhSolver_type, Vortin, psi, A_handle, P_handle, S_handle)
       USE params
       USE variables
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer HuynhSolver_type
       real*8, dimension(Knod, Knod, Knod, Nel) :: Vortin, psi
       real*8, dimension(Knod**NDIM, Nel) :: temp

       integer el, bel, i, j, k, l, m, ij, lm, jm1, mm1, km1s, jk, ijk, jx, jy, lx, ly
       integer colc, cole, colw, coln, cols, nnz, nrows, Knsq, K4, ierr
       real*8 residual, volume

       real*8, dimension(:), allocatable :: b, x

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

  
       nrows = Knod**NDIM * Nel
       Allocate (x(nrows), b(nrows))

       Knsq = Knod * Knod

       temp = 0.d0
       DO bel = 1, NelB
         el = boundaryPointsElemID(bel)
         DO k = 1, Knod
           km1s = (k-1) * Knsq
           DO j = 1, Knod
             jk = (j-1) * Knod + km1s
             DO i = 1, Knod
               ijk = i + jk
               DO m = 1, Knod
                 temp(ijk,el) = temp(ijk,el) + DOT_PRODUCT( BndrySrc(1:Knod,m,ijk,bel) , BC_Values(1:Knod,m,bel) )
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       residual = 0.d0
       volume = 0.d0
       DO el = 1, Nel
         DO k = 1, Knod
!           km1s = (k-1) * Knsq
           DO j = 1, Knod
!             jm1 = (j-1) * Knod
             DO i = 1, Knod
!               ij = jm1 + i
               residual = residual + wgt(i) * wgt(j) * wgt(k) * Vortin(i,j,k,el) * Vol_Jac(i,j,k,el)
               volume = volume + wgt(i) * wgt(j) * wgt(k) * Vol_Jac(i,j,k,el)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       print *,'Res, Vol, eps ', residual, volume, (1+residual)/volume

       x = 0.d0
       b = 0.d0

       nrows = 0
       DO el = 1, Nel
         DO k = 1, Knod
           km1s = (k-1) * Knsq
           DO j = 1, Knod
             jk = (j-1) * Knod + km1s
             DO i = 1, Knod
               ijk = i + jk
               nrows = nrows + 1
!             b(nrows) = - Vortin(i,j,el) - BndrySrc(ij,el)
!             b(nrows) = - Vol_Jac(i,j,el) * Vortin(i,j,el) - BndrySrc(ij,el)
               if(abs((1-residual)/volume) > .1) then
                b(nrows) = - Vol_Jac(i,j,k,el) * Vortin(i,j,k,el) - temp(ijk,el)
               else
                b(nrows) = - Vol_Jac(i,j,k,el) * Vortin(i,j,k,el) - temp(ijk,el)
!              b(nrows) = - Vol_Jac(i,j,el) * (-(1+residual)/volume + Vortin(i,j,el)) - temp(ij,el)
               endif
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       ierr = APLLES_Solve(A_handle, x, b, S_handle, P_handle)

       nrows = 0
       DO el = 1, Nel
         DO k = 1, Knod
           DO j = 1, Knod
             DO i = 1, Knod
               nrows = nrows + 1
               psi(i,j,k,el) = x(nrows)
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       deAllocate (x, b)

END SUBROUTINE GetLaplacian


SUBROUTINE GetLaplacGrads(HuynhSolver_type, Phi, getCurl)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type, getCurl
       real*8, dimension(Knod, Knod, Knod, Nel) :: Phi

       integer i, j, k, jx, jy, jz, el, eln, ijP, ijPm, idir, ibnd
       real*8, dimension(Knod, NDIM) :: loclPhi
       ! Phi and GradPhi at the L/R (0:1), S/N (2:3), and B/T (4:5) boundaries of a mesh
       real*8, dimension(Knod, Knod, 0:FC, Nel) :: bndrPhi, bndrGradPhi
       ! Common Phi and GradPhi at the L/R (0:1), S/N (2:3), and B/T (4:5) boundaries of a mesh
       real*8, dimension(Knod, Knod, 0:FC, Nel) :: comPhi, comGradPhi

       real*8, dimension(NDIM, NDIM) :: Dxsi_iDx_j
       real*8, dimension(Knod, Knod) :: Acoef, Bcoef, Ccoef
       real*8, dimension(Knod, Knod, 0:FC) :: FaceA, FaceB, FaceC, NormA
       real*8, dimension(Knod, Knod, Knod, NDIM) :: SolnA, SolnB, SolnC

       real*8, dimension(Knod, Knod, Knod) :: f_tilda, g_tilda, h_tilda
       real*8, dimension(Knod, Knod, 0:1) :: f_tildaB, g_tildaB, h_tildaB

       real*8 du_dxsi, du_deta, du_dzta

       real*8 tmpx, tmpy, tmpz, Ovr_Jac

   INTERFACE

     PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

     END FUNCTION Cross_Product

   END INTERFACE

       comPhi = 0.d0
       comGradPhi = 0.d0

       ! Extrapolate the unknown, uin, and its derivative to the mesh boundaries using Lagrange polynomials of order Knod-1
       DO el = 1, Nel

         DO k = 1, Knod
           DO j = 1, Knod

            !!! Local 1D tensor notation in the x (=1), y (=2), and z (=3) directions
             !! IMPORTANT: the indices (i,j) have been flipped for the y direction because they will be used
             !!            more efficiently in the dot-product operation (also allows for streamlining)
             !!            the same approach is applied to the z direction
             DO i = 1, Knod
               ! Phi and Velocity in local coordinates (along horizontal (1:Knod,j) solution points)
               loclPhi(i,1) =  Phi(i,j,k,el)
               ! Phi and Velocity in local coordinates (along vertical (j,1:Knod) solution points)
               loclPhi(i,2) =  Phi(j,i,k,el)
               ! Phi and Velocity in local coordinates (along horizontal (1:Knod,j) solution points)
               loclPhi(i,3) =  Phi(j,k,i,el)
             ENDDO

             ijP = 0
             ! Extrapolation operations in x (=1) and y (=2) directions
             DO idir = 1, NDIM
               ! Extraploated boundary values of Phi and GradPhi to left/south/bottom (=0) and right/north/top (=1)
               DO ibnd = 0, 1
                 bndrPhi(j,k,ijP,el) = dot_product(loclPhi(1:Knod,idir), SolnBndryLgrangeBasis(1:Knod,ibnd))
                 bndrGradPhi(j,k,ijP,el) = dot_product(loclPhi(1:Knod,idir), SolnBndryGradLgrangeBasis(1:Knod,ibnd))
                 ijP = ijP + 1
               ENDDO
             ENDDO

           ENDDO
         ENDDO

       ENDDO

       IF (HuynhSolver_type .eq. 2) THEN

!         gLprime(0) = -0.5d0*Knod*Knod
!         gLprime(1) =  0.5d0*Knod*knod
         gLprime(0) = -0.5d0*Knod*(Knod+1)
         gLprime(1) =  0.5d0*Knod*(Knod+1)

         ! Get the common values of Phi at the mesh interfaces (we use simple averaging in this case)
         ! We definitely need a more efficient strategy - right now we're saving com-s twice, and can probably save on bndrGradPhi as well
         ! We probably can and should write a more compact method so that we don't repeat the same stuff for NEWS (which becomes messier for NEWSBT in 3D)
         DO el = 1, Nel

           ijP = 0
           ! Extrapolation operations in x (=1), y (=2), and z (=3) directions
           DO idir = 1, NDIM
             DO ibnd = 0, 1
               ! mesh to the left/south of left/south face or right/north of right/north face
               eln = elemID(i2f(ijP),el)
               IF (eln .lt. 0) THEN

                 CALL GetLapBndryComVals(-eln, ibnd, Face_Acoef(1:Knod,1:Knod,ijP,el), Face_Bcoef(1:Knod,1:Knod,ijP,el),  &
                                          Face_Ccoef(1:Knod,1:Knod,ijP,el), Face_Norm(1:Knod,1:Knod,ijP,el),              &
                                          Face_Jac(1:Knod,1:Knod,ijP,el),                                                 &
                                          bndrPhi(1:Knod,1:Knod,ijP,el), bndrGradPhi(1:Knod,1:Knod,ijP,el),               &
                                          comPhi(1:Knod,1:Knod,ijP,el), comGradPhi(1:Knod,1:Knod,ijP,el))

                 IF (BC_Switch(-eln) .eq. DirichletBC) THEN

                   VelocJump(1:Knod,1:Knod,-eln) = comGradPhi(1:Knod,1:Knod,ijP,el) / Face_Norm(1:Knod,1:Knod,ijP,el)

                 ELSEIF (BC_Switch(-eln) .eq. NeumannBC) THEN

!!                   VelocJump(1:Knod,-eln) = comPhi(1:Knod,ijP,el)
!!!  NEED TO FIX THIS
!                   DO j = 1, Knod
!                     VelocJump(j,-eln) = dot_product(comPhi(1:Knod,ijP,el), SolnNodesGradLgrangeBasis(1:Knod,j)) / &
!                                                     Face_Norm(j,ijP,el)
!                     IF (getCurl .eq. 1 .and. idir .eq. 2) VelocJump(j,-eln) = - VelocJump(j,-eln)
!                   ENDDO

                 ENDIF

               ELSEIF (eln .gt. 0) THEN

                 ijPm = nbr(ijP)
                 CALL GetLapIntrnComVals_Meth2(ibnd, Face_Acoef(1:Knod,1:Knod,ijP,el), Face_Bcoef(1:Knod,1:Knod,ijP,el),  &
                                               Face_Ccoef(1:Knod,1:Knod,ijP,el), Face_Jac(1:Knod,1:Knod,ijP,el),          &
                                               bndrPhi(1:Knod,1:Knod,ijP,el), bndrPhi(1:Knod,1:Knod,ijPm,eln),            &
                                               bndrGradPhi(1:Knod,1:Knod,ijP,el), comPhi(1:Knod,1:Knod,ijP,el))

               ENDIF
               ijP = ijP + 1
             ENDDO
           ENDDO

         ENDDO

         DO el = 1, Nel

           DO ijP = 0, FC
             eln = elemID(i2f(ijP),el)
             ! common GradPhi(el,right) is Average of bndrGradPhi(el,right) + (el+1,left); same works for Left, South, and North
             IF (eln .gt. 0) comGradPhi(1:Knod,1:Knod,ijP,el) = &
                                0.5d0*(bndrGradPhi(1:Knod,1:Knod,ijP,el) + bndrGradPhi(1:Knod,1:Knod,nbr(ijP),eln))
           ENDDO

         ENDDO

       ELSEIF (HuynhSolver_type .eq. 11) THEN

         gLprime(0) = -0.5d0*Knod*(Knod+1)
         gLprime(1) =  0.5d0*Knod*(Knod+1)

         DO el = 1, Nel

           ijP = 0
           ! Extrapolation operations in x (=1) and y (=2) directions
           DO idir = 1, NDIM
             DO ibnd = 0, 1
               ! mesh to the left/south of left/south face or right/north of right/north face
               eln = elemID(i2f(ijP),el)
               IF (eln .lt. 0) THEN

!                 CALL GetLapBndryComVals(-eln, ibnd, Face_Acoef(1:Knod,ijP,el), Face_Bcoef(1:Knod,ijP,el), &
!                                          Face_Norm(1:Knod,ijP,el), Face_Jac(1:Knod,ijP,el), &
!                                          bndrPhi(1:Knod,ijP,el), bndrGradPhi(1:Knod,ijP,el), &
!                                          comPhi(1:Knod,ijP,el), comGradPhi(1:Knod,ijP,el))

                 IF (BC_Switch(-eln) .eq. DirichletBC) THEN

!                   VelocJump(1:Knod,-eln) = comGradPhi(1:Knod,ijP,el) / Face_Norm(1:Knod,ijP,el)

                 ELSEIF (BC_Switch(-eln) .eq. NeumannBC) THEN
!!                   VelocJump(1:Knod,-eln) = comPhi(1:Knod,ijP,el)
!                   DO j = 1, Knod
!                     VelocJump(j,-eln) = dot_product(comPhi(1:Knod,ijP,el), SolnNodesGradLgrangeBasis(1:Knod,j)) / &
!                                                     Face_Norm(j,ijP,el)
!                     IF (getCurl .eq. 1 .and. idir .eq. 2) VelocJump(j,-eln) = - VelocJump(j,-eln)
!                   ENDDO

                 ENDIF

               ELSEIF (eln .gt. 0) THEN

                 IF (ibnd .eq. 1) THEN

                   ijPm = nbr(ijP)
!                   CALL GetLapIntrnComVals_Meth11(ibnd, Face_Acoef(1:Knod,ijP,el), Face_Bcoef(1:Knod,ijP,el), &
!                                                  Face_Jac(1:Knod,ijP,el), bndrPhi(1:Knod,ijP,el), bndrGradPhi(1:Knod,ijP,el), &
!                                                  comPhi(1:Knod,ijP,el), comGradPhi(1:Knod,ijP,el), &
!                                                  Face_Acoef(1:Knod,ijPm,eln), Face_Bcoef(1:Knod,ijPm,eln), &
!                                                  Face_Jac(1:Knod,ijPm,eln), bndrPhi(1:Knod,ijPm,eln), &
!                                                  bndrGradPhi(1:Knod,ijPm,eln), comPhi(1:Knod,ijPm,eln), &
!                                                  comGradPhi(1:Knod,ijPm,eln))

                 ENDIF

               ENDIF
               ijP = ijP + 1

             ENDDO
           ENDDO

         ENDDO

       ENDIF


       DO el = 1, Nel

         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
               Ovr_Jac = 1.d0 / Vol_Jac(jx,jy,jz,el)

               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el))
               ! Jac * dzta / d_j
               Dxsi_iDx_j(:,3) = Cross_Product(Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el))

               ! Axx, Axy, Axz in notes
               SolnA(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnB(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,2)) * Ovr_Jac
               SolnC(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,3)) * Ovr_Jac

               ! Ayy in notes
               ! Saving in cyclic order makes it possible to lump operations in x, y, and z directions in one uniform call
! *** May actually have to do jy,jx,jz for better performance (looks that way)
               SolnA(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,2)) * Ovr_Jac
               SolnB(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnC(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,3)) * Ovr_Jac

               ! Azz in notes
               SolnA(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,3)) * Ovr_Jac
               SolnB(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnC(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,2)) * Ovr_Jac
             ENDDO
           ENDDO
         ENDDO

         ! Get first derivatives of the variable
         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod

               ! Get grads of unknown uin along xsi (for eta=const) and along eta (for xsi=const
               du_dxsi = dot_product(Phi(1:Knod,jy,jz,el), SolnNodesGradLgrangeBasis(1:Knod,jx)) + &
                         dot_product(comPhi(jy,jz,0:1,el) - bndrPhi(jy,jz,0:1,el), NodesGradRadau(jx,0:1))
               du_deta = dot_product(Phi(jx,1:Knod,jz,el), SolnNodesGradLgrangeBasis(1:Knod,jy)) + &
                         dot_product(comPhi(jx,jz,2:3,el) - bndrPhi(jx,jz,2:3,el), NodesGradRadau(jy,0:1))
               du_dzta = dot_product(Phi(jx,jy,1:Knod,el), SolnNodesGradLgrangeBasis(1:Knod,jz)) + &
                         dot_product(comPhi(jx,jy,4:5,el) - bndrPhi(jx,jy,4:5,el), NodesGradRadau(jz,0:1))

               ! Get f~ and g~ as per Huynh's paper
               f_tilda(jx,jy,jz) = du_dxsi * SolnA(jx,jy,jz,1) + du_deta * SolnB(jx,jy,jz,1) + du_dzta * SolnC(jx,jy,jz,1)
               g_tilda(jy,jx,jz) = du_deta * SolnA(jy,jx,jz,2) + du_dxsi * SolnB(jy,jx,jz,2) + du_dzta * SolnC(jy,jx,jz,2)
               h_tilda(jz,jx,jy) = du_dzta * SolnA(jz,jx,jy,3) + du_dxsi * SolnB(jz,jx,jy,3) + du_deta * SolnC(jz,jx,jy,3)

             ENDDO
           ENDDO
         ENDDO

         ! Get ?_tildas at the mesh boundaries
         DO ibnd = 0, 1
           DO k = 1, Knod
             DO j = 1, Knod
               f_tildaB(j,k,ibnd) = comGradPhi(j,k,ibnd  ,el) - dot_product(f_tilda(1:Knod,j,k), SolnBndryLgrangeBasis(1:Knod,ibnd))      ! Left of mesh - derivative
               g_tildaB(j,k,ibnd) = comGradPhi(j,k,ibnd+2,el) - dot_product(g_tilda(1:Knod,j,k), SolnBndryLgrangeBasis(1:Knod,ibnd))      ! South of mesh - derivative
               h_tildaB(j,k,ibnd) = comGradPhi(j,k,ibnd+4,el) - dot_product(h_tilda(1:Knod,j,k), SolnBndryLgrangeBasis(1:Knod,ibnd))      ! Bottom of mesh - derivative
             ENDDO
           ENDDO
         ENDDO

         ! Get potential velocity
         IF (getCurl .eq. 0) THEN
           DO jz = 1, Knod
             DO jy = 1, Knod
               DO jx = 1, Knod
                 tmpx = ( f_tilda(jx,jy,jz) + dot_product(f_tildaB(jy,jz,0:1), NodesRadau(jx,0:1)) ) / Dxsi_iDx_j(1,1)
                 tmpy = ( g_tilda(jy,jx,jz) + dot_product(g_tildaB(jx,jz,0:1), NodesRadau(jy,0:1)) ) / Dxsi_iDx_j(2,2)
                 tmpz = ( h_tilda(jz,jx,jy) + dot_product(h_tildaB(jx,jy,0:1), NodesRadau(jz,0:1)) ) / Dxsi_iDx_j(3,3)
                 Uvel(jx,jy,jz,el) = tmpx
                 Vvel(jx,jy,jz,el) = tmpy
                 Wvel(jx,jy,jz,el) = tmpz
!                 Uvel(jx,jy,jz,el) = tmpx * Dxsi_iDx_j(1,1) * Ovr_Jac
!                 Uvel(jx,jy,el) = Uvel(jx,jy,el) + ( tmpx * Vol_Dx_iDxsi_j(jx,jy,1,1,el) + tmpy * Vol_Dx_iDxsi_j(jx,jy,1,2,el) ) / &
!                                                     Vol_Jac(jx,jy,el)
!                 Vvel(jx,jy,el) = Vvel(jx,jy,el) + ( tmpx * Vol_Dx_iDxsi_j(jx,jy,2,1,el) + tmpy * Vol_Dx_iDxsi_j(jx,jy,2,2,el) ) / &
!                                                     Vol_Jac(jx,jy,el)
               ENDDO
             ENDDO
           ENDDO
         ! Get vortical velocity
         ELSEIF (getCurl .eq. 1) THEN
           DO jz = 1, Knod
             DO jy = 1, Knod
               DO jx = 1, Knod
                 tmpx = ( f_tilda(jx,jy,jz) + dot_product(f_tildaB(jy,jz,0:1), NodesRadau(jx,0:1)) ) / Dxsi_iDx_j(1,1)
                 tmpy = ( g_tilda(jy,jx,jz) + dot_product(g_tildaB(jx,jz,0:1), NodesRadau(jy,0:1)) ) / Dxsi_iDx_j(2,2)
                 tmpz = ( h_tilda(jz,jx,jy) + dot_product(h_tildaB(jx,jy,0:1), NodesRadau(jz,0:1)) ) / Dxsi_iDx_j(3,3)
                 Uvel(jx,jy,jz,el) =  tmpy
                 Vvel(jx,jy,jz,el) = -tmpx
                 Wvel(jx,jy,jz,el) = 0.d0
!               tmpx = f_tilda(jx,jy) + dot_product(f_tildaB(jy,0:1), NodesRadau(jx,0:1))
!!               tmpy = g_tilda(jx,jy) + dot_product(g_tildaB(jx,0:1), NodesRadau(jy,0:1))
!               tmpy = g_tilda(jy,jx) + dot_product(g_tildaB(jx,0:1), NodesRadau(jy,0:1))
!               Uvel(jx,jy,el) =  ( tmpx * Vol_Dx_iDxsi_j(jx,jy,2,1,el) + tmpy * Vol_Dx_iDxsi_j(jx,jy,2,2,el) ) / &
!                                                   Vol_Jac(jx,jy,el)
!               Vvel(jx,jy,el) = -( tmpx * Vol_Dx_iDxsi_j(jx,jy,1,1,el) + tmpy * Vol_Dx_iDxsi_j(jx,jy,1,2,el) ) / &
!                                                   Vol_Jac(jx,jy,el)
               ENDDO
             ENDDO
           ENDDO
         ELSE
           print *,'ERROR: getCurl is either 0 (for potential velocity) or 1 (for vortical velocity)!"'
           stop
         ENDIF
       ENDDO

END SUBROUTINE GetLaplacGrads


SUBROUTINE GetLapBndryComVals(eln, ijPos, Acoef, Bcoef, Ccoef, Anorm, Jac, bndrPhi, bndrDPhi, comPhi, comDPhi)
       USE params
       USE variables

       implicit NONE
       integer eln, ijPos
       ! dimension of calling routine for the following is (1:Knod,ijPos,el)
       ! boundary metric stuff; extrapolated boundary values (per mesh) of Phi and its grad, as well as their common values (com)
       real*8, dimension(Knod, Knod) :: Acoef, BCoef, Ccoef, Anorm, Jac, bndrPhi, bndrDPhi, comPhi, comDPhi

       integer j, k, l, km1, lk
       ! cross derivatives using common values comPhi
       real*8, dimension(Knod, Knod) :: crossDPhi_j, crossDPhi_k
       ! small dense matrix to obtain comPhi (per element) for a given Neumann BC
       real*8 NeuMatrix(Knod**2, Knod**2), NeuRHS(Knod**2)


       IF (BC_Switch(eln) .eq. DirichletBC) THEN

         comPhi(1:Knod,1:Knod) = BC_Values(1:Knod,1:Knod,eln)

         ! Get the corrected values of grad(Phi) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
         bndrDPhi(1:Knod,1:Knod) = bndrDPhi(1:Knod,1:Knod) + (comPhi(1:Knod,1:Knod)-bndrPhi(1:Knod,1:Knod))*gLprime(ijPos) 

         ! Get the corrected values of grad(Phi) along the mesh interface;
         DO k = 1, Knod
           DO j = 1, Knod
             crossDPhi_j(j,k) = dot_product(comPhi(1:Knod,k), SolnNodesGradLgrangeBasis(1:Knod,j))
             crossDPhi_k(j,k) = dot_product(comPhi(j,1:Knod), SolnNodesGradLgrangeBasis(1:Knod,k))
           ENDDO
         ENDDO

         ! Get the interface f  (dimension(Knod,Knod))
         comDPhi = (Acoef*bndrDPhi + Bcoef*crossDPhi_j + Ccoef*crossDPhi_k) / Jac

       ELSEIF (BC_Switch(eln) .eq. NeumannBC) THEN

         comDPhi(1:Knod,1:Knod) = BC_Values(1:Knod,1:Knod,eln) * Anorm(1:Knod,1:Knod)

!         NeuRHS(1:Knod) = Jac(1:Knod)*comDPhi(1:Knod) - Acoef(1:Knod)* (bndrDPhi(1:Knod) - bndrPhi(1:Knod)*gLprime(ijPos))
!         DO j = 1, Knod
!           NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
!           NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime(ijPos)
!         ENDDO


                     NeuMatrix = 0.d0
                 DO k = 1, Knod
                   km1 = (k-1) * Knod
                   DO l = 1, Knod
                     lk = l + km1
!                    DO j = 1, Knod
!                      jl = j + lm1
!                       jk = j + km1
!                       NeuMatrix_Orig(jk,lk) = NeuMatrix_Orig(jk,lk) + FaceB(j,k,ijP) * SolnNodesGradLgrangeBasis(l,j)
!                       NeuMatrix_Orig(jk,jl) = NeuMatrix_Orig(jk,jl) + FaceC(j,k,ijP) * SolnNodesGradLgrangeBasis(l,k)
!                     ENDDO
         NeuRHS(lk) = Jac(l,k)*comDPhi(l,k) - Acoef(l,k)* (bndrDPhi(l,k) - bndrPhi(l,k)*gLprime(ijPos))
                     NeuMatrix(lk,lk) = NeuMatrix(lk,lk) + Acoef(l,k)*gLprime(ijPos)
                   ENDDO
                 ENDDO
         CALL Gauss_Solver(Knod**2, NeuMatrix, comPhi, NeuRHS)

       ENDIF

END SUBROUTINE GetLapBndryComVals


SUBROUTINE GetLapIntrnComVals_Meth2(ijPos, Acoef, Bcoef, Ccoef, Jac, bndrPhi, bndrNbrPhi, bndrDPhi, comPhi)
       USE params
       USE variables

       implicit NONE
       integer ijPos
       ! dimension of calling routine is (1:Knod,ijPos,el)
       ! boundary metric stuff, extrapolated boundary values (per mesh) of Phi and its grad, and common interface value (com)
       real*8, dimension(Knod, Knod) :: Acoef, BCoef, Ccoef, Jac, bndrPhi, bndrNbrPhi, bndrDPhi, comPhi

       integer j, k
       ! cross derivatives using common values comPhi
       real*8, dimension(Knod, Knod) :: crossDPhi_j, crossDPhi_k


       ! ComPhi: Average of PhiBndry(el,right) + PhiBndy(el+1,left)
       comPhi = 0.5d0*(bndrPhi + bndrNbrPhi)

       ! Get the corrected values of grad(Phi) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
       bndrDPhi(1:Knod,1:Knod) = bndrDPhi(1:Knod,1:Knod) + (comPhi(1:Knod,1:Knod)-bndrPhi(1:Knod,1:Knod))*gLprime(ijPos) 

       ! Get the corrected values of grad(Phi) along the mesh interface;
       DO k = 1, Knod
         DO j = 1, Knod
           crossDPhi_j(j,k) = dot_product(comPhi(1:Knod,k), SolnNodesGradLgrangeBasis(1:Knod,j))
           crossDPhi_k(j,k) = dot_product(comPhi(j,1:Knod), SolnNodesGradLgrangeBasis(1:Knod,k))
         ENDDO
       ENDDO

       ! Get the interface f
       bndrDPhi = (Acoef*bndrDPhi + Bcoef*crossDPhi_j + Ccoef*crossDPhi_k) / Jac

END SUBROUTINE GetLapIntrnComVals_Meth2


SUBROUTINE Gauss_Solver(n,A,x,b)

! Gauss Elimination solver without pivoting (for now) copied from 
       implicit NONE

       integer n
       real*8 A(n,n), x(n), b(n)

       integer i, j, k
       real*8 c

! Decomposition (Elimination)
       DO k = 1, n-1
         DO i = k+1, n
           c = a(i,k) / a(k,k)
           a(i,k) = 0.d0
           b(i) = b(i) - c*b(k)
           DO j = k+1, n
             a(i,j) = a(i,j) - c*a(k,j)
           ENDDO
         ENDDO
       ENDDO

! Back-substitution
       x(n) = b(n)/a(n,n)
       DO i = n-1, 1, -1
         c = 0.d0
         DO j = i+1, n
           c = c + a(i,j)*x(j)
         ENDDO 
         x(i) = (b(i)- c) / a(i,i)
       ENDDO 

END SUBROUTINE Gauss_Solver


FUNCTION inv(A) RESULT(Ainv)

       integer dp
       parameter (dp=8)
       real(dp), dimension(:,:), intent(in) :: A
       real(dp), dimension(size(A,1),size(A,2)) :: Ainv

       real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
       integer, dimension(size(A,1)) :: ipiv   ! pivot indices
       integer :: n, info

       ! External procedures defined in LAPACK
       external DGETRF
       external DGETRI

       ! Store A in Ainv to prevent it from being overwritten by LAPACK
       Ainv = A
       n = size(A,1)

       ! DGETRF computes an LU factorization of a general M-by-N matrix A
       ! using partial pivoting with row interchanges.
       call DGETRF(n, n, Ainv, n, ipiv, info)

       if (info /= 0) then
         stop 'Matrix is numerically singular!'
       end if

       ! DGETRI computes the inverse of a matrix using the LU factorization
       ! computed by DGETRF.
       call DGETRI(n, Ainv, n, ipiv, work, n, info)

       if (info /= 0) then
         stop 'Matrix inversion failed!'
       end if

END FUNCTION inv


SUBROUTINE SolveConvectionDiffusion(HuynhSolver_type, tIntegrator_type, Reyn, dt, numStep, dumpFreq, prob_type)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer HuynhSolver_type, tIntegrator_type, numStep, dumpFreq, prob_type
       real*8 Reyn, dt

       integer nt, i, j, k, el, kx, ky, kz, lx, ly, lz
       real*8 x, y, z, tmp
       real*8 dto2, dto3, dto6
       real*8, dimension(Knod, Knod, Knod, Nel) :: Vort_tmp, Vort_str, f_of_Vort, psi, k1, k2, k3, k4
       real*8 t1, t2, t3, t4, t5, tlap, tcrl, tdif, tcon

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

  INTERFACE

    FUNCTION itoa(i) RESULT(res)

       character(11) res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
    END FUNCTION itoa

  END INTERFACE

       Vort = Vort0             ! leaving Vort0 alone, in case needed elsewhere; generally, we don't need to save in real problems
       Uvel = 0.d0
       Vvel = 0.d0
       Wvel = 0.d0
       VelocJump = 0.d0

       dto2 = dt / 2.d0   ! of course this is lazy programming :)
       dto3 = dt / 3.d0
       dto6 = dt / 6.d0

       CALL SetLaplacian(HuynhSolver_type, A_handle, P_handle, S_handle)

 tlap = 0.d0
 tcrl = 0.d0
 tdif = 0.d0
 tcon = 0.d0

       DO nt = 1, numStep
        IF (MOD(nt,1) .eq. 0) print *,'timestep ',nt
         Vort_tmp = Vort

         DO el = 1, Nel
         DO k = 1, Knod
         DO j = 1, Knod
         DO i = 1, Knod
!           IF (Abs(Vort(i,j,k,el)) < 1.d-10) Vort(i,j,k,el) = 0.d0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         IF (tIntegrator_type .eq. 1) THEN

           CALL  EulerTimeIntegrate(dt, Reyn, Vort, k1, HuynhSolver_type, A_handle, P_handle, S_handle)
           Vort = Vort + k1

         ELSEIF (tIntegrator_type .eq. 2) THEN

           CALL  EulerTimeIntegrate(dt, Reyn, Vort, k1, HuynhSolver_type, A_handle, P_handle, S_handle)
           CALL  EulerTimeIntegrate(dt, Reyn, Vort + k1, k2, HuynhSolver_type, A_handle, P_handle, S_handle)
           Vort = Vort + 0.5d0 * (k1 + k2)

if (.false.) then
           Vort_str = Vort_tmp + dto2*f_of_Vort
         t1 = omp_get_wtime()
!           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 2
         t2 = omp_get_wtime()
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
!           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
         t3 = omp_get_wtime()
     print *,'stage 2 ',MaxVal(Sqrt(Uvel**2 + Vvel**2 + Wvel**2))
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
         t4 = omp_get_wtime()
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
         t5 = omp_get_wtime()
 tlap = tlap + t2 - t1
 tcrl = tcrl + t3 - t2
 tdif = tdif + t4 - t3
 tcon = tcon + t5 - t4
           Vort = Vort + dt*f_of_Vort
endif
         ELSEIF (tIntegrator_type .eq. 4) THEN

           ! I have no idea where I copied this RK4 from; seems incorrect
!           CALL  EulerTimeIntegrate(dt, Reyn, Vort, k1, HuynhSolver_type, A_handle, P_handle, S_handle)
!           CALL  EulerTimeIntegrate(dto2, Reyn, Vort + k1/2.d0, k2, HuynhSolver_type, A_handle, P_handle, S_handle)
!           CALL  EulerTimeIntegrate(dto2, Reyn, Vort + k2, k3, HuynhSolver_type, A_handle, P_handle, S_handle)
!           CALL  EulerTimeIntegrate(dt, Reyn, Vort + 2.d0*k3, k4, HuynhSolver_type, A_handle, P_handle, S_handle)
!           Vort = Vort + (k1 + 4.d0*k2 + 4.d0*k3 + k4) / 6.d0

           ! We don't need to save all k's; we can definitely improve below or use a low-storage SSP method
           CALL  EulerTimeIntegrate(dt, Reyn, Vort, k1, HuynhSolver_type, A_handle, P_handle, S_handle)
           CALL  EulerTimeIntegrate(dt, Reyn, Vort + k1/2.d0, k2, HuynhSolver_type, A_handle, P_handle, S_handle)
           CALL  EulerTimeIntegrate(dt, Reyn, Vort + k2/2.d0, k3, HuynhSolver_type, A_handle, P_handle, S_handle)
           CALL  EulerTimeIntegrate(dt, Reyn, Vort + k3, k4, HuynhSolver_type, A_handle, P_handle, S_handle)
           Vort = Vort + (k1 + 2.d0*k2 + 2.d0*k3 + k4) / 6.d0

if (.false.) then
           Vort = Vort + dto6*f_of_Vort
           Vort_str = Vort_tmp + dto2*f_of_Vort
!           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 2
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
!           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)

           Vort = Vort + dto3*f_of_Vort
           Vort_str = Vort_tmp + dto2*f_of_Vort
!           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 3
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
!           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)

           Vort = Vort + dto3*f_of_Vort
           Vort_str = Vort_tmp + dt*f_of_Vort
!           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 4
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
!           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)

           Vort = Vort + dto6*f_of_Vort
 endif

         ENDIF

!         IF (MOD(nt,dumpFreq) .eq. 0) CALL dumpResult(nt, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)
!         IF (MOD(nt,dumpFreq) .eq. 0) CALL Paraview_Dump(HuynhSolver_type, nt, nt*dt)
         IF (MOD(nt,dumpFreq) .eq. 0) CALL Paraview_Output(nt)

       ENDDO

!print *,'avg tlap ',0.5*tlap/numStep
!print *,'avg tcrl ',0.5*tcrl/numStep
!print *,'avg tdif ',0.5*tdif/numStep
!print *,'avg tcon ',0.5*tcon/numStep

!       IF (MOD(numStep,dumpFreq) .ne. 0) &
!           CALL dumpResult(numStep, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)

!       IF (MOD(numStep,dumpFreq) .ne. 0) CALL Paraview_Dump(HuynhSolver_type, numStep, numStep*dt)
       IF (MOD(numStep,dumpFreq) .ne. 0) CALL Paraview_Output(numStep)

CONTAINS

  SUBROUTINE EulerTimeIntegrate(dt, Reyn, VortIn, VortOut, HuynhSolver_type, A_handle, P_handle, S_handle)
       USE params
       USE variables
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer HuynhSolver_type
       real*8 Reyn, dt
       real*8, dimension(Knod, Knod, Knod, Nel) :: VortIn, VortOut, f_of_Vort, psi

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

       integer i, j

print *,'**** CONVECTION **** '

         BC_Values = BC_Psi
         BC_Switch = BC_Switch_Lap
         CALL GetLaplacian(HuynhSolver_type, VortIn, psi, A_handle, P_handle, S_handle)                                ! Stage 1
         CALL GetLaplacGrads(HuynhSolver_type, psi, 1)

         BC_Values = BC_VelNorm
         CALL GetConvectedFlux(HuynhSolver_type, VortIn, f_of_Vort)

         VortOut =  dt*f_of_Vort

print *,'Post Conv Max Vort ',minval(vortin+vortOut),maxval(vortin+vortOut)

print *,'**** DIFFUSION **** '

 if (.false.) then
         VortIn = VortOut
         BC_Values = BC_Psi
         CALL GetLaplacian(HuynhSolver_type, VortIn, psi, A_handle, P_handle, S_handle)                                ! Stage 1
         CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
 endif

print *,'Velocity ',minval(sqrt(uvel**2 + vvel**2)),maxval(sqrt(uvel**2 + vvel**2))

print *,'Jump ',minval(VelocJump),maxval(VelocJump)

         DO i = 1, NBndry
           IF (NoSlip(i)) THEN
             DO j = BndryNum(i-1)+1, BndryNum(i)
!             DO j = 1,nelb
               VelocJump(1:Knod,1:Knod,j) = VelocJump(1:Knod,1:Knod,j) - BC_VelParl(1:Knod,1:Knod,j)
               BC_Values(1:Knod,1:Knod,j) = VelocJump(1:Knod,1:Knod,j) * Reyn / dt
             ENDDO
           ELSE
             DO j = BndryNum(i-1)+1, BndryNum(i)
               BC_Values(1:Knod,1:Knod,j) = 0.d0
             ENDDO
           ENDIF
         ENDDO

         BC_Switch = NeumannBC
         CALL GetDiffusedFlux(HuynhSolver_type, VortIn, f_of_Vort)

         VortOut = VortOut + f_of_Vort * dt / Reyn

!print *,'Post Diff Max Vort ',minval(vortin+vortOut),maxval(vortin+vortOut)

         RETURN

  END SUBROUTINE EulerTimeIntegrate

END SUBROUTINE SolveConvectionDiffusion


SUBROUTINE GetDiffusedFlux(HuynhSolver_type, Phi, LapPhi)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type
       real*8, dimension(Knod, Knod, Knod, Nel) :: Phi, LapPhi

       integer i, j, k, jx, jy, jz, el, eln, ijP, ijPm, idir, ibnd
       real*8, dimension(Knod, NDIM) :: loclPhi
       ! Phi and GradPhi at the L/R (0:1), S/N (2:3), and B/T (4:5) boundaries of a mesh
       real*8, dimension(Knod, Knod, 0:FC, Nel) :: bndrPhi, bndrGradPhi
       ! Common Phi and GradPhi at the L/R (0:1), S/N (2:3), and B/T (4:5) boundaries of a mesh
       real*8, dimension(Knod, Knod, 0:FC, Nel) :: comPhi, comGradPhi

       real*8, dimension(NDIM, NDIM) :: Dxsi_iDx_j
       real*8, dimension(Knod, Knod) :: Acoef, Bcoef, Ccoef
       real*8, dimension(Knod, Knod, 0:FC) :: FaceA, FaceB, FaceC, NormA
       real*8, dimension(Knod, Knod, Knod, NDIM) :: SolnA, SolnB, SolnC

       real*8, dimension(Knod, Knod, Knod) :: f_tilda, g_tilda, h_tilda
       real*8, dimension(Knod, Knod, 0:1) :: f_tildaB, g_tildaB, h_tildaB

       real*8 du_dxsi, du_deta, du_dzta, Ovr_Jac

   INTERFACE

     PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

     END FUNCTION Cross_Product

   END INTERFACE

       comPhi = 0.d0
       comGradPhi = 0.d0

       ! Extrapolate the unknown, uin, and its derivative to the mesh boundaries using Lagrange polynomials of order Knod-1
       DO el = 1, Nel

         DO k = 1, Knod
           DO j = 1, Knod

            !!! Local 1D tensor notation in the x (=1), y (=2), and z (=3) directions
             !! IMPORTANT: the indices (i,j) have been flipped for the y direction because they will be used
             !!            more efficiently in the dot-product operation (also allows for streamlining)
             !!            the same approach is applied to the z direction
             DO i = 1, Knod
               ! Phi and Velocity in local coordinates (along horizontal (1:Knod,j) solution points)
               loclPhi(i,1) =  Phi(i,j,k,el)
               ! Phi and Velocity in local coordinates (along vertical (j,1:Knod) solution points)
               loclPhi(i,2) =  Phi(j,i,k,el)
               ! Phi and Velocity in local coordinates (along horizontal (1:Knod,j) solution points)
               loclPhi(i,3) =  Phi(j,k,i,el)
             ENDDO

             ijP = 0
             ! Extrapolation operations in x (=1) and y (=2) directions
             DO idir = 1, NDIM
               ! Extraploated boundary values of Phi and GradPhi to left/south/bottom (=0) and right/north/top (=1)
               DO ibnd = 0, 1
                 bndrPhi(j,k,ijP,el) = dot_product(loclPhi(1:Knod,idir), SolnBndryLgrangeBasis(1:Knod,ibnd))
                 bndrGradPhi(j,k,ijP,el) = dot_product(loclPhi(1:Knod,idir), SolnBndryGradLgrangeBasis(1:Knod,ibnd))
                 ijP = ijP + 1
               ENDDO
             ENDDO

           ENDDO
         ENDDO

       ENDDO

       IF (HuynhSolver_type .eq. 2) THEN

!         gLprime(0) = -0.5d0*Knod*Knod
!         gLprime(1) =  0.5d0*Knod*knod
         gLprime(0) = -0.5d0*Knod*(Knod+1)
         gLprime(1) =  0.5d0*Knod*(Knod+1)

         ! Get the common values of Phi at the mesh interfaces (we use simple averaging in this case)
         ! We definitely need a more efficient strategy - right now we're saving com-s twice, and can probably save on bndrGradPhi as well
         ! We probably can and should write a more compact method so that we don't repeat the same stuff for NEWS (which becomes messier for NEWSBT in 3D)
         DO el = 1, Nel

           ijP = 0
           ! Extrapolation operations in x (=1) and y (=2) directions
           DO idir = 1, NDIM
             DO ibnd = 0, 1
               ! mesh to the left/south of left/south face or right/north of right/north face
               eln = elemID(i2f(ijP),el)
               IF (eln .lt. 0) THEN

                 CALL GetLapBndryComVals(-eln, ibnd, Face_Acoef(1:Knod,1:Knod,ijP,el), Face_Bcoef(1:Knod,1:Knod,ijP,el),  &
                                          Face_Ccoef(1:Knod,1:Knod,ijP,el), Face_Norm(1:Knod,1:Knod,ijP,el),              &
                                          Face_Jac(1:Knod,1:Knod,ijP,el),                                                 &
                                          bndrPhi(1:Knod,1:Knod,ijP,el), bndrGradPhi(1:Knod,1:Knod,ijP,el),               &
                                          comPhi(1:Knod,1:Knod,ijP,el), comGradPhi(1:Knod,1:Knod,ijP,el))

                 VelocJump(1:Knod,1:Knod,-eln) = comPhi(1:Knod,1:Knod,ijP,el)

               ELSEIF (eln .gt. 0) THEN

                 ijPm = nbr(ijP)
                 CALL GetLapIntrnComVals_Meth2(ibnd, Face_Acoef(1:Knod,1:Knod,ijP,el), Face_Bcoef(1:Knod,1:Knod,ijP,el),  &
                                               Face_Ccoef(1:Knod,1:Knod,ijP,el), Face_Jac(1:Knod,1:Knod,ijP,el),          &
                                               bndrPhi(1:Knod,1:Knod,ijP,el), bndrPhi(1:Knod,1:Knod,ijPm,eln),            &
                                               bndrGradPhi(1:Knod,1:Knod,ijP,el), comPhi(1:Knod,1:Knod,ijP,el))

               ENDIF
               ijP = ijP + 1
             ENDDO
           ENDDO

         ENDDO

         DO el = 1, Nel

           DO ijP = 0, FC
             eln = elemID(i2f(ijP),el)
             ! common GradPhi(el,right) is Average of bndrGradPhi(el,right) + (el+1,left); same works for Left, South, and North
             IF (eln .gt. 0) comGradPhi(1:Knod,1:Knod,ijP,el) = &
                                0.5d0*(bndrGradPhi(1:Knod,1:Knod,ijP,el) + bndrGradPhi(1:Knod,1:Knod,nbr(ijP),eln))
           ENDDO

         ENDDO

       ELSEIF (HuynhSolver_type .eq. 11) THEN

         gLprime(0) = -0.5d0*Knod*(Knod+1)
         gLprime(1) =  0.5d0*Knod*(Knod+1)

         DO el = 1, Nel

           ijP = 0
           ! Extrapolation operations in x (=1) and y (=2) directions
           DO idir = 1, NDIM
             DO ibnd = 0, 1
               ! mesh to the left/south of left/south face or right/north of right/north face
               eln = elemID(i2f(ijP),el)
               IF (eln .lt. 0) THEN

!                 CALL GetLapBndryComVals(-eln, ibnd, Face_Acoef(1:Knod,ijP,el), Face_Bcoef(1:Knod,ijP,el), &
!                                          Face_Norm(1:Knod,ijP,el), Face_Jac(1:Knod,ijP,el), &
!                                          bndrPhi(1:Knod,ijP,el), bndrGradPhi(1:Knod,ijP,el), &
!                                          comPhi(1:Knod,ijP,el), comGradPhi(1:Knod,ijP,el))

!                 VelocJump(1:Knod,-eln) = comPhi(1:Knod,ijP,el)

               ELSEIF (eln .gt. 0) THEN

                 IF (ibnd .eq. 1) THEN

                   ijPm = nbr(ijP)
!                   CALL GetLapIntrnComVals_Meth11(ibnd, Face_Acoef(1:Knod,ijP,el), Face_Bcoef(1:Knod,ijP,el), &
!                                                  Face_Jac(1:Knod,ijP,el), bndrPhi(1:Knod,ijP,el), bndrGradPhi(1:Knod,ijP,el), &
!                                                  comPhi(1:Knod,ijP,el), comGradPhi(1:Knod,ijP,el), &
!                                                  Face_Acoef(1:Knod,ijPm,eln), Face_Bcoef(1:Knod,ijPm,eln), &
!                                                  Face_Jac(1:Knod,ijPm,eln), bndrPhi(1:Knod,ijPm,eln), &
!                                                  bndrGradPhi(1:Knod,ijPm,eln), comPhi(1:Knod,ijPm,eln), &
!                                                  comGradPhi(1:Knod,ijPm,eln))

                 ENDIF

               ENDIF
               ijP = ijP + 1

             ENDDO
           ENDDO

         ENDDO

       ENDIF


       DO el = 1, Nel

         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
               Ovr_Jac = 1.d0 / Vol_Jac(jx,jy,jz,el)

               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el))
               ! Jac * dzta / d_j
               Dxsi_iDx_j(:,3) = Cross_Product(Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el))

               ! Axx, Axy, Axz in notes
               SolnA(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnB(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,2)) * Ovr_Jac
               SolnC(jx,jy,jz,1) = Dot_Product(Dxsi_iDx_j(:,1),Dxsi_iDx_j(:,3)) * Ovr_Jac

               ! Ayy in notes
               ! Saving in cyclic order makes it possible to lump operations in x, y, and z directions in one uniform call
! *** May actually have to do jy,jx,jz for better performance (looks that way)
               SolnA(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,2)) * Ovr_Jac
               SolnB(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnC(jy,jx,jz,2) = Dot_Product(Dxsi_iDx_j(:,2),Dxsi_iDx_j(:,3)) * Ovr_Jac

               ! Azz in notes
               SolnA(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,3)) * Ovr_Jac
               SolnB(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,1)) * Ovr_Jac
               SolnC(jz,jx,jy,3) = Dot_Product(Dxsi_iDx_j(:,3),Dxsi_iDx_j(:,2)) * Ovr_Jac
             ENDDO
           ENDDO
         ENDDO

         ! Get first derivatives of the variable
         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod

               ! Get grads of unknown uin along xsi (for eta=const) and along eta (for xsi=const
               du_dxsi = dot_product(Phi(1:Knod,jy,jz,el), SolnNodesGradLgrangeBasis(1:Knod,jx)) + &
                         dot_product(comPhi(jy,jz,0:1,el) - bndrPhi(jy,jz,0:1,el), NodesGradRadau(jx,0:1))
               du_deta = dot_product(Phi(jx,1:Knod,jz,el), SolnNodesGradLgrangeBasis(1:Knod,jy)) + &
                         dot_product(comPhi(jx,jz,2:3,el) - bndrPhi(jx,jz,2:3,el), NodesGradRadau(jy,0:1))
               du_dzta = dot_product(Phi(jx,jy,1:Knod,el), SolnNodesGradLgrangeBasis(1:Knod,jz)) + &
                         dot_product(comPhi(jx,jy,4:5,el) - bndrPhi(jx,jy,4:5,el), NodesGradRadau(jz,0:1))

               ! Get f~ and g~ as per Huynh's paper
               f_tilda(jx,jy,jz) = du_dxsi * SolnA(jx,jy,jz,1) + du_deta * SolnB(jx,jy,jz,1) + du_dzta * SolnC(jx,jy,jz,1)
               g_tilda(jy,jx,jz) = du_deta * SolnA(jy,jx,jz,2) + du_dxsi * SolnB(jy,jx,jz,2) + du_dzta * SolnC(jy,jx,jz,2)
               h_tilda(jz,jx,jy) = du_dzta * SolnA(jz,jx,jy,3) + du_dxsi * SolnB(jz,jx,jy,3) + du_deta * SolnC(jz,jx,jy,3)

             ENDDO
           ENDDO
         ENDDO

         ! Get ?_tildas at the mesh boundaries
         DO ibnd = 0, 1
           DO k = 1, Knod
             DO j = 1, Knod
               ! Left of mesh - derivative
               f_tildaB(j,k,ibnd) = comGradPhi(j,k,ibnd  ,el) - dot_product(f_tilda(1:Knod,j,k), SolnBndryLgrangeBasis(1:Knod,ibnd))
               ! South of mesh - derivative
               g_tildaB(j,k,ibnd) = comGradPhi(j,k,ibnd+2,el) - dot_product(g_tilda(1:Knod,j,k), SolnBndryLgrangeBasis(1:Knod,ibnd))
               ! Bottom of mesh - derivative
               h_tildaB(j,k,ibnd) = comGradPhi(j,k,ibnd+4,el) - dot_product(h_tilda(1:Knod,j,k), SolnBndryLgrangeBasis(1:Knod,ibnd))
             ENDDO
           ENDDO
         ENDDO

         ! Now get the Laplacian
         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
               LapPhi(jx,jy,jz,el) = (dot_product(f_tilda(1:Knod,jy,jz), SolnNodesGradLgrangeBasis(1:Knod,jx))   + &
                                      dot_product(f_tildaB(jy,jz,0:1), NodesGradRadau(jx,0:1))                   + &
                                      dot_product(g_tilda(1:Knod,jx,jz), SolnNodesGradLgrangeBasis(1:Knod,jy))   + &
                                      dot_product(g_tildaB(jx,jz,0:1), NodesGradRadau(jy,0:1))                   + &
                                      dot_product(h_tilda(1:Knod,jx,jy), SolnNodesGradLgrangeBasis(1:Knod,jz))   + &
                                      dot_product(h_tildaB(jx,jy,0:1), NodesGradRadau(jz,0:1)))  / Vol_Jac(jx,jy,jz,el)
             ENDDO
           ENDDO
         ENDDO

       ENDDO
  
END SUBROUTINE GetDiffusedFlux


SUBROUTINE GetConvectedFlux(HuynhSolver_type, Phi, dPhi)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type
       real*8, dimension(Knod, Knod, Knod, Nel) :: Phi, dPhi

       integer, parameter, dimension(0:FC) :: sgn = (/-1, 1, -1, 1, -1, 1/)
       integer i, j, k, jx, jy, jz, el, eln, ijP, ijPm, idir, ibnd
       real*8 bndrVel
       real*8, dimension(Knod, NDIM) :: loclPhi, loclVel
       real*8, dimension(Knod, Knod, Knod, NDIM, Nel) :: discFlx            ! nodal values of Discontinous Fluxes in X and Y dirs
       ! Phi, interface flux, and flux jumps at the L/R (0:1) and S/N (2:3) boundaries of a mesh
       real*8, dimension(Knod, Knod, 0:FC, Nel) :: bndrPhi, bndrFlx, upwndFlx

       real*8, dimension(NDIM, NDIM) :: Dxsi_iDx_j

       real*8, dimension(Knod, Knod, 0:FC) :: bndrDiscFlx

   INTERFACE

     PURE FUNCTION Cross_Product(x,y)
       implicit NONE

       real*8, dimension(3), intent(in) :: x, y
       real*8, dimension(3) :: Cross_Product

     END FUNCTION Cross_Product

   END INTERFACE

       bndrPhi = 0.d0
       bndrFlx = 0.d0
       discFlx = 0.d0
       upwndFlx = 0.d0

       ! Extrapolate the unknown, Phi and the Flux to the mesh boundaries using Lagrange polynomials of order Knod-1
       DO el = 1, Nel

         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
               ! Jac * dxsi / dx_j
               Dxsi_iDx_j(:,1) = Cross_Product(Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el))
               ! Jac * deta / dx_j
               Dxsi_iDx_j(:,2) = Cross_Product(Vol_Dx_iDxsi_j(:,3,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el))
               ! Jac * dzta / d_j
               Dxsi_iDx_j(:,3) = Cross_Product(Vol_Dx_iDxsi_j(:,1,jx,jy,jz,el),Vol_Dx_iDxsi_j(:,2,jx,jy,jz,el))
             ENDDO
           ENDDO
         ENDDO

         DO k = 1, Knod
           DO j = 1, Knod

           !!! Local 1D tensor notation in the x (=1) and y (=2) directions
             !! IMPORTANT: the indices (i,j) have been flipped for the y direction because they will be used
             !!            more efficiently in the dot-product operation (also allows for streamlining)
             DO i = 1, Knod
               ! Phi and Velocity in local coordinates (along horizontal (1:Knod,j) solution points)
               loclPhi(i,1) =  Phi(i,j,k,el)
               loclVel(i,1) =  Uvel(i,j,k,el)*Dxsi_iDx_j(1,1) + Vvel(i,j,k,el)*Dxsi_iDx_j(2,1) + &
                               Wvel(i,j,k,el)*Dxsi_iDx_j(3,1)
               ! Phi and Velocity in local coordinates (along vertical (j,1:Knod) solution points)
               loclPhi(i,2) =  Phi(j,i,k,el)
               loclVel(i,2) =  Uvel(j,i,k,el)*Dxsi_iDx_j(1,2) + Vvel(j,i,k,el)*Dxsi_iDx_j(2,2) + &
                               Wvel(j,i,k,el)*Dxsi_iDx_j(3,2)
               ! Phi and Velocity in local coordinates (along horizontal (1:Knod,j) solution points)
               loclPhi(i,3) =  Phi(j,k,i,el)
               loclVel(i,3) =  Uvel(j,k,i,el)*Dxsi_iDx_j(1,3) + Vvel(j,k,i,el)*Dxsi_iDx_j(2,3) + &
                               Wvel(j,k,i,el)*Dxsi_iDx_j(3,3)
             ENDDO

             ijP = 0
             ! Extrapolation operations in x (=1) and y (=2) directions
             DO idir = 1, NDIM

             ! Extraploated boundary values of Phi and Velocity*Phi to left/south (=0) and right/north (=1)
               DO ibnd = 0, 1
                 bndrPhi(j,k,ijP,el) = dot_product(loclPhi(1:Knod,idir), SolnBndryLgrangeBasis(1:Knod,ibnd))

                 ! mesh to the left of left face or right of right face
                 eln = elemID(i2f(ijP),el)
                 IF (eln .lt. 0) THEN
!           IF (BC_Switch(-eln) .eq. DirichletBC) THEN
                   bndrFlx(j,k,ijP,el) = bndrPhi(j,k,ijP,el) * BC_Values(j,k,-eln)
!           ELSEIF (BC_Switch(-eln) .eq. NeumannBC) THEN
!           ENDIF
                 ELSE
                   bndrFlx(j,k,ijP,el) = bndrPhi(j,k,ijP,el) * dot_product(loclVel(1:Knod,idir), &
                                                                           SolnBndryLgrangeBasis(1:Knod,ibnd))
                 ENDIF
                 ijP = ijP + 1
               ENDDO

               ! Discontinuous flux values on (1:Knod,j) solution points for horizontal direction
               ! and on (j,1:Knod) solution points for vertical direction
               DO i = 1, Knod
                 discFlx(i,j,k,idir,el) = loclVel(i,idir) * loclPhi(i,idir)
               ENDDO

             ENDDO

           ENDDO
         ENDDO

       ENDDO


       DO el = 1, Nel

         DO ijP = 0, FC
           ijPm = nbr(ijP)
           ! mesh to the left/south of left/south face or right/north of right/north face
           eln = elemID(i2f(ijP),el)
           IF (eln .gt. 0) THEN
             DO k = 1, Knod
               DO j = 1, Knod
                 IF (Abs(bndrPhi(j,k,ijP,el) - bndrPhi(j,k,ijPm,eln)) .gt. 1.0d-6) THEN
                   bndrVel = (bndrFlx(j,k,ijP,el) - bndrFlx(j,k,ijPm,eln)) / (bndrPhi(j,k,ijP,el) - bndrPhi(j,k,ijPm,eln))
                   upwndFlx(j,k,ijP,el) = 0.5d0*( bndrFlx(j,k,ijP,el) + bndrFlx(j,k,ijPm,eln) + sgn(ijP) * Abs(bndrVel) * &
                                                ( bndrPhi(j,k,ijP,el) - bndrPhi(j,k,ijPm,eln) ) )
                 ELSE
                   upwndFlx(j,k,ijP,el) = 0.5d0*( bndrFlx(j,k,ijP,el) + bndrFlx(j,k,ijPm,eln) )
                 ENDIF
               ENDDO
             ENDDO
           ELSEIF (eln .lt. 0) THEN
           ! Dirichlet BC
!           IF (BC_Switch(-eln) .eq. DirichletBC) THEN
             DO k = 1, Knod
               DO j = 1, Knod
                 upwndFlx(j,k,ijP,el) = bndrFlx(j,k,ijP,el)
               ENDDO
             ENDDO
           ! Neumann BC
!           ELSEIF (BC_Switch(-eln) .eq. NeumannBC) THEN
!           ENDIF
           ENDIF
         ENDDO

       ENDDO


       DO el = 1, Nel

         ijP = 0
         DO idir = 1, NDIM
           DO ibnd = 0, 1
             DO k = 1, Knod
               DO j = 1, Knod
                 bndrDiscFlx(j,k,ijP) = upwndFlx(j,k,ijP,el) - &
                                        dot_product(discFlx(1:Knod,j,k,idir,el), SolnBndryLgrangeBasis(1:Knod,ibnd))
               ENDDO
             ENDDO
             ijP = ijP + 1
           ENDDO
         ENDDO

         ! Now get the convection
         DO jz = 1, Knod
           DO jy = 1, Knod
             DO jx = 1, Knod
!             dPhi(jx,jy,el) = dPhi(jx,jy,el) -  &
               dPhi(jx,jy,jz,el) = - &
                               (dot_product(discFlx(1:Knod,jy,jz,1,el), SolnNodesGradLgrangeBasis(1:Knod,jx))   + &
                                dot_product(bndrDiscFlx(jy,jz,0:1), NodesGradRadau(jx,0:1))                     + &
                                dot_product(discFlx(1:Knod,jx,jz,2,el), SolnNodesGradLgrangeBasis(1:Knod,jy))   + &
                                dot_product(bndrDiscFlx(jx,jz,2:3), NodesGradRadau(jy,0:1))                     + &
                                dot_product(discFlx(1:Knod,jx,jy,3,el), SolnNodesGradLgrangeBasis(1:Knod,jz))   + &
                                dot_product(bndrDiscFlx(jx,jy,4:5), NodesGradRadau(jz,0:1))) / Vol_Jac(jx,jy,jz,el)
             ENDDO
           ENDDO
         ENDDO

       ENDDO
  
END SUBROUTINE GetConvectedFlux
