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
       integer, parameter :: DefaultBC = 0, DirichletBC = 1, NeumannBC = 2

       integer Knod, Lnod, Nel, NelB

       logical exact
END MODULE params


MODULE variables
       implicit NONE
                                              
       real*8, allocatable :: xcoord(:), ycoord(:)  ! nodal positions of the geometry elements (structured elements for now)
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
       real*8, allocatable :: VelocJump(:,:)
       real*8, allocatable :: gps(:)                ! scaled/parametric coordinates of the geometric nodes (in 1-D)
       real*8, allocatable :: sps(:)                ! scaled/parametric coordinates of the solution points (in 1-D)
       real*8, allocatable :: wgt(:)                ! Gaussian weights at the solution points (in 1-D)

       real*8, allocatable :: SJacs(:,:)
       real*8, allocatable :: Vort0(:,:,:), Vort(:,:,:)   ! initial Vort0(i,kx,ky) and convected+diffused Vort(i,kx,ky) values of the unknowns
                                                    ! per element i and at (kx,ky) nodes per elem
       real*8, allocatable :: Uvel(:,:,:), Vvel(:,:,:)    ! u(i,kx,ky) and v(i,kx,ky) velocity components
                                                    ! per element i and at (kx,ky) nodes per elem
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
       ! Neumann or Dirichlet BC. The patches are numbered from the bottom/south and go counter-clockwise
       ! S=1; E=2; N=3; W=4
       ! For now, memory is over-allocated based on the original value of NelB (NelB gets modified properly in MeshSetup)
       real*8, allocatable :: BC_VelNorm(:,:), BC_VelParl(:,:)
       real*8, allocatable ::  BC_Values(:,:)                  ! BC values for the Knod points of NelB boundary elements
       integer, allocatable ::  BC_Switch(:)                   ! The BC switch for these patches is Dirichlet = 1; Neumann = 2

       real*8, allocatable :: BndrySrc(:,:,:)                  ! Poisson Solver's RHS term to be dotted by BC_Values

       character(len=30) :: aplles_solver_name, aplles_precon_name

END MODULE variables


MODULE CnvrtTensor2FemIndices

       ! We use tensor product of 1-D arrays, which is incosistent with FEM index notation
       ! This index conversion is the easiest way to handle this, especially when dealing
       ! w/ higher order elements.  BUT, we have to manually create the indices for each case

       ! index t2f(i,j) converts node (i,j) to traditional CCW FEM nodes, or other FEM formats (to be included)

       USE params

       PRIVATE

       INTERFACE T2F
         MODULE PROCEDURE T2F1D, T2F2D
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
         print *,'Only up to quadratic finite elements (Lnod = 3) are supported!'
         STOP
       ENDIF
  END FUNCTION T2F1D

  FUNCTION T2F2D(i,j)
       implicit NONE

       integer t2f2d
       integer i,j

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

END MODULE CnvrtTensor2FemIndices


PROGRAM TwoD_Vorticity_Transport
       USE params
       USE variables
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer Nelx, Nely, Lnod_in, prob_type, HuynhSolver_type, tIntegrator_type, numStep, dumpFreq, fast
       real*8 Reyn, fac, dt
       character(132) dum

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle
 
       !!!! User Input Data
       open(unit=2, file='input.dat', status='old')
         read(2,'(a)') dum
         read(2,*) Nelx, Nely           ! num. of meshes/elements in x, y dirs (structured grids; must generalize to unstructured)
         read(2,'(a)') dum
         read(2,*) Knod                 ! solution nodes per elem (must generalize to include x, y dirs)
         read(2,'(a)') dum
         read(2,*) Lnod_in              ! geometry parametric order per elem (must generalize to include x, y dirs)
         read(2,'(a)') dum
         read(2,'(L)') exact            ! Exact geometry for concentric cylinders
         read(2,'(a)') dum
         read(2,*) Reyn                 ! Reynolds number
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
         IF (prob_type .eq. 2 .or. prob_type .eq. 6 .or. prob_type .eq. 5 .or. prob_type > 8) THEN
           print *,'problem type unavailable '
           stop
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
       ! Assuming a Square Domain
       Nel = Nelx * Nely
       NelB = 2 * (Nelx + Nely)

       ! Nothing to do here: Just allocate arrays
       CALL Allocate_Arrays(Nelx, Nely, Nel, Nelb, Knod, Lnod)

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
       CALL SetupMesh(Nelx, Nely, fac, prob_type)

       ! Nothing to do here: UNLESS we want to add new problem types
       !                     Set up IC, BC, and Source for a given problem
       ! NOTE: Make sure to add the corresponding exact solution for the problem type in Function u_xct
       CALL SetupICBCandSrc(Nel, NelB, Knod, Lnod, prob_type)

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

END PROGRAM TwoD_Vorticity_Transport


SUBROUTINE Allocate_Arrays(Nx, Ny, N, NB, K, L)
       USE variables

       implicit NONE
       integer Nx, Ny, N, NB, K, L
       integer Lm1, Ksq

       Lm1 = L - 1
       Ksq = K*K
 
       allocate (xcoord((Nx*Lm1+1)*(Ny*Lm1+1)), ycoord((Nx*Lm1+1)*(Ny*Lm1+1)))
       allocate (nodeID(L*L,N))
       allocate (elemID(4,N))
       allocate (BoundaryNodeID(L,NB), BC_Values(K,NB), BC_VelNorm(K,NB), BC_VelParl(K,NB))
       allocate (BC_Switch(NB))
       allocate (boundaryPointsElemID(NB))
       allocate (VelocJump(k,NB))
       allocate (BndrySrc(K,Ksq,NB))
       allocate (sps(K), wgt(K))
       allocate (gps(L))
       allocate (SJacs(Ksq,N))
       allocate (Vort0(K,K,N), Vort(K,K,N), Uvel(K,K,N), Vvel(K,K,N))
       allocate (SolnBndryLgrangeBasis(K,0:1), SolnBndryGradLgrangeBasis(K,0:1))
       allocate (NodesRadau(K,0:1), NodesGradRadau(K,0:1))
       allocate (SolnNodesGradLgrangeBasis(K,K))
       allocate (GeomBndryLgrangeBasis(L,0:1), GeomBndryGradLgrangeBasis(L,0:1))
       allocate (GeomNodesLgrangeBasis(L,K), GeomNodesGradLgrangeBasis(L,K))

END SUBROUTINE Allocate_Arrays


SUBROUTINE deAllocate_Arrays
       USE variables

       deallocate (xcoord, ycoord)
       deallocate (nodeID)
       deallocate (elemID)
       deallocate (BoundaryNodeID, BC_Values, BC_VelNorm, BC_VelParl)
       deallocate (BC_Switch)
       deallocate (boundaryPointsElemID)
       deallocate (VelocJump)
       deallocate (BndrySrc)
       deallocate (sps, wgt)
       deallocate (gps)
       deallocate (SJacs)
       deallocate (Vort0, Vort, Uvel, Vvel)
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


SUBROUTINE SetupMesh(Nelx, Nely, fac, prob_type)
       USE params
       USE variables

       implicit NONE
       integer Nelx, Nely, prob_type
       real*8 fac

       ! A few grid types
       ! square box w/ random perturbation = 1
       ! square box w/ sinusoidal perturbation = 2
       ! a pair of concentric circles = 3
       ! square boxes 1 and 2 for periodic BC
       integer grid_type

       integer i, j, Lm1, nb, nl, nr, ns, nn, nc, nd, nd_sav
       real*8 xL, yL, dx, dy, rnd_num(8*2*Nel)
       real*8 R_in, R_out, theta_in, theta_out, dr, dtheta, R_j, theta_i, percent
       logical periodic, left_boundary, right_boundary, south_boundary, north_boundary

       periodic = .FALSE.
       CALL RANDOM_NUMBER(rnd_num)

       xL = 2.249d0*2*pi
       yL = 3.d0
       xL = 2.25d0
       yL = 2*2.25d0

       xL = 1.d0
       yL = 1.d0

       grid_type = 2
       IF (prob_type .le. 2) periodic = .TRUE.
       IF (prob_type .gt. 6) grid_type = 3

       ! Quick Paraview viz of geometry (useful when checking effect of fac randomization of grids)
       OPEN(unit = 8, file = 'Geometry.vtk', status = 'unknown')
       write(8,'(a)') '# vtk DataFile Version 3.0'
       write(8,'(a)') '2D Unstructured Grid of Quads'
       write(8,'(a)') 'ASCII'
       write(8,'(a)') ' '
       write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'

       ! poor man's attempt at creating unstructured mesh using simple structured data; just for testing
       IF (grid_type .eq. 1) THEN

         IF (Lnod .gt. 2) THEN
           print *,'Only bilinear elements allowed for grid_type 1! '
           STOP
         ENDIF

         Lm1 = Lnod - 1
         dx = xL / (Nelx * Lm1)
         dy = yL / (Nely * Lm1)

         ! Set geometry nodes (not the same as solution nodes)
         !j = 1; i = 1
         nd = 1
         xcoord(nd) = 0.d0
         ycoord(nd) = 0.d0

         DO i = 2, Nelx
           nd = nd + 1
           xcoord(nd) = (i-1) * dx
           ycoord(nd) = 0.d0
         ENDDO

         !j = 1; i = Nelx + 1
         nd = nd + 1
         xcoord(nd) = xL
         ycoord(nd) = 0.d0

         DO j = 2, Nely
           nd = nd + 1
           xcoord(nd) = 0.d0
           ycoord(nd) = (j-1) * dy
           DO i = 2, Nelx
             nd = nd + 1
             xcoord(nd) = (i-1) * dx
             ycoord(nd) = (j-1) * dy
             IF (fac .lt. 1.d0) THEN
               xcoord(nd) = xcoord(nd) + fac*(0.5-rnd_num(nd))*dx
               ycoord(nd) = ycoord(nd) + fac*(0.5-rnd_num(Nel+nd))*dy
             ENDIF
           ENDDO
           nd = nd + 1
           xcoord(nd) = xL
           ycoord(nd) = (j-1) * dy
         ENDDO

         !j = Nely+1; i = 1
         nd = nd + 1
         xcoord(nd) = 0.d0
         ycoord(nd) = yL

         DO i = 2, Nelx
           nd = nd + 1
           xcoord(nd) = (i-1) * dx
           ycoord(nd) = yL
         ENDDO

         !j = Nely+1; i = Nelx+1
         nd = nd + 1
         xcoord(nd) = xL
         ycoord(nd) = yL

       ELSEIF (grid_type .eq. 2) THEN

         Lm1 = Lnod - 1
         dx = xL / (Nelx * Lm1)
         dy = yL / (Nely * Lm1)

         nd = 0
         DO j = 1, Lm1*Nely
           DO i = 1, Lm1*Nelx
             nd = nd + 1
             xcoord(nd) = (i-1)*dx + sin((i-1)*2*pi*dx) * sin((j-1)*2*pi*dy) * fac
             ycoord(nd) = (j-1)*dy + sin((i-1)*2*pi*dx) * sin((j-1)*2*pi*dy) * fac
           ENDDO
           nd = nd + 1
           xcoord(nd) = xL
           ycoord(nd) = (j-1)*dy
         ENDDO

         DO i = 1, Lm1*Nelx
           nd = nd + 1
           xcoord(nd) = (i-1)*dx
           ycoord(nd) = yL ! (j-1)*dy
         ENDDO
         nd = nd + 1
         xcoord(nd) = xL
         ycoord(nd) = yL

       ELSEIF (grid_type .eq. 3) THEN

         Lm1 = Lnod - 1
         R_in = 0.5d0
         R_out = 1.d0
         dr = (R_out - R_in) / (Nely * Lm1)
         print *,'enter percent rotation '
         read *, percent
         theta_in = 0.d0 + percent*2.d0*pi
         theta_out = 2.d0*pi *(1.0d0 + percent)
         print *,theta_out - theta_in,2.d0*pi
         dtheta = (theta_out - theta_in) / (Nelx * Lm1)

         nd = 0
         DO j = 1, Lm1*Nely+1
           DO i = 1, Lm1*Nelx
             nd = nd + 1
             theta_i = theta_in + (i-1)*dtheta
             R_j = R_in + (j-1)*dr
             if (j .eq. Lm1*Nely+1) R_j = R_out
             xcoord(nd) = R_j * cos(-theta_i)
             ycoord(nd) = R_j * sin(-theta_i)
           ENDDO
         ENDDO

       ELSE
         print *,'selected grid_type is not supported! '
         STOP
       ENDIF

       ! Paraview viz
       write(8,'(a,i8,a)') 'POINTS ', nd, ' float'
       DO i = 1, nd
         write(8,*) real(xcoord(i)),real(ycoord(i)),' 0.0'
       ENDDO
       nd_sav = nd

       ! Set mesh connectivity to nodes (nodeID); also get connectivity to neighboring meshes (elemID)
       ! nodeID starts from SW corners and goes CCW; elemID starts with S face and goes CCW
       ! NOTE: This is of course a poor man's (structured grid) approach to the general problem
       !       Also, the -ve numbers refer to boundary faces, which are assigned in the BC_IC setup routine
       nb = 0
       nc = 0
       nd = 0
       DO j = 1, Nely
         DO i = 1, Nelx
           nc = nc + 1
           nd = nd + 1

           left_boundary = .false.
           right_boundary = .false.
           south_boundary = .false.
           north_boundary = .false.

           IF (j .eq. 1) THEN
!             IF (periodic) THEN
!               elemID(1,nc) = nc + Nelx * (Nely - 1)
!             ELSE
               south_boundary = .true.
               nb = nb + 1
               ns = nb
               elemID(1,nc) = -nb
               boundaryPointsElemID(nb) = nc
!             ENDIF
           ELSE
             elemID(1,nc) = nc - Nelx
           ENDIF
           IF (j .eq. Nely) THEN
!             IF (periodic) THEN
!               elemID(3,nc) = nc - Nelx * (Nely - 1)
!             ELSE
               north_boundary = .true.
               nb = nb + 1
               nn = nb
               elemID(3,nc) = -nb
               boundaryPointsElemID(nb) = nc
!             ENDIF
           ELSE
             elemID(3,nc) = nc + Nelx
           ENDIF
           ! for concentric circles there is no boundary in the i direction
           ! also true of periodic BC
           IF (i .eq. 1) THEN
             IF (grid_type .eq. 3 .or. periodic) THEN
!             IF (grid_type .eq. 3) THEN
                elemID(4,nc) = nc + Nelx - 1
             ELSE
               left_boundary = .true.
               nb = nb + 1
               nl = nb
               elemID(4,nc) = -nb
               boundaryPointsElemID(nb) = nc
             ENDIF
           ELSE             
             elemID(4,nc) = nc - 1
           ENDIF
           IF (i .eq. Nelx) THEN
             IF (grid_type .eq. 3 .or. periodic) THEN
!             IF (grid_type .eq. 3) THEN
                elemID(2,nc) = nc + 1 - Nelx
             ELSE
               right_boundary = .true.
               nb = nb + 1
               nr = nb
               elemID(2,nc) = -nb
               boundaryPointsElemID(nb) = nc
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

           IF (left_boundary) THEN
             ! CCW rotation for meshing
!             BoundaryNodeID(2,nl) = nodeID(1,nc)
!             BoundaryNodeID(1,nl) = nodeID(4,nc)
             ! For now, I'm assuming all sides are in positive x/y direction 
             ! NOTE: need to account for direction later
             BoundaryNodeID(2,nl) = nodeID(4,nc)
             BoundaryNodeID(1,nl) = nodeID(1,nc)
           ENDIF

           IF (right_boundary) THEN
             BoundaryNodeID(1,nr) = nodeID(2,nc)
             BoundaryNodeID(2,nr) = nodeID(3,nc)
           ENDIF

           IF (south_boundary) THEN
             BoundaryNodeID(1,ns) = nodeID(1,nc)
             BoundaryNodeID(2,ns) = nodeID(2,nc)
           ENDIF

           IF (north_boundary) THEN
             ! CCW rotation for meshing
!             BoundaryNodeID(1,nn) = nodeID(3,nc)
!             BoundaryNodeID(2,nn) = nodeID(4,nc)
             ! For now, I'm assuming all sides are in positive x/y direction 
             ! NOTE: need to account for direction later
             BoundaryNodeID(1,nn) = nodeID(4,nc)
             BoundaryNodeID(2,nn) = nodeID(3,nc)
           ENDIF

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
         ENDDO

         ! accounting for the node Nelx+1, except for grid_type 3
         IF (grid_type .ne. 3) nd = nd + 1
         ! doing same for higher order meshes
         IF (Lm1 .eq. 2) THEN
           nd = nd + Lm1*Nelx + 1
           IF (grid_type .eq. 3) nd = nd - 1
         ENDIF
         IF (Lm1 .eq. 3) THEN
           nd = nd + 2 * (Lm1*Nelx + 1)
           IF (grid_type .eq. 3) nd = nd - 2
         ENDIF
       ENDDO

!       DO i = 1, nc 
!        print *,i,nodeID(1:(Lm1+1)**2,i)
!       ENDDO

       ! Need this for the periodic case
       NelB = nb

       IF (Lnod > 3) RETURN

       ! Paraview viz
       write(8,'(a)') ' '
       write(8,'(a,i8,1x,i8)') 'CELLS ', nc, (Lnod*Lnod+1)*nc
       DO i = 1, nc
         write(8,*) Lnod*Lnod, nodeID(1:Lnod*Lnod,i)-1
       ENDDO
       write(8,'(a)') ' '
       write(8,'(a,i8)') 'CELL_TYPES ',nc
       DO i = 1, nc
         write(8,*) Lnod**3 + 1
       ENDDO
       CLOSE(unit = 8)


END SUBROUTINE SetupMesh


SUBROUTINE SetupICBCandSrc(N, NB, K, L, prob_type)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer N, NB, K, L, prob_type

       integer i, kx, ky, jx, jy
       real*8 x, y, xy, eps, ys, xm, ym

       real*8 theta2,theta1,thetam,thetad,rad2,rad1,radm,radd
       real*8, dimension(Lnod,Lnod) :: xloc, yloc

       eps = 1.d-4

       OPEN(UNIT = 3, file = 'Solution_Nodes.csv', status = 'unknown')

       ! For now, we are assuming a constant value of BC across a single patch at the boundary
       ! That is, the geometry is topologically a square, each side of which has a constant
       ! Neumann or Dirichlet BC (or periodic). The patches are numbered from the bottom and go counter-clockwise
       ! S=1; E=2; N=3; W=4.  The BC switch for these patches is Dirichlet=1, or Neumann=2

       BC_Switch(1:NB) = NeumannBC
       BC_Values(1:K,1:NB) = 0.d0   ! This is legacy and should be removed later

       BC_VelNorm(1:K,1:NB) = 0.d0   ! Wall normal velocity
       BC_VelParl(1:K,1:NB) = 0.d0   ! Wall parallel velocity

       IF (prob_type .eq. 3) THEN

          ! Square Cavity Problem; top velocity is 1

          BC_Switch(1:NB) = DirichletBC
          DO i = 1, NB
            IF (ABS(ycoord(BoundaryNodeID(1,i)) - 1.d0) < 1.d-6 .and. ABS(ycoord(BoundaryNodeID(2,i)) - 1.d0) < 1.d-6) THEN
              BC_VelParl(1:K,i) = 1.d0   ! Wall parallel velocity
            ENDIF
          ENDDO

         Vort0 = 0.d0

       ELSEIF (prob_type .eq. 4) THEN

       ELSE
         print *,'only prob_type = 3 (square cavity) is allowed '
         stop
       ENDIF

       CLOSE(UNIT = 3)

END SUBROUTINE SetupICBCandSrc


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

       integer nt, i, kx,ky,jx,jy
       real*8 x,y
       real*8 dto2, dto3, dto6
       real*8, dimension(Knod, Knod, Nel) :: Vort_tmp, Vort_str, f_of_Vort, psi
       real*8 t1, t2, t3, t4, t5, tlap, tcrl, tdif, tcon

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

       Vort = Vort0             ! leaving Vort0 alone, in case needed elsewhere; generally, we don't need to save in real problems
       Uvel = 0.d0
       Vvel = 0.d0
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
        print *,'timestep ',nt
         Vort_tmp = Vort
         t1 = omp_get_wtime()
         t2 = omp_get_wtime()
print *,'**** DIFFUSION **** '
!!         IF (nt > 1) THEN
           BC_Values = 0.d0
!!           BC_Values = BC_VelParl
           CALL GetLaplacian(HuynhSolver_type, Vort_tmp, psi, A_handle, P_handle, S_handle)                                ! Stage 1
           CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
print *,'Velocity ',minval(sqrt(uvel**2 + vvel**2)),maxval(sqrt(uvel**2 + vvel**2))
!!         ENDIF
!!         BC_Values = BC_VelNorm - VelocJump

!!         CALL GetLaplacian(HuynhSolver_type, Vort0, psi, A_handle, P_handle, S_handle)                                ! Stage 1
!!         CALL GetLaplacGrads(HuynhSolver_type, psi, 0)

print *,'Jump ',minval(VelocJump),maxval(VelocJump)
         VelocJump = VelocJump - BC_VelParl
print *,'Jump ',minval(VelocJump),maxval(VelocJump)
         BC_Values = (VelocJump) * Reyn / dt
!         CALL GetLaplacCurlLo(HuynhSolver_type, psi)
!vort = psi
!CALL dumpResult(nt, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)
!stop
         t3 = omp_get_wtime()
!    print *,'stage 1 ',MaxVal(Sqrt(Uvel**2 + Vvel**2))
         BC_Switch = NeumannBC
         CALL GetDiffusedFlux(HuynhSolver_type, Vort_tmp, f_of_Vort)
         BC_Switch = DirichletBC

   vort = vort + f_of_Vort * dt / Reyn
   vort_tmp = vort
print *,'Post Diff Max Vort ',minval(vort),maxval(vort)

print *,'**** CONVECTION **** '
! CALL dumpResult(nt, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)

if (.true.) then
         t4 = omp_get_wtime()
!         f_of_Vort = f_of_Vort / Reyn
!         IF (prob_type .eq. 2) f_of_Vort = 0.d0
!         BC_Values = BC_VelNorm
           BC_Values = 0.d0
!!           BC_Values = BC_VelParl
           CALL GetLaplacian(HuynhSolver_type, Vort_tmp, psi, A_handle, P_handle, S_handle)                                ! Stage 1
           CALL GetLaplacGrads(HuynhSolver_type, psi, 1)

!!         BC_Values = BC_VelNorm - VelocJump
!!         CALL GetLaplacian(HuynhSolver_type, Vort0, psi, A_handle, P_handle, S_handle)                                ! Stage 1
!!         CALL GetLaplacGrads(HuynhSolver_type, psi, 0)

         BC_Values = BC_VelNorm
         CALL GetConvectedFlux(HuynhSolver_type, Vort_tmp, f_of_Vort)

endif
         t5 = omp_get_wtime()
 tlap = tlap + t2 - t1
 tcrl = tcrl + t3 - t2
 tdif = tdif + t4 - t3
 tcon = tcon + t5 - t4
         IF (tIntegrator_type .eq. 1) THEN

           Vort = Vort + dt*f_of_Vort
print *,'Post Conv Max Vort ',minval(vort),maxval(vort)

         ELSEIF (tIntegrator_type .eq. 2) THEN

           Vort_str = Vort_tmp + dto2*f_of_Vort
         t1 = omp_get_wtime()
           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 2
         t2 = omp_get_wtime()
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
         t3 = omp_get_wtime()
     print *,'stage 2 ',MaxVal(Sqrt(Uvel**2 + Vvel**2))
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

         ELSEIF (tIntegrator_type .eq. 4) THEN

           Vort = Vort + dto6*f_of_Vort
           Vort_str = Vort_tmp + dto2*f_of_Vort
           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 2
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)

           Vort = Vort + dto3*f_of_Vort
           Vort_str = Vort_tmp + dto2*f_of_Vort
           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 3
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)

           Vort = Vort + dto3*f_of_Vort
           Vort_str = Vort_tmp + dt*f_of_Vort
           CALL GetLaplacian(HuynhSolver_type, Vort_str, psi, A_handle, P_handle, S_handle)                              ! Stage 4
           IF (nt > 1) CALL GetLaplacGrads(HuynhSolver_type, psi, 1)
           CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
           CALL GetDiffusedFlux(HuynhSolver_type, Vort_str, f_of_Vort)
           f_of_Vort = f_of_Vort / Reyn
!           IF (prob_type .eq. 2) f_of_Vort = 0.d0
           CALL GetConvectedFlux(HuynhSolver_type, Vort_str, f_of_Vort)

           Vort = Vort + dto6*f_of_Vort

         ENDIF

         IF (MOD(nt,dumpFreq) .eq. 0) CALL dumpResult(nt, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)

       ENDDO
print *,'avg tlap ',0.5*tlap/numStep
print *,'avg tcrl ',0.5*tcrl/numStep
print *,'avg tdif ',0.5*tdif/numStep
print *,'avg tcon ',0.5*tcon/numStep
       IF (MOD(numStep,dumpFreq) .ne. 0) &
           CALL dumpResult(numStep, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)

       OPEN(unit = 8, file = 'Vorticity.vtk', status = 'unknown')
       write(8,'(a)') '# vtk DataFile Version 3.0'
       write(8,'(a)') '2D Unstructured Grid of Quads'
       write(8,'(a)') 'ASCII'
       write(8,'(a)') ' '
!       write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'
       write(8,'(a)') 'DATASET POLYDATA'

       write(8,'(a,i8,a)') 'POINTS ', Knod*Knod*Nel, ' float'
       DO i = 1, Nel
         DO ky = 1, Knod
           DO kx = 1, Knod
             ! interpolate x/y-coord of sol pt at (kx,ky) using nodal coordinates of element i
             x = 0.d0
             y = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                           ! along x-dir                  along y-dir                    x-coord of geom
                 x = x + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * xcoord(nodeID(t2f(jx,jy),i))
                 y = y + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * ycoord(nodeID(t2f(jx,jy),i))
               ENDDO
             ENDDO
             write(8,*) real(x),real(y),' 0.0'
           ENDDO
         ENDDO
       ENDDO

 if(.false.) then
       write(8,'(a)') ' '
       write(8,'(a,i8,1x,i8)') 'CELLS ', Knod*Knod*Nel, 2*Knod*Knod*Nel
       DO i = 1, Knod*Knod*Nel
         write(8,*) "1",i-1
       ENDDO
       write(8,'(a)') ' '
       write(8,'(a,i8)') 'CELL_TYPES ', Knod*Knod*Nel
       DO i = 1, Knod*Knod*Nel
         write(8,*)  "1"
       ENDDO
endif
Vort_tmp = vort
       write(8,'(a)') ' '
       write(8,'(a,i8,a)') 'POINT_DATA ', Knod*Knod*Nel
       write(8,'(a)') 'SCALARS vorticity float 1'
       write(8,'(a)') 'LOOKUP_TABLE default'
       DO i = 1, Nel
         DO ky = 1, Knod
           DO kx = 1, Knod
             write(8,*) real(Vort_tmp(kx,ky,i))
           ENDDO
         ENDDO
       ENDDO
if(.true.) then
       BC_Values = BC_VelParl
       BC_Values = 0.d0
       CALL GetLaplacian(HuynhSolver_type, Vort_tmp, psi, A_handle, P_handle, S_handle)                                ! Stage 1
       CALL GetLaplacGrads(HuynhSolver_type, psi, 1)

!!       BC_Values = BC_VelNorm - VelocJump
!!       CALL GetLaplacian(HuynhSolver_type, Vort0, psi, A_handle, P_handle, S_handle)                                ! Stage 1
!!       CALL GetLaplacGrads(HuynhSolver_type, psi, 0)
endif

       write(8,'(a)') 'VECTORS velocity float'

       DO i = 1, Nel
         DO ky = 1, Knod
           DO kx = 1, Knod
             write(8,*) real(Uvel(kx,ky,i)),real(Vvel(kx,ky,i)),' 0.0'
           ENDDO
         ENDDO
       ENDDO
       CLOSE(unit = 8)

END SUBROUTINE SolveConvectionDiffusion


SUBROUTINE SetLaplacian(HuynhSolver_type, A_handle, P_handle, S_handle)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices
       USE APLLES_Solvers_Module
       USE omp_lib

       USE iso_c_binding

       implicit NONE
       integer HuynhSolver_type

       integer i, j, l, m, p, el, bel, jx, jy, lx, ly, jm1, ij, mj, im, mm1, lm, ml
       integer colc, cole, colw, coln, cols, nnz, nrows, Ksq, K4, ierr
       real*8 gL, tmp, tmpx, tmpy
       real*8 DxDxsi, DyDxsi, DxDeta, DyDeta, FaceJac
       real*8, dimension(Knod) :: Acoef, Bcoef
       real*8, dimension(Knod, Knod) :: SBLB_i0_NGR_j0, SBLB_i0_NGR_j1, SBLB_i1_NGR_j1, SBLB_i1_NGR_j0
       real*8, dimension(Knod, Knod) :: SNGLB_SBLBdotNGR, SNGLB_2xSBLBdotNGR
       real*8, dimension(Knod, Knod) :: SBGLB_SBLB_i0_NGR_j0, SBGLB_SBLB_i0_NGR_j1, SBGLB_SBLB_i1_NGR_j1, SBGLB_SBLB_i1_NGR_j0
       real*8, dimension(Lnod, Lnod) :: xloc, yloc
       real*8, dimension(Knod, Nel) :: FaceAE, FaceBE, FaceAW, FaceBW, FaceAN, FaceBN, FaceAS, FaceBS
       real*8, dimension(Knod, Nel) :: NormAE, NormAW, NormAN, NormAS
       real*8, dimension(Knod, Knod, Nel) :: SolnJac, SolnAxx, SolnAxy, SolnAyy
       real*8, dimension(Knod, Knod) :: NeuMatrix_Orig, NeuMatrix  ! small dense matrix to obtain comx (per element) for a given Neumann BC

       real*8, allocatable, dimension(:,:,:) :: LaplacianCenter, LaplacianEast, LaplacianWest, LaplacianNorth, LaplacianSouth

       integer, dimension(:), allocatable :: rowptr, colidx
       real*8, dimension(:), allocatable :: values

       real*8, dimension(:), allocatable :: b, x

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

       integer, dimension(:), allocatable :: rowptr_filtered, colidx_filtered
       real*8, dimension(:), allocatable :: values_filtered

       real*8 xx,yy,theta2,theta1,thetam,thetad,rad2,rad1,radm,radd

interface
function inv(A) result(Ainv)
  integer dp
  parameter (dp=8)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info
end function inv
end interface

       Ksq = Knod * Knod
       K4 = Ksq * Ksq
       allocate (LaplacianCenter(Ksq,Ksq,Nel), LaplacianEast(Ksq,Ksq,Nel), LaplacianWest(Ksq,Ksq,Nel), &
                                               LaplacianNorth(Ksq,Ksq,Nel), LaplacianSouth(Ksq,Ksq,Nel))

       gL = 0.25d0 * (Knod*Knod)
!       gL = 0.25d0 * (Knod*(Knod+1))

       BndrySrc = 0.d0
       LaplacianCenter = 0.d0
       LaplacianEast = 0.d0
       LaplacianWest = 0.d0
       LaplacianNorth = 0.d0
       LaplacianSouth = 0.d0

       DO i = 1, Knod
         DO j = 1, Knod
           ! Note: SBLB_i0_NGR_j0(i , j) + SBLB_i1_NGR_j1(Knod+1-i , Knod+1-j) = 0
           !       SBLB_i0_NGR_j1(i , j) + SBLB_i1_NGR_j0(Knod+1-i , Knod+1-j) = 0
           !       But we don't take advantage of this, if there is any

           ! E00(i,j), E01(i,j), E11(i,j), E10(i,j) in notes
           SBLB_i0_NGR_j0(i,j) = 0.5d0 * SolnBndryLgrangeBasis(i,0) * NodesGradRadau(j,0)
           SBLB_i0_NGR_j1(i,j) = 0.5d0 * SolnBndryLgrangeBasis(i,0) * NodesGradRadau(j,1)
           SBLB_i1_NGR_j1(i,j) = 0.5d0 * SolnBndryLgrangeBasis(i,1) * NodesGradRadau(j,1)
           SBLB_i1_NGR_j0(i,j) = 0.5d0 * SolnBndryLgrangeBasis(i,1) * NodesGradRadau(j,0)
           ! C(i,j) in notes
           SNGLB_SBLBdotNGR(i,j) = SolnNodesGradLgrangeBasis(i,j) - (SBLB_i0_NGR_j0(i,j) + SBLB_i1_NGR_j1(i,j))
           ! D(i,j) in notes
           SNGLB_2xSBLBdotNGR(i,j) = SolnNodesGradLgrangeBasis(i,j) - 2.d0*(SBLB_i0_NGR_j0(i,j) + SBLB_i1_NGR_j1(i,j))
           ! F00(i,j), F01(i,j), F11(i,j), F10(i,j) in notes
           SBGLB_SBLB_i0_NGR_j0(i,j) = 0.5d0 * (SolnBndryGradLgrangeBasis(i,0) + gL*SolnBndryLgrangeBasis(i,0)) * &
                                                NodesGradRadau(j,0)
           SBGLB_SBLB_i0_NGR_j1(i,j) = 0.5d0 * (SolnBndryGradLgrangeBasis(i,0) + gL*SolnBndryLgrangeBasis(i,0)) * &
                                                NodesGradRadau(j,1)
           SBGLB_SBLB_i1_NGR_j1(i,j) = 0.5d0 * (SolnBndryGradLgrangeBasis(i,1) - gL*SolnBndryLgrangeBasis(i,1)) * &
                                                NodesGradRadau(j,1)
           SBGLB_SBLB_i1_NGR_j0(i,j) = 0.5d0 * (SolnBndryGradLgrangeBasis(i,1) - gL*SolnBndryLgrangeBasis(i,1)) * &
                                                NodesGradRadau(j,0)
         ENDDO
       ENDDO

       DO el = 1, Nel

         DO ly = 1, Lnod
           DO lx = 1, Lnod
             xloc(lx,ly) = xcoord(nodeID(t2f(lx,ly),el))
             yloc(lx,ly) = ycoord(nodeID(t2f(lx,ly),el))
           ENDDO
         ENDDO

         DO jy = 1, Knod
           ! Geometric metric stuff at all solution nodes, Soln (per element)
           DO jx = 1, Knod
             DxDxsi = 0.d0
             DyDxsi = 0.d0
             DxDeta = 0.d0
             DyDeta = 0.d0
!xx = 0.d0
!yy = 0.d0
             DO ly = 1, Lnod
               DO lx = 1, Lnod
                 DxDxsi = DxDxsi + &
                              ! grad at x-dir                    no grad at y-dir               x/y-coord of geom
                              GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
                 DyDxsi = DyDxsi + &
                              GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
                 DxDeta = DxDeta + &
                              ! no grad at x-dir             grad at y-dir                      x/y-coord of geom
                              GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
                 DyDeta = DyDeta + &
                              GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
!  xx = xx + GeomNodesLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
!  yy = yy + GeomNodesLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
               ENDDO
             ENDDO

 if (exact) then
! print *,el,jx,jy,xx,yy,sqrt(xx**2+yy**2),atan2(yy,xx),DxDxsi,DyDxsi,DxDeta,DyDeta, &
! print *,el,jx,jy,DxDxsi,DyDxsi,DxDeta,DyDeta,(DxDxsi * DyDeta - DxDeta * DyDxsi)
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             DxDxsi = -thetad * (radm + sps(jy)*radd) * sin(thetam + sps(jx)*thetad)
             DyDxsi =  thetad * (radm + sps(jy)*radd) * cos(thetam + sps(jx)*thetad)
             DxDeta =  radd * cos(thetam + sps(jx)*thetad)
             DyDeta =  radd * sin(thetam + sps(jx)*thetad)

! print *,el,jx,jy,-thetad*(radm+sps(jy)*radd)*sin(thetam+sps(jx)*thetad), &
!                   thetad*(radm+sps(jy)*radd)*cos(thetam+sps(jx)*thetad), &
!                   radd*cos(thetam+sps(jx)*thetad), &
!                   radd*sin(thetam+sps(jx)*thetad), &
!                  -radd*thetad*(radm+sps(jy)*radd)
! print *,el,jx,jy,100*(1.d0-DxDxsi/(-thetad*(radm+sps(jy)*radd)*sin(thetam+sps(jx)*thetad))), &
!                  100*(1.d0-DyDxsi/(thetad*(radm+sps(jy)*radd)*cos(thetam+sps(jx)*thetad))), &
!                  100*(1.d0-DxDeta/(radd*cos(thetam+sps(jx)*thetad))), &
!                  100*(1.d0-DyDeta/(radd*sin(thetam+sps(jx)*thetad))), &
!                  100*(1.d0-(DxDxsi * DyDeta - DxDeta * DyDxsi)/(-radd*thetad*(radm+sps(jy)*radd)))
! print *,el,jx,jy,xloc(1,1),yloc(1,1),xloc(Lnod,1),yloc(Lnod,1),theta2,theta1,rad2,rad1
! print *,el,jx,jy,xloc(1,1),yloc(1,1),xloc(Lnod,1),yloc(Lnod,1),theta2,theta1,thetad   
! print *,el,jx,jy,xx,yy,sqrt(xx**2+yy**2),atan2(yy,xx),radm+sps(jy)*radd,thetam+sps(jx)*thetad
! print *,' '
 endif

             SolnJac(jx,jy,el) = 1.d0 / (DxDxsi * DyDeta - DxDeta * DyDxsi)
             ! Ax in notes
             SolnAxx(jx,jy,el) = (DxDeta * DxDeta + DyDeta * DyDeta) * SolnJac(jx,jy,el)
             ! B in notes   --  should try Axy and Ayx in (jx,jy) and (jy,jx) storage form
             SolnAxy(jx,jy,el) = (DxDxsi * DxDeta + DyDxsi * DyDeta) * SolnJac(jx,jy,el)
             ! Ay in notes  --  should try saving in (jy,jx) form to speed up a bit
             SolnAyy(jx,jy,el) = (DxDxsi * DxDxsi + DyDxsi * DyDxsi) * SolnJac(jx,jy,el)
           ENDDO

           ! Right Face: Geometric metric stuff (per element)
           DxDxsi = 0.d0
           DyDxsi = 0.d0
           DxDeta = 0.d0
           DyDeta = 0.d0
           DO ly = 1, Lnod
             DO lx = 1, Lnod
               DxDxsi = DxDxsi + &
                                 ! grad at x-dir                   no grad along y-dir            x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(lx,1) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
               DyDxsi = DyDxsi + &
                                 GeomBndryGradLgrangeBasis(lx,1) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
               DxDeta = DxDeta + &
                                 ! no grad at x-dir            grad along y-dir                   x/y-coord of geom
                                 GeomBndryLgrangeBasis(lx,1) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
               DyDeta = DyDeta + &
                                 GeomBndryLgrangeBasis(lx,1) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
             ENDDO
           ENDDO

!print *,'E el,jx,jy, dx,dy ',el,jx,jy,dxdxsi,dydxsi,dxdeta,dydeta

 if (exact) then
             DxDxsi = -thetad * (radm + sps(jy)*radd) * sin(theta2)
             DyDxsi =  thetad * (radm + sps(jy)*radd) * cos(theta2)
             DxDeta =  radd * cos(theta2)
             DyDeta =  radd * sin(theta2)
 endif

           FaceJac = 1.d0 / (DxDxsi * DyDeta - DxDeta * DyDxsi)
           FaceAE(jy,el) = (DxDeta * DxDeta + DyDeta * DyDeta) * FaceJac
           FaceBE(jy,el) = (DxDxsi * DxDeta + DyDxsi * DyDeta) * FaceJac
           NormAE(jy,el) = Sqrt(DxDeta * DxDeta + DyDeta * DyDeta)

           ! Left Face: Geometric metric stuff (per element)
           DxDxsi = 0.d0
           DyDxsi = 0.d0
           DxDeta = 0.d0
           DyDeta = 0.d0
           DO ly = 1, Lnod
             DO lx = 1, Lnod
               DxDxsi = DxDxsi + &
                                 ! grad at x-dir                   no grad along y-dir            x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(lx,0) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
               DyDxsi = DyDxsi + &
                                 GeomBndryGradLgrangeBasis(lx,0) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
               DxDeta = DxDeta + &
                                 ! no grad at x-dir            grad along y-dir                   x/y-coord of geom
                                 GeomBndryLgrangeBasis(lx,0) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
               DyDeta = DyDeta + &
                                 GeomBndryLgrangeBasis(lx,0) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
             ENDDO
           ENDDO
!print *,'W el,jx,jy, dx,dy ',el,jx,jy,dxdxsi,dydxsi,dxdeta,dydeta

 if (exact) then
             DxDxsi = -thetad * (radm + sps(jy)*radd) * sin(theta1)
             DyDxsi =  thetad * (radm + sps(jy)*radd) * cos(theta1)
             DxDeta =  radd * cos(theta1)
             DyDeta =  radd * sin(theta1)
 endif

           FaceJac = 1.d0 / (DxDxsi * DyDeta - DxDeta * DyDxsi)
           FaceAW(jy,el) = (DxDeta * DxDeta + DyDeta * DyDeta) * FaceJac
           FaceBW(jy,el) = (DxDxsi * DxDeta + DyDxsi * DyDeta) * FaceJac
           NormAW(jy,el) = Sqrt(DxDeta * DxDeta + DyDeta * DyDeta)

           ! North Face: Geometric metric stuff (per element)
           DxDxsi = 0.d0
           DyDxsi = 0.d0
           DxDeta = 0.d0
           DyDeta = 0.d0
           DO ly = 1, Lnod
             DO lx = 1, Lnod
               DxDeta = DxDeta + &
                                 ! grad at x-dir                   no grad along y-dir            x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(ly,1) * GeomNodesLgrangeBasis(lx,jy) * xloc(lx,ly)
               DyDeta = DyDeta + &
                                 GeomBndryGradLgrangeBasis(ly,1) * GeomNodesLgrangeBasis(lx,jy) * yloc(lx,ly)
               DxDxsi = DxDxsi + &
                                 ! no grad at x-dir            grad along y-dir                   x/y-coord of geom
                                 GeomBndryLgrangeBasis(ly,1) * GeomNodesGradLgrangeBasis(lx,jy) * xloc(lx,ly)
               DyDxsi = DyDxsi + &
                                 GeomBndryLgrangeBasis(ly,1) * GeomNodesGradLgrangeBasis(lx,jy) * yloc(lx,ly)
             ENDDO
           ENDDO

 if (exact) then
             DxDxsi = -thetad * (rad2) * sin(thetam + sps(jy)*thetad)
             DyDxsi =  thetad * (rad2) * cos(thetam + sps(jy)*thetad)
             DxDeta =  radd * cos(thetam + sps(jy)*thetad)
             DyDeta =  radd * sin(thetam + sps(jy)*thetad)
 endif

           FaceJac = 1.d0 / (DxDxsi * DyDeta - DxDeta * DyDxsi)
           FaceAN(jy,el) = (DxDxsi * DxDxsi + DyDxsi * DyDxsi) * FaceJac
           FaceBN(jy,el) = (DxDxsi * DxDeta + DyDxsi * DyDeta) * FaceJac
           NormAN(jy,el) = Sqrt(DxDxsi * DxDxsi + DyDxsi * DyDxsi)

           ! South Face: Geometric metric stuff (per element)
           DxDxsi = 0.d0
           DyDxsi = 0.d0
           DxDeta = 0.d0
           DyDeta = 0.d0
           DO ly = 1, Lnod
             DO lx = 1, Lnod
               DxDeta = DxDeta + &
                                 ! grad at x-dir                   no grad along y-dir            x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(ly,0) * GeomNodesLgrangeBasis(lx,jy) * xloc(lx,ly)
               DyDeta = DyDeta + &
                                 GeomBndryGradLgrangeBasis(ly,0) * GeomNodesLgrangeBasis(lx,jy) * yloc(lx,ly)
               DxDxsi = DxDxsi + &
                                 ! no grad at x-dir            grad along y-dir                   x/y-coord of geom
                                 GeomBndryLgrangeBasis(ly,0) * GeomNodesGradLgrangeBasis(lx,jy) * xloc(lx,ly)
               DyDxsi = DyDxsi + &
                                 GeomBndryLgrangeBasis(ly,0) * GeomNodesGradLgrangeBasis(lx,jy) * yloc(lx,ly)
             ENDDO
           ENDDO
!print *,'S el,jx,jy, dx,dy ',el,jx,jy,dxdxsi,dydxsi,dxdeta,dydeta

 if (exact) then
             DxDxsi = -thetad * (rad1) * sin(thetam + sps(jy)*thetad)
             DyDxsi =  thetad * (rad1) * cos(thetam + sps(jy)*thetad)
             DxDeta =  radd * cos(thetam + sps(jy)*thetad)
             DyDeta =  radd * sin(thetam + sps(jy)*thetad)
 endif

           FaceJac = 1.d0 / (DxDxsi * DyDeta - DxDeta * DyDxsi)
           FaceAS(jy,el) = (DxDxsi * DxDxsi + DyDxsi * DyDxsi) * FaceJac
           FaceBS(jy,el) = (DxDxsi * DxDeta + DyDxsi * DyDeta) * FaceJac
           NormAS(jy,el) = Sqrt(DxDxsi * DxDxsi + DyDxsi * DyDxsi)

         ENDDO
!print *,' '
       ENDDO


       DO el = 1, Nel
         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
             DO m = 1, Knod
               mm1 = (m-1) * Knod
               mj = jm1 + m
               im = mm1 + i
               tmpx = 0.d0
               tmpy = 0.d0
               DO l = 1, Knod
                 tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SNGLB_SBLBdotNGR(m,l)
                 tmpy = tmpy + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SNGLB_SBLBdotNGR(m,l)
               ENDDO
               LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) + tmpx
               LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) + tmpy
               DO l = 1, Knod
                 lm = mm1 + l
                 LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) - ( &
                                               SolnAxy(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SNGLB_SBLBdotNGR(m,j) + &
                                               SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * SNGLB_SBLBdotNGR(l,i) )
               ENDDO
             ENDDO
           ENDDO
         ENDDO

         ! Right Face
         IF (elemID(2,el) .lt. 0) THEN

           IF (BC_Switch(-elemID(2,el)) .eq. DirichletBC) THEN

             bel = -elemID(2,el)
             DO j = 1, Knod
               jm1 = (j-1) * Knod
!               tmp = 2.d0 * gL * FaceAE(j,el) * BC_Values(j,-elemID(2,el))
               tmp = 2.d0 * gL * FaceAE(j,el)
               DO i = 1, Knod
                 ij = jm1 + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + tmp * NodesGradRadau(i,1)
                 BndrySrc(j,ij,bel) = BndrySrc(j,ij,bel) + tmp * NodesGradRadau(i,1)
                 DO m = 1, Knod
                   mm1 = (m-1) * Knod
                   mj = jm1 + m
!                   BndrySrc(ij,el) = BndrySrc(ij,el) + &
!                                       SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,1) &
!                                                                                 * BC_Values(j,-elemID(2,el)) - &
!                                       NodesGradRadau(i,1) * BC_Values(m,-elemID(2,el)) &
!                                                           * ( SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) + &
!                                                               FaceBE(j,el) * SolnNodesGradLgrangeBasis(m,j) )
                   BndrySrc(j,ij,bel) = BndrySrc(j,ij,bel) + &
                                        SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,1)
                   BndrySrc(m,ij,bel) = BndrySrc(m,ij,bel) - NodesGradRadau(i,1) &
                                                        * ( SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) + &
                                                            FaceBE(j,el) * SolnNodesGradLgrangeBasis(m,j) )
                   tmpx = 0.d0
                   DO l = 1, Knod
                     tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i1_NGR_j1(m,l)
                   ENDDO
                   LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) - tmpx + 2.d0 * FaceAE(j,el) * &
                                                 ( SBGLB_SBLB_i1_NGR_j1(m,i) - gL * SBLB_i1_NGR_j1(m,i) )
                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j)
                   DO l = 1, Knod
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) + tmpy * SBLB_i1_NGR_j1(l,i)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ELSEIF (BC_Switch(-elemID(2,el)) .eq. NeumannBC) THEN

             DO j = 1, Knod
               DO i = 1, Knod
                 NeuMatrix_Orig(i,j) = FaceBE(i,el) * SolnNodesGradLgrangeBasis(j,i)
               ENDDO
               NeuMatrix_Orig(j,j) = NeuMatrix_Orig(j,j) - 2.d0 * gL * FaceAE(j,el)
             ENDDO

             NeuMatrix = inv(NeuMatrix_Orig)

             bel = -elemID(2,el)
             DO j = 1, Knod
               jm1 = (j-1) * Knod
               DO i = 1, Knod
                 ij = jm1 + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + NodesGradRadau(i,1) * BC_Values(j,-elemID(2,el)) * NormAE(j,el)
                 BndrySrc(j,ij,bel) = BndrySrc(j,ij,bel) + NodesGradRadau(i,1) * NormAE(j,el)
                 DO m = 1, Knod
                   mm1 = (m-1) * Knod
                   mj = jm1 + m

!                   tmpx = 0.d0
!                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(i,1)
!                   DO l = 1, Knod
!                     tmpx = tmpx + NeuMatrix(j,l) * BC_Values(l,-elemID(2,el)) * NormAE(l,el)
!                     BndrySrc(ij,el) = BndrySrc(ij,el) + tmpy * NeuMatrix(m,l) * BC_Values(l,-elemID(2,el)) * NormAE(l,el)
!                   ENDDO
!                   BndrySrc(ij,el) = BndrySrc(ij,el) - tmpx * SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,1)
                   tmpx = SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,1)
                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(i,1)
                   DO l = 1, Knod
                     BndrySrc(l,ij,bel) = BndrySrc(l,ij,bel) + (tmpy * NeuMatrix(m,l) - tmpx * NeuMatrix(j,l)) * NormAE(l,el)
                   ENDDO

                   tmpx = 0.d0
                   DO l = 1, Knod
                     tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i1_NGR_j1(m,l)
                   ENDDO
                   LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) - tmpx

                   tmpy =0.d0
                   DO l = 1, Knod
                     tmpy = tmpy + SolnAxy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * NeuMatrix(l,m)
                   ENDDO
                   DO l = 1, Knod
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) - 2.d0 * tmpy * FaceAE(m,el) * &
                                                                            (SBGLB_SBLB_i1_NGR_j1(l,i) - gL * SBLB_i1_NGR_j1(l,i))
                   ENDDO

                   DO l = 1, Knod
                     tmpx = 0.d0
                     DO p = 1, Knod
                       tmpx = tmpx + SolnAxx(p,j,el) * SNGLB_2xSBLBdotNGR(p,i) * &
                                                                            (SBGLB_SBLB_i1_NGR_j1(l,p) - gL * SBLB_i1_NGR_j1(l,p))
                     ENDDO
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) + 2.d0 * tmpx * FaceAE(m,el) * NeuMatrix(j,m)
                   ENDDO

                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j)
                   DO l = 1, Knod
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) + tmpy * SBLB_i1_NGR_j1(l,i)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ENDIF

         ELSE

           DO j = 1, Knod
             Acoef(j) = FaceAW(j,elemID(2,el))
             Bcoef(j) = 0.5d0 * (FaceBE(j,el) + FaceBW(j,elemID(2,el)))
           ENDDO

           DO j = 1, Knod
             jm1 = (j-1) * Knod
             DO i = 1, Knod
               ij = jm1 + i
               DO m = 1, Knod
                 mm1 = (m-1) * Knod
                 mj = jm1 + m
                 tmpx = 0.d0
                 DO l = 1, Knod
                   tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i0_NGR_j1(m,l)
                 ENDDO
                 LaplacianEast(mj,ij,el) = LaplacianEast(mj,ij,el) + tmpx + &
                                               Acoef(j) * SBGLB_SBLB_i0_NGR_j1(m,i) + &
                                               gL * FaceAE(j,el) * SBLB_i0_NGR_j1(m,i)
                 LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) + ( &
                                               FaceAE(j,el) * SBGLB_SBLB_i1_NGR_j1(m,i) - &
                                               gL * Acoef(j) * SBLB_i1_NGR_j1(m,i) )
                 tmp = Bcoef(j) * SolnNodesGradLgrangeBasis(m,j)
                 tmpy = tmp + SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j)
                 DO l = 1, Knod
                   lm = mm1 + l
                   LaplacianEast(lm,ij,el) = LaplacianEast(lm,ij,el) - tmpy * SBLB_i0_NGR_j1(l,i) 
                   LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) - tmp * SBLB_i1_NGR_j1(l,i)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

         ENDIF

         ! Left Face
         IF (elemID(4,el) .lt. 0) THEN

           IF (BC_Switch(-elemID(4,el)) .eq. DirichletBC) THEN

             bel = -elemID(4,el)
             DO j = 1, Knod
               jm1 = (j-1) * Knod
!               tmp = - 2.d0 * gL * FaceAW(j,el) * BC_Values(j,-elemID(4,el))
               tmp = - 2.d0 * gL * FaceAW(j,el)
               DO i = 1, Knod
                 ij = jm1 + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + tmp * NodesGradRadau(i,0)
                 BndrySrc(j,ij,bel) = BndrySrc(j,ij,bel) + tmp * NodesGradRadau(i,0)
                 DO m = 1, Knod
                   mm1 = (m-1) * Knod
                   mj = jm1 + m
!                   BndrySrc(ij,el) = BndrySrc(ij,el) + &
!                                       SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,0) &
!                                                                                 * BC_Values(j,-elemID(4,el)) - &
!                                       NodesGradRadau(i,0) * BC_Values(m,-elemID(4,el)) &
!                                                           * ( SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) + &
!                                                               FaceBW(j,el) * SolnNodesGradLgrangeBasis(m,j) )
                   BndrySrc(j,ij,bel) = BndrySrc(j,ij,bel) + &
                                        SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,0)
                   BndrySrc(m,ij,bel) = BndrySrc(m,ij,bel) - NodesGradRadau(i,0) &
                                                        * ( SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) + &
                                                            FaceBW(j,el) * SolnNodesGradLgrangeBasis(m,j) )
                   tmpx = 0.d0
                   DO l = 1, Knod
                     tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i0_NGR_j0(m,l)
                   ENDDO
                   LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) - tmpx + 2.d0 * FaceAW(j,el) * &
                                                 ( SBGLB_SBLB_i0_NGR_j0(m,i) + gL * SBLB_i0_NGR_j0(m,i) )
                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j)
                   DO l = 1, Knod
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) + tmpy * SBLB_i0_NGR_j0(l,i)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ELSEIF (BC_Switch(-elemID(4,el)) .eq. NeumannBC) THEN

             DO j = 1, Knod
               DO i = 1, Knod
                 NeuMatrix_Orig(i,j) = FaceBW(i,el) * SolnNodesGradLgrangeBasis(j,i)
               ENDDO
               NeuMatrix_Orig(j,j) = NeuMatrix_Orig(j,j) + 2.d0 * gL * FaceAW(j,el)
             ENDDO

             NeuMatrix = inv(NeuMatrix_Orig)

             bel = -elemID(4,el)
             DO j = 1, Knod
               jm1 = (j-1) * Knod
               DO i = 1, Knod
                 ij = jm1 + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + NodesGradRadau(i,0) * BC_Values(j,-elemID(4,el)) * NormAW(j,el)
                 BndrySrc(j,ij,bel) = BndrySrc(j,ij,bel) + NodesGradRadau(i,0) * NormAW(j,el)
                 DO m = 1, Knod
                   mm1 = (m-1) * Knod
                   mj = jm1 + m

!                   tmpx = 0.d0
!                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(i,0)
!                   DO l = 1, Knod
!                     tmpx = tmpx + NeuMatrix(j,l) * BC_Values(l,-elemID(4,el)) * NormAW(l,el)
!                     BndrySrc(ij,el) = BndrySrc(ij,el) + tmpy * NeuMatrix(m,l) * BC_Values(l,-elemID(4,el)) * NormAW(l,el)
!                   ENDDO
!                   BndrySrc(ij,el) = BndrySrc(ij,el) - tmpx * SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,0)
                   tmpx = SolnAxx(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(m,0)
                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(i,0)
                   DO l = 1, Knod
                     BndrySrc(l,ij,bel) = BndrySrc(l,ij,bel) + (tmpy * NeuMatrix(m,l) - tmpx * NeuMatrix(j,l)) * NormAW(l,el)
                   ENDDO

                   tmpx = 0.d0
                   DO l = 1, Knod
                     tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i0_NGR_j0(m,l)
                   ENDDO
                   LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) - tmpx

                   tmpy =0.d0
                   DO l = 1, Knod
                     tmpy = tmpy + SolnAxy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * NeuMatrix(l,m)
                   ENDDO
                   DO l = 1, Knod
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) - 2.d0 * tmpy * FaceAW(m,el) * &
                                                                            (SBGLB_SBLB_i0_NGR_j0(l,i) + gL * SBLB_i0_NGR_j0(l,i))
                   ENDDO

                   DO l = 1, Knod
                     tmpx = 0.d0
                     DO p = 1, Knod
                       tmpx = tmpx + SolnAxx(p,j,el) * SNGLB_2xSBLBdotNGR(p,i) * &
                                                                            (SBGLB_SBLB_i0_NGR_j0(l,p) + gL * SBLB_i0_NGR_j0(l,p))
                     ENDDO
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) + 2.d0 * tmpx * FaceAW(m,el) * NeuMatrix(j,m)
                   ENDDO

                   tmpy = SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j)
                   DO l = 1, Knod
                     lm = mm1 + l
                     LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) + tmpy * SBLB_i0_NGR_j0(l,i)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ENDIF

         ELSE

           DO j = 1, Knod
             Acoef(j) = FaceAE(j,elemID(4,el))
             Bcoef(j) = 0.5d0 * (FaceBW(j,el) + FaceBE(j,elemID(4,el)))
           ENDDO

           DO j = 1, Knod
             jm1 = (j-1) * Knod
             DO i = 1, Knod
               ij = jm1 + i
               DO m = 1, Knod
                 mm1 = (m-1) * Knod
                 mj = jm1 + m
                 tmpx = 0.d0
                 DO l = 1, Knod
                   tmpx = tmpx + SolnAxx(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * SBLB_i1_NGR_j0(m,l)
                 ENDDO
                 LaplacianWest(mj,ij,el) = LaplacianWest(mj,ij,el) + tmpx + &
                                               Acoef(j) * SBGLB_SBLB_i1_NGR_j0(m,i) - &
                                               gL * FaceAW(j,el) * SBLB_i1_NGR_j0(m,i)
                 LaplacianCenter(mj,ij,el) = LaplacianCenter(mj,ij,el) + ( &
                                               FaceAW(j,el) * SBGLB_SBLB_i0_NGR_j0(m,i) + &
                                               gL * Acoef(j) * SBLB_i0_NGR_j0(m,i) )
                 tmp = Bcoef(j) * SolnNodesGradLgrangeBasis(m,j)
                 tmpy = tmp + SolnAxy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j)
                 DO l = 1, Knod
                   lm = mm1 + l
                   LaplacianWest(lm,ij,el) = LaplacianWest(lm,ij,el) - tmpy * SBLB_i1_NGR_j0(l,i) 
                   LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) - tmp * SBLB_i0_NGR_j0(l,i)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

         ENDIF

         ! North Face
         IF (elemID(3,el) .lt. 0) THEN

           IF (BC_Switch(-elemID(3,el)) .eq. DirichletBC) THEN

             bel = -elemID(3,el)
             DO i = 1, Knod
!               tmp = 2.d0 * gL * FaceAN(i,el) * BC_Values(i,-elemID(3,el))
               tmp = 2.d0 * gL * FaceAN(i,el)
               DO j = 1, Knod
                 ij = (j-1) * Knod + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + tmp * NodesGradRadau(j,1)
                 BndrySrc(i,ij,bel) = BndrySrc(i,ij,bel) + tmp * NodesGradRadau(j,1)
                 DO m = 1, Knod
                   im = (m-1) * Knod + i
!                   BndrySrc(ij,el) = BndrySrc(ij,el) + &
!                                       SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,1) &
!                                                                                 * BC_Values(i,-elemID(3,el)) - &
!                                       NodesGradRadau(j,1) * BC_Values(m,-elemID(3,el)) &
!                                                           * ( SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) + &
!                                                               FaceBN(i,el) * SolnNodesGradLgrangeBasis(m,i) )
                   BndrySrc(i,ij,bel) = BndrySrc(i,ij,bel) + &
                                       SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,1)
                   BndrySrc(m,ij,bel) = BndrySrc(m,ij,bel) - NodesGradRadau(j,1) &
                                                       * ( SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) + &
                                                           FaceBN(i,el) * SolnNodesGradLgrangeBasis(m,i) )
                   tmpy = 0.d0
                   DO l = 1, Knod
                     tmpy = tmpy + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SBLB_i1_NGR_j1(m,l)
                   ENDDO
                   LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) - tmpy + 2.d0 * FaceAN(i,el) * &
                                                 ( SBGLB_SBLB_i1_NGR_j1(m,j) - gL * SBLB_i1_NGR_j1(m,j) )
                   tmpx = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i)
                   DO l = 1, Knod
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) + tmpx * SBLB_i1_NGR_j1(l,j)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ELSEIF (BC_Switch(-elemID(3,el)) .eq. NeumannBC) THEN

             DO j = 1, Knod
               DO i = 1, Knod
                 NeuMatrix_Orig(i,j) = FaceBN(i,el) * SolnNodesGradLgrangeBasis(j,i)
               ENDDO
               NeuMatrix_Orig(j,j) = NeuMatrix_Orig(j,j) - 2.d0 * gL * FaceAN(j,el)
             ENDDO

             NeuMatrix = inv(NeuMatrix_Orig)

             bel = -elemID(3,el)
             DO i = 1, Knod
               DO j = 1, Knod
                 ij = (j-1) * Knod + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + NodesGradRadau(j,1) * BC_Values(i,-elemID(3,el)) * NormAN(i,el)
                 BndrySrc(i,ij,bel) = BndrySrc(i,ij,bel) + NodesGradRadau(j,1) * NormAN(i,el)
                 DO m = 1, Knod
                   im = (m-1) * Knod + i

!                   tmpx = 0.d0
!                   tmpy = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(j,1)
!                   DO l = 1, Knod
!                     tmpx = tmpx + NeuMatrix(i,l) * BC_Values(l,-elemID(3,el)) * NormAN(l,el)
!                     BndrySrc(ij,el) = BndrySrc(ij,el) + tmpy * NeuMatrix(m,l) * BC_Values(l,-elemID(3,el)) * NormAN(l,el)
!                   ENDDO
!                   BndrySrc(ij,el) = BndrySrc(ij,el) - tmpx * SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,1)
                   tmpx = SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,1)
                   tmpy = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(j,1)
                   DO l = 1, Knod
                     BndrySrc(l,ij,bel) = BndrySrc(l,ij,bel) + (tmpy * NeuMatrix(m,l) - tmpx * NeuMatrix(i,l)) * NormAN(l,el)
                   ENDDO

                   tmpx = 0.d0
                   DO l = 1, Knod
                     tmpx = tmpx + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SBLB_i1_NGR_j1(m,l)
                   ENDDO
                   LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) - tmpx

                   tmpy =0.d0
                   DO l = 1, Knod
                     tmpy = tmpy + SolnAxy(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * NeuMatrix(l,m)
                   ENDDO
                   DO l = 1, Knod
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) - 2.d0 * tmpy * FaceAN(m,el) * &
                                                                            (SBGLB_SBLB_i1_NGR_j1(l,j) - gL * SBLB_i1_NGR_j1(l,j))
                   ENDDO

                   DO l = 1, Knod
                     tmpx = 0.d0
                     DO p = 1, Knod
                       tmpx = tmpx + SolnAyy(i,p,el) * SNGLB_2xSBLBdotNGR(p,j) * &
                                                                            (SBGLB_SBLB_i1_NGR_j1(l,p) - gL * SBLB_i1_NGR_j1(l,p))
                     ENDDO
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) + 2.d0 * tmpx * FaceAN(m,el) * NeuMatrix(i,m)
                   ENDDO

                   tmpy = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i)
                   DO l = 1, Knod
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) + tmpy * SBLB_i1_NGR_j1(l,j)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ENDIF

         ELSE

           DO j = 1, Knod
             Acoef(j) = FaceAS(j,elemID(3,el))
             Bcoef(j) = 0.5d0 * (FaceBN(j,el) + FaceBS(j,elemID(3,el)))
           ENDDO

           DO i = 1, Knod
             DO j = 1, Knod
               ij = (j-1) * Knod + i
               DO m = 1, Knod
                 im = (m-1) * Knod + i
                 tmpy = 0.d0
                 DO l = 1, Knod
                   tmpy = tmpy + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SBLB_i0_NGR_j1(m,l)
                 ENDDO
                 LaplacianNorth(im,ij,el) = LaplacianNorth(im,ij,el) + tmpy + &
                                               Acoef(i) * SBGLB_SBLB_i0_NGR_j1(m,j) + &
                                               gL * FaceAN(i,el) * SBLB_i0_NGR_j1(m,j)
                 LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) + ( &
                                               FaceAN(i,el) * SBGLB_SBLB_i1_NGR_j1(m,j) - &
                                               gL * Acoef(i) * SBLB_i1_NGR_j1(m,j) )
                 tmp = Bcoef(i) * SolnNodesGradLgrangeBasis(m,i)
                 tmpx = tmp + SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i)
                 DO l = 1, Knod
                   ml = (l-1) * Knod + m
                   LaplacianNorth(ml,ij,el) = LaplacianNorth(ml,ij,el) - tmpx * SBLB_i0_NGR_j1(l,j) 
                   LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) - tmp * SBLB_i1_NGR_j1(l,j)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

         ENDIF

         ! South Face
         IF (elemID(1,el) .lt. 0) THEN

           IF (BC_Switch(-elemID(1,el)) .eq. DirichletBC) THEN

             bel = -elemID(1,el)
             DO i = 1, Knod
!               tmp = - 2.d0 * gL * FaceAS(i,el) * BC_Values(i,-elemID(1,el))
               tmp = - 2.d0 * gL * FaceAS(i,el)
               DO j = 1, Knod
                 ij = (j-1) * Knod + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + tmp * NodesGradRadau(j,0)
                 BndrySrc(i,ij,bel) = BndrySrc(i,ij,bel) + tmp * NodesGradRadau(j,0)
                 DO m = 1, Knod
                   im = (m-1) * Knod + i
!                   BndrySrc(ij,el) = BndrySrc(ij,el) + &
!                                       SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,0) &
!                                                                                 * BC_Values(i,-elemID(1,el)) - &
!                                       NodesGradRadau(j,0) * BC_Values(m,-elemID(1,el)) &
!                                                           * ( SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) + &
!                                                               FaceBS(i,el) * SolnNodesGradLgrangeBasis(m,i) )
                   BndrySrc(i,ij,bel) = BndrySrc(i,ij,bel) + &
                                        SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,0)
                   BndrySrc(m,ij,bel) = BndrySrc(m,ij,bel) - NodesGradRadau(j,0) &
                                                        * ( SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) + &
                                                            FaceBS(i,el) * SolnNodesGradLgrangeBasis(m,i) )
                   tmpy = 0.d0
                   DO l = 1, Knod
                     tmpy = tmpy + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SBLB_i0_NGR_j0(m,l)
                   ENDDO
                   LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) - tmpy + 2.d0 * FaceAS(i,el) * &
                                                 ( SBGLB_SBLB_i0_NGR_j0(m,j) + gL * SBLB_i0_NGR_j0(m,j) )
                   tmpx = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i)
                   DO l = 1, Knod
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) + tmpx * SBLB_i0_NGR_j0(l,j)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ELSEIF (BC_Switch(-elemID(1,el)) .eq. NeumannBC) THEN

             DO j = 1, Knod
               DO i = 1, Knod
                 NeuMatrix_Orig(i,j) = FaceBS(i,el) * SolnNodesGradLgrangeBasis(j,i)
               ENDDO
               NeuMatrix_Orig(j,j) = NeuMatrix_Orig(j,j) + 2.d0 * gL * FaceAS(j,el)
             ENDDO

             NeuMatrix = inv(NeuMatrix_Orig)

             bel = -elemID(1,el)
             DO i = 1, Knod
               DO j = 1, Knod
                 ij = (j-1) * Knod + i
!                 BndrySrc(ij,el) = BndrySrc(ij,el) + NodesGradRadau(j,0) * BC_Values(i,-elemID(1,el)) * NormAS(i,el)
                 BndrySrc(i,ij,bel) = BndrySrc(i,ij,bel) + NodesGradRadau(j,0) * NormAS(i,el)
                 DO m = 1, Knod
                   im = (m-1) * Knod + i

!                   tmpx = 0.d0
!                   tmpy = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(j,0)
!                   DO l = 1, Knod
!                     tmpx = tmpx + NeuMatrix(i,l) * BC_Values(l,-elemID(1,el)) * NormAS(l,el)
!                     BndrySrc(ij,el) = BndrySrc(ij,el) + tmpy * NeuMatrix(m,l) * BC_Values(l,-elemID(1,el)) * NormAS(l,el)
!                   ENDDO
!                   BndrySrc(ij,el) = BndrySrc(ij,el) - tmpx * SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,0)
                   tmpx = SolnAyy(i,m,el) * SNGLB_2xSBLBdotNGR(m,j) * NodesGradRadau(m,0)
                   tmpy = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i) * NodesGradRadau(j,0)
                   DO l = 1, Knod
                     BndrySrc(l,ij,bel) = BndrySrc(l,ij,bel) + (tmpy * NeuMatrix(m,l) - tmpx * NeuMatrix(i,l)) * NormAS(l,el)
                   ENDDO

                   tmpx = 0.d0
                   DO l = 1, Knod
                     tmpx = tmpx + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SBLB_i0_NGR_j0(m,l)
                   ENDDO
                   LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) - tmpx

                   tmpy =0.d0
                   DO l = 1, Knod
                     tmpy = tmpy + SolnAxy(l,j,el) * SNGLB_2xSBLBdotNGR(l,i) * NeuMatrix(l,m)
                   ENDDO
                   DO l = 1, Knod
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) - 2.d0 * tmpy * FaceAS(m,el) * &
                                                                            (SBGLB_SBLB_i0_NGR_j0(l,j) + gL * SBLB_i0_NGR_j0(l,j))
                   ENDDO

                   DO l = 1, Knod
                     tmpx = 0.d0
                     DO p = 1, Knod
                       tmpx = tmpx + SolnAyy(i,p,el) * SNGLB_2xSBLBdotNGR(p,j) * &
                                                                            (SBGLB_SBLB_i0_NGR_j0(l,p) + gL * SBLB_i0_NGR_j0(l,p))
                     ENDDO
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) + 2.d0 * tmpx * FaceAS(m,el) * NeuMatrix(i,m)
                   ENDDO

                   tmpy = SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i)
                   DO l = 1, Knod
                     ml = (l-1) * Knod + m
                     LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) + tmpy * SBLB_i0_NGR_j0(l,j)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ENDIF

         ELSE

           DO j = 1, Knod
             Acoef(j) = FaceAN(j,elemID(1,el))
             Bcoef(j) = 0.5d0 * (FaceBS(j,el) + FaceBN(j,elemID(1,el)))
           ENDDO

           DO i = 1, Knod
             DO j = 1, Knod
               ij = (j-1) * Knod + i
               DO m = 1, Knod
                 im = (m-1) * Knod + i
                 tmpy = 0.d0
                 DO l = 1, Knod
                   tmpy = tmpy + SolnAyy(i,l,el) * SNGLB_2xSBLBdotNGR(l,j) * SBLB_i1_NGR_j0(m,l)
                 ENDDO
                 LaplacianSouth(im,ij,el) = LaplacianSouth(im,ij,el) + tmpy + &
                                               Acoef(i) * SBGLB_SBLB_i1_NGR_j0(m,j) - &
                                               gL * FaceAS(i,el) * SBLB_i1_NGR_j0(m,j)
                 LaplacianCenter(im,ij,el) = LaplacianCenter(im,ij,el) + ( &
                                               FaceAS(i,el) * SBGLB_SBLB_i0_NGR_j0(m,j) + &
                                               gL * Acoef(i) * SBLB_i0_NGR_j0(m,j) )
                 tmp = Bcoef(i) * SolnNodesGradLgrangeBasis(m,i)
                 tmpx = tmp + SolnAxy(m,j,el) * SNGLB_2xSBLBdotNGR(m,i)
                 DO l = 1, Knod
                   ml = (l-1) * Knod + m
                   LaplacianSouth(ml,ij,el) = LaplacianSouth(ml,ij,el) - tmpx * SBLB_i1_NGR_j0(l,j) 
                   LaplacianCenter(ml,ij,el) = LaplacianCenter(ml,ij,el) - tmp * SBLB_i0_NGR_j0(l,j)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

         ENDIF

         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
!             BndrySrc(1:Knod,ij,el) = BndrySrc(1:Knod,ij,el) !* SolnJac(i,j,el)
             SJacs(ij,el) = 1.d0 / SolnJac(i,j,el)
             DO m = 1, Knod
               mm1 = (m-1) * Knod
               DO l = 1, Knod
                 lm = mm1 + l
                 LaplacianCenter(lm,ij,el) = LaplacianCenter(lm,ij,el) !* SolnJac(i,j,el)
                 LaplacianEast(lm,ij,el) = LaplacianEast(lm,ij,el) !* SolnJac(i,j,el)
                 LaplacianWest(lm,ij,el) = LaplacianWest(lm,ij,el) !* SolnJac(i,j,el)
                 LaplacianNorth(lm,ij,el) = LaplacianNorth(lm,ij,el) !* SolnJac(i,j,el)
                 LaplacianSouth(lm,ij,el) = LaplacianSouth(lm,ij,el) !* SolnJac(i,j,el)
               ENDDO
             ENDDO
           ENDDO
         ENDDO

       ENDDO

!      Everything above here is used in conjunction with LaplacianP to get fast diffusion; for now, below is to solve
!      the Poisson equation, but it needs to be modified to replace Laplacian* variables into matrix A
       nrows = Ksq * Nel
       nnz = K4 * Nel
       DO el = 1, Nel
         IF (elemID(1,el) > 0) nnz = nnz + K4
         IF (elemID(2,el) > 0) nnz = nnz + K4
         IF (elemID(3,el) > 0) nnz = nnz + K4
         IF (elemID(4,el) > 0) nnz = nnz + K4
       ENDDO

       Allocate (rowptr(nrows+1), colidx(nnz), values(nnz))

       rowptr(nrows+1) = nnz + 1

       nrows = 0
       nnz = 0
       DO el = 1, Nel
         colc = (el - 1) * Ksq
         cole = (elemID(2,el) - 1) * Ksq
         colw = (elemID(4,el) - 1) * Ksq
         coln = (elemID(3,el) - 1) * Ksq
         cols = (elemID(1,el) - 1) * Ksq
         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
             nrows = nrows + 1
             rowptr(nrows) = nnz + 1
             DO m = 1, Knod
               mm1 = (m-1) * Knod
               DO l = 1, Knod
                 lm = mm1 + l
                 nnz = nnz + 1
                 values(nnz) = LaplacianCenter(lm,ij,el)
                 colidx(nnz) = colc + lm
               ENDDO
             ENDDO
             IF (elemID(2,el) > 0) THEN
               DO m = 1, Knod
                 mm1 = (m-1) * Knod
                 DO l = 1, Knod
                   lm = mm1 + l
                   nnz = nnz + 1
                   values(nnz) = LaplacianEast(lm,ij,el)
                   colidx(nnz) = cole + lm
                 ENDDO
               ENDDO
             ENDIF
             IF (elemID(4,el) > 0) THEN
               DO m = 1, Knod
                 mm1 = (m-1) * Knod
                 DO l = 1, Knod
                   lm = mm1 + l
                   nnz = nnz + 1
                   values(nnz) = LaplacianWest(lm,ij,el)
                   colidx(nnz) = colw + lm
                 ENDDO
               ENDDO
             ENDIF
             IF (elemID(3,el) > 0) THEN
               DO m = 1, Knod
                 mm1 = (m-1) * Knod
                 DO l = 1, Knod
                   lm = mm1 + l
                   nnz = nnz + 1
                   values(nnz) = LaplacianNorth(lm,ij,el)
                   colidx(nnz) = coln + lm
                 ENDDO
               ENDDO
             ENDIF
             IF (elemID(1,el) > 0) THEN
               DO m = 1, Knod
                 mm1 = (m-1) * Knod
                 DO l = 1, Knod
                   lm = mm1 + l
                   nnz = nnz + 1
                   values(nnz) = LaplacianSouth(lm,ij,el)
                   colidx(nnz) = cols + lm
                 ENDDO
               ENDDO
             ENDIF
           ENDDO
         ENDDO
       ENDDO
           
       ! Filter out zero's
!       print *, 'Zero Count: ', count( values(1:nnz) == 0.0d0 ), nnz, rowptr(nrows+1)
       print *, 'Zero Count: ', count( ABS(values(1:nnz)) < 1.0d-14 ), nnz, rowptr(nrows+1)
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
       if (.true.) then
         ierr = APLLES_Setup_Matrix_CSR(nrows,rowptr,colidx,values,A_handle)
       else
         ierr = APLLES_Setup_Matrix_CSR(nrows,rowptr_filtered,colidx_filtered,values_filtered,A_handle)
       endif

       ierr = APLLES_Matrix_Copy_CSR_To_Other(A_handle, "CSR", A_handle)

       ierr = APLLES_Setup_Solver(A_handle,trim(Aplles_solver_name),S_handle)

       ierr = APLLES_Setup_Precon(A_handle,trim(Aplles_precon_name),P_handle)

       deAllocate (rowptr, colidx, values)
       deAllocate (LaplacianCenter, LaplacianEast, LaplacianWest, LaplacianNorth, LaplacianSouth)

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
       real*8, dimension(Knod, Knod, Nel) :: Vortin, psi
       real*8, dimension(Knod * Knod, Nel) :: temp

       integer el, bel, i, j, l, m, ij, lm, jm1, mm1, jx, jy, lx, ly
       integer colc, cole, colw, coln, cols, nnz, nrows, Ksq, K4, ierr
       real*8 residual, volume

       real*8, dimension(:), allocatable :: b, x

       type(APLLES_MatrixHandleType) :: A_handle
       type(APLLES_PreconHandleType) :: P_handle
       type(APLLES_SolverHandleType) :: S_handle

  
       nrows = Knod * Knod * Nel
       Allocate (x(nrows), b(nrows))

       temp = 0.d0
       DO bel = 1, NelB
         el = boundaryPointsElemID(bel)
         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
             temp(ij,el) = temp(ij,el) + DOT_PRODUCT( BndrySrc(1:Knod,ij,bel) , BC_Values(1:Knod,bel) )
           ENDDO
         ENDDO
       ENDDO

       residual = 0.d0
       volume = 0.d0
       DO el = 1, Nel
         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
             residual = residual + wgt(i) *wgt(j) * Vortin(i,j,el) * SJacs(ij,el)
             volume = volume + wgt(i) *wgt(j) * SJacs(ij,el)
           ENDDO
         ENDDO
       ENDDO
       print *,'Res, Vol, eps ', residual, volume, (1+residual)/volume

       x = 0.d0
       b = 0.d0

       nrows = 0
       DO el = 1, Nel
         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
             nrows = nrows + 1
!             b(nrows) = - Vortin(i,j,el) - BndrySrc(ij,el)
!             b(nrows) = - SJacs(ij,el) * Vortin(i,j,el) - BndrySrc(ij,el)
            if(abs((1-residual)/volume) > .1) then
              b(nrows) = - SJacs(ij,el) * Vortin(i,j,el) - temp(ij,el)
            else
              b(nrows) = - SJacs(ij,el) * Vortin(i,j,el) - temp(ij,el)
!              b(nrows) = - SJacs(ij,el) * (-(1+residual)/volume + Vortin(i,j,el)) - temp(ij,el)
            endif
           ENDDO
         ENDDO
       ENDDO

       ierr = APLLES_Solve(A_handle, x, b, S_handle, P_handle)

       nrows = 0
       DO el = 1, Nel
         DO j = 1, Knod
           DO i = 1, Knod
             nrows = nrows + 1
             psi(i,j,el) = x(nrows)
           ENDDO
         ENDDO
       ENDDO

       deAllocate (x, b)

END SUBROUTINE GetLaplacian


SUBROUTINE GetLaplacGrads(HuynhSolver_type, uin, getCurl)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type, getCurl
       real*8, dimension(Knod, Knod, Nel) :: uin

       integer i, j, ip, jx, jy, lx, ly
       real*8, dimension(Knod, 0:1, Nel) :: uBx,  uBy,  comx,  comy    ! u at the L/R (Bx) and S/N (By) boundaries of a mesh
                                                                       ! com(mon) values at the interface along x and y
       real*8, dimension(Knod, 0:1, Nel):: dxuBx, dyuBy                ! derivatives of uBx/uBy along x (dx) and y (dy) directions
       real*8, dimension(Knod, 0:1, Nel):: dcomx, dcomy                ! derivatives of common values

       real*8, dimension(Knod) :: crossDuB                             ! cross derivatives using common values com
       real*8, dimension(Knod) :: Jac, dx_dxsi, dx_deta, dy_dxsi, dy_deta  ! coord transformation stuff
       real*8, dimension(Knod, Knod) :: sdx_dxsi, sdx_deta, sdy_dxsi, sdy_deta, sJac, f_tilda, g_tilda
       real*8, dimension(Knod, 0:1) :: f_tildaB, g_tildaB
       real*8 NeuMatrix(Knod, Knod), NeuRHS(Knod), Acoef(Knod), Bcoef(Knod) ! small dense matrix to obtain comx (per element) for a given Neumann BC
       real*8 gLprime, sdu_dxsi, sdu_deta

       real*8 tmpx, tmpy
       real*8, dimension(Lnod, Lnod) :: xloc, yloc

       real*8 theta2,theta1,thetam,thetad,rad2,rad1,radm,radd

       comx=0.d0
       comy=0.d0
       dcomx=0.d0
       dcomy=0.d0
       ubx=0.d0
       uby=0.d0
       dxubx=0.d0
       dyuby=0.d0

       ! Extrapolate the unknown, uin, and its derivative to the mesh boundaries using Lagrange polynomials of order Knod-1
       DO i = 1, Nel
         DO j = 1, Knod
           uBx(j,0,i) = dot_product(uin(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,0))            ! Left of mesh  - value
           uBx(j,1,i) = dot_product(uin(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,1))            ! Right of mesh - value
           uBy(j,0,i) = dot_product(uin(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,0))            ! South of mesh  - value
           uBy(j,1,i) = dot_product(uin(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,1))            ! North of mesh  - value

           dxuBx(j,0,i) = dot_product(uin(1:Knod,j,i), SolnBndryGradLgrangeBasis(1:Knod,0))      ! Left of mesh  - derivative along x
           dxuBx(j,1,i) = dot_product(uin(1:Knod,j,i), SolnBndryGradLgrangeBasis(1:Knod,1))      ! Right of mesh - derivative along x
           dyuBy(j,0,i) = dot_product(uin(j,1:Knod,i), SolnBndryGradLgrangeBasis(1:Knod,0))      ! South of mesh  - derivative along y
           dyuBy(j,1,i) = dot_product(uin(j,1:Knod,i), SolnBndryGradLgrangeBasis(1:Knod,1))      ! North of mesh - derivative along y
         ENDDO
       ENDDO

!       IF (HuynhSolver_type .eq. 2) THEN
       IF (.true.) THEN

         gLprime = 0.5d0*Knod*Knod
!         gLprime = 0.5d0*Knod*(Knod+1)

         ! Get the common values of uin at the mesh interfaces (we use simple averaging in this case)
         ! We definitely need a more efficient strategy - right now we're saving com-s twice, and can probably save on d?uB? as well
         ! We probably can and should write a more compact method so that we don't repeat the same stuff for NEWS (which becomes messier for NEWSBT in 3D)
         DO i = 1, Nel

           DO jy = 1, Lnod
             DO jx = 1, Lnod
               xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
               yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
             ENDDO
           ENDDO

           ! Right face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             dx_dxsi(1:Knod) = -thetad * (radm + sps(1:Knod)*radd) * sin(theta2)
             dy_dxsi(1:Knod) =  thetad * (radm + sps(1:Knod)*radd) * cos(theta2)
             dx_deta(1:Knod) =  radd * cos(theta2)
             dy_deta(1:Knod) =  radd * sin(theta2)
 endif

           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(2,i) .gt. 0) THEN

             ! Right Comx: Average of uBx(i,right) + uBx(i+1,left)
             comx(1:Knod,1,i) = 0.5d0*(uBx(1:Knod,1,i) + uBx(1:Knod,0,elemID(2,i)))

             ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
             dxuBx(1:Knod,1,i) = dxuBx(1:Knod,1,i) + (comx(1:Knod,1,i) - uBx(1:Knod,1,i))*gLprime   ! Right corrected derivative (7.17a) in paper

             ! Get the corrected values of grad(u) along the mesh interface;
             DO j = 1, Knod
               crossDuB(j) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             ! Get the interface f
             dxuBx(1:Knod,1,i) = (dxuBx(1:Knod,1,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod)

           ELSEIF (elemID(2,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(2,i)) .eq. DirichletBC) THEN

               comx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i))

               ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
               dxuBx(1:Knod,1,i) = dxuBx(1:Knod,1,i) + (comx(1:Knod,1,i) - uBx(1:Knod,1,i))*gLprime   ! Right corrected derivative (7.17a) in paper

               ! Get the corrected values of grad(u) along the mesh interface;
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               ! Get the interface f
               dxuBx(1:Knod,1,i) = (dxuBx(1:Knod,1,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod)

               dcomx(1:Knod,1,i) = dxuBx(1:Knod,1,i)

               Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
               VelocJump(1:Knod,-elemID(2,i)) = dcomx(1:Knod,1,i) / Sqrt(Acoef(1:Knod))

             ELSEIF (BC_Switch(-elemID(2,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomx(1:Knod,1,i) - Acoef(1:Knod)*(dxuBx(1:Knod,1,i) - uBx(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
!                 NeuMatrix(j,1:Knod) = -Bcoef(1:Knod) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,1,i), NeuRHS)

!               VelocJump(1:Knod,-elemID(2,i)) = comx(1:Knod,1,i)
               DO j = 1, Knod
                 VelocJump(j,-elemID(2,i)) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j)) / Sqrt(Acoef(j))
               ENDDO

             ENDIF

           ENDIF

           ! Left face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             dx_dxsi(1:Knod) = -thetad * (radm + sps(1:Knod)*radd) * sin(theta1)
             dy_dxsi(1:Knod) =  thetad * (radm + sps(1:Knod)*radd) * cos(theta1)
             dx_deta(1:Knod) =  radd * cos(theta1)
             dy_deta(1:Knod) =  radd * sin(theta1)
 endif

           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(4,i) .gt. 0) THEN

             ! Saving Left Comx is unnecessary ;I'm doing it here cuz it's easier to see the approach at this development stage
             comx(1:Knod,0,i) = 0.5d0*(uBx(1:Knod,0,i) + uBx(1:Knod,1,elemID(4,i))) 
             
             dxuBx(1:Knod,0,i) = dxuBx(1:Knod,0,i) - (comx(1:Knod,0,i) - uBx(1:Knod,0,i))*gLprime   ! Left corrected derivative (7.17b) in paper

             DO j = 1, Knod
               crossDuB(j) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             dxuBx(1:Knod,0,i) = (dxuBx(1:Knod,0,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod)       
            
           ELSEIF (elemID(4,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(4,i)) .eq. DirichletBC) THEN

               comx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i))

               ! Left corrected derivative (7.17b) in paper             
               dxuBx(1:Knod,0,i) = dxuBx(1:Knod,0,i) - (comx(1:Knod,0,i) - uBx(1:Knod,0,i))*gLprime

               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               dxuBx(1:Knod,0,i) = (dxuBx(1:Knod,0,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod)       

               dcomx(1:Knod,0,i) = dxuBx(1:Knod,0,i)

               Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
               VelocJump(1:Knod,-elemID(4,i)) = dcomx(1:Knod,0,i) / Sqrt(Acoef(1:Knod))
            
             ELSEIF (BC_Switch(-elemID(4,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomx(1:Knod,0,i) - Acoef(1:Knod)*(dxuBx(1:Knod,0,i) + uBx(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
!                 NeuMatrix(j,1:Knod) = -Bcoef(1:Knod) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,0,i), NeuRHS)

!               VelocJump(1:Knod,-elemID(4,i)) = comx(1:Knod,0,i)
               DO j = 1, Knod
                 VelocJump(j,-elemID(4,i)) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j)) / Sqrt(Acoef(j))
               ENDDO

             ENDIF

           ENDIF

           ! North face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             dx_dxsi(1:Knod) = -thetad * (rad2) * sin(thetam + sps(1:Knod)*thetad)
             dy_dxsi(1:Knod) =  thetad * (rad2) * cos(thetam + sps(1:Knod)*thetad)
             dx_deta(1:Knod) =  radd * cos(thetam + sps(1:Knod)*thetad)
             dy_deta(1:Knod) =  radd * sin(thetam + sps(1:Knod)*thetad)
 endif

           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(3,i) .gt. 0) THEN

             ! North Comy: Average of uBy(i,north) + uBy(i+1,south)
             comy(1:Knod,1,i) = 0.5d0*(uBy(1:Knod,1,i) + uBy(1:Knod,0,elemID(3,i)))
             
             ! North corrected derivative (7.17a) in paper
             dyuBy(1:Knod,1,i) = dyuBy(1:Knod,1,i) + (comy(1:Knod,1,i) - uBy(1:Knod,1,i))*gLprime

             DO j = 1, Knod
               crossDuB(j) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             dyuBy(1:Knod,1,i) = (dyuBy(1:Knod,1,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod)            

           ELSEIF (elemID(3,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(3,i)) .eq. DirichletBC) THEN

               comy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i))

               ! North corrected derivative (7.17a) in paper             
               dyuBy(1:Knod,1,i) = dyuBy(1:Knod,1,i) + (comy(1:Knod,1,i) - uBy(1:Knod,1,i))*gLprime

               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               dyuBy(1:Knod,1,i) = (dyuBy(1:Knod,1,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod)            

               dcomy(1:Knod,1,i) = dyuBy(1:Knod,1,i)

               Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
               VelocJump(1:Knod,-elemID(3,i)) = dcomy(1:Knod,1,i) / Sqrt(Acoef(1:Knod))

             ELSEIF (BC_Switch(-elemID(3,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomy(1:Knod,1,i) - Acoef(1:Knod)*(dyuBy(1:Knod,1,i) - uBy(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,1,i), NeuRHS)

!               VelocJump(1:Knod,-elemID(3,i)) = comy(1:Knod,1,i)
               DO j = 1, Knod
                 VelocJump(j,-elemID(3,i)) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j)) / Sqrt(Acoef(j))
                 IF (getCurl .eq. 1) VelocJump(j,-elemID(3,i)) = - VelocJump(j,-elemID(3,i))
               ENDDO

             ENDIF

           ENDIF

           ! South face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             dx_dxsi(1:Knod) = -thetad * (rad1) * sin(thetam + sps(1:Knod)*thetad)
             dy_dxsi(1:Knod) =  thetad * (rad1) * cos(thetam + sps(1:Knod)*thetad)
             dx_deta(1:Knod) =  radd * cos(thetam + sps(1:Knod)*thetad)
             dy_deta(1:Knod) =  radd * sin(thetam + sps(1:Knod)*thetad)
 endif

           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(1,i) .gt. 0) THEN

             ! Saving South Comy is unnecessary
             comy(1:Knod,0,i) = 0.5d0*(uBy(1:Knod,0,i) + uBy(1:Knod,1,elemID(1,i)))

             ! South corrected derivative (7.17b) in paper
             dyuBy(1:Knod,0,i) = dyuBy(1:Knod,0,i) - (comy(1:Knod,0,i) - uBy(1:Knod,0,i))*gLprime

             DO j = 1, Knod
               crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             dyuBy(1:Knod,0,i) = (dyuBy(1:Knod,0,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod) 

           ELSEIF (elemID(1,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(1,i)) .eq. DirichletBC) THEN

               comy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i))

               ! South corrected derivative (7.17b) in paper
               dyuBy(1:Knod,0,i) = dyuBy(1:Knod,0,i) - (comy(1:Knod,0,i) - uBy(1:Knod,0,i))*gLprime

               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               dyuBy(1:Knod,0,i) = (dyuBy(1:Knod,0,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod) 

               dcomy(1:Knod,0,i) = dyuBy(1:Knod,0,i)

               Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
               VelocJump(1:Knod,-elemID(1,i)) = dcomy(1:Knod,0,i) / Sqrt(Acoef(1:Knod))

             ELSEIF (BC_Switch(-elemID(1,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomy(1:Knod,0,i) - Acoef(1:Knod)*(dyuBy(1:Knod,0,i) + uBy(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,0,i), NeuRHS)

!               VelocJump(1:Knod,-elemID(1,i)) = comy(1:Knod,0,i)
               DO j = 1, Knod
                 VelocJump(j,-elemID(1,i)) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j)) / Sqrt(Acoef(j))
                 IF (getCurl .eq. 1) VelocJump(j,-elemID(1,i)) = - VelocJump(j,-elemID(1,i))
               ENDDO

             ENDIF

           ENDIF

         ENDDO

         DO i = 1, Nel
           ! Right face
           IF (elemID(2,i) .gt. 0) dcomx(1:Knod,1,i) = 0.5d0*(dxuBx(1:Knod,1,i) + dxuBx(1:Knod,0,elemID(2,i)))   ! Right dcomx: Average of dxuBx(i,right) + dxuBx(i+1,left)
           ! Left face
           IF (elemID(4,i) .gt. 0) dcomx(1:Knod,0,i) = 0.5d0*(dxuBx(1:Knod,0,i) + dxuBx(1:Knod,1,elemID(4,i)))   ! Left dcomx
           ! North face
           IF (elemID(3,i) .gt. 0) dcomy(1:Knod,1,i) = 0.5d0*(dyuBy(1:Knod,1,i) + dyuBy(1:Knod,0,elemID(3,i)))   ! North dcomy: Average of dyuBy(i,north) + dyuBy(i+1,south)
           ! South face
           IF (elemID(1,i) .gt. 0) dcomy(1:Knod,0,i) = 0.5d0*(dyuBy(1:Knod,0,i) + dyuBy(1:Knod,1,elemID(1,i)))   ! South dcomy
         ENDDO

!       ELSEIF (HuynhSolver_type .eq. 11) THEN
       ELSE

         gLprime = 0.5*Knod*(Knod+1)

         DO i = 1, Nel

           DO jy = 1, Lnod
             DO jx = 1, Lnod
               xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
               yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
             ENDDO
           ENDDO

           ! Right Face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
           Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
           Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

           ! Boundary face
           IF (elemID(2,i) .lt. 0) THEN
             ! Dirichlet BC
             IF (BC_Switch(-elemID(2,i)) .eq. DirichletBC) THEN
               comx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i))

               ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
               dxuBx(1:Knod,1,i) = dxuBx(1:Knod,1,i) + (comx(1:Knod,1,i) - uBx(1:Knod,1,i))*gLprime   ! Right corrected derivative (7.17a) in paper

               ! Get the corrected values of grad(u) along the mesh interface;
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               ! Get the interface f
               dcomx(1:Knod,1,i) = (Acoef(1:Knod)*dxuBx(1:Knod,1,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod)

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(2,i)) .eq. NeumannBC) THEN
               dcomx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomx(1:Knod,1,i) - Acoef(1:Knod)*(dxuBx(1:Knod,1,i) - uBx(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,1,i), NeuRHS)

             ENDIF
           ! Internal right face commoned with the left face of its right neighbor
           ELSE
             ! Right Face Stuff
             Acoef(1:Knod) = Acoef(1:Knod) / Jac(1:Knod)
             Bcoef(1:Knod) = Bcoef(1:Knod) / Jac(1:Knod)
             NeuRHS(1:Knod) = - Acoef(1:Knod) * (dxuBx(1:Knod,1,i) - uBx(1:Knod,1,i)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! Left Face Stuff from the Right neighbor
             ip = elemID(2,i)
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! grad at x-dir                   no grad along y-dir              x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! no grad at x-dir            grad along y-dir                     x/y-coord of geom
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) / Jac(1:Knod)
             Bcoef(1:Knod) = (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)) / Jac(1:Knod)

             NeuRHS(1:Knod) =  NeuRHS(1:Knod) + Acoef(1:Knod)*(dxuBx(1:Knod,0,ip) + uBx(1:Knod,0,ip)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = NeuMatrix(j,1:Knod) + Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! Get the common face value of the unknown
             CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,1,i), NeuRHS)
             comx(1:Knod,0,ip) = comx(1:Knod,1,i)

             ! Get the common f~ value at the face
             DO j = 1, Knod
               crossDuB(j) = dot_product(comx(1:Knod,0,ip), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO
             dcomx(1:Knod,0,ip) = Acoef(1:Knod) * (dxuBx(1:Knod,0,ip) - (comx(1:Knod,0,ip) - uBx(1:Knod,0,ip))*gLprime) -  &
                                  Bcoef(1:Knod) * crossDuB(1:Knod)
             dcomx(1:Knod,1,i) = dcomx(1:Knod,0,ip)
           ENDIF

           ! Left Face
           IF (elemID(4,i) .lt. 0) THEN
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
             Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

             ! Dirichlet BC
             IF (BC_Switch(-elemID(4,i)) .eq. DirichletBC) THEN
               comx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i))
             
               dxuBx(1:Knod,0,i) = dxuBx(1:Knod,0,i) - (comx(1:Knod,0,i) - uBx(1:Knod,0,i))*gLprime   ! Left corrected derivative (7.17b) in paper
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO
               dcomx(1:Knod,0,i) = (Acoef(1:Knod)*dxuBx(1:Knod,0,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod)       

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(4,i)) .eq. NeumannBC) THEN
               dcomx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomx(1:Knod,0,i) - Acoef(1:Knod)*(dxuBx(1:Knod,0,i) + uBx(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,0,i), NeuRHS)

             ENDIF
           ENDIF

           ! North Face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
           Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
           Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

           ! Boundary face
           IF (elemID(3,i) .lt. 0) THEN
             ! Dirichlet BC
             IF (BC_Switch(-elemID(3,i)) .eq. DirichletBC) THEN
               comy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i))

               ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
               dyuBy(1:Knod,1,i) = dyuBy(1:Knod,1,i) + (comy(1:Knod,1,i) - uBy(1:Knod,1,i))*gLprime   ! North corrected derivative (7.17a) in paper

               ! Get the corrected values of grad(u) along the mesh interface;
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               ! Get the interface g
               dcomy(1:Knod,1,i) = (Acoef(1:Knod)*dyuBy(1:Knod,1,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod)

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(3,i)) .eq. NeumannBC) THEN
               dcomy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomy(1:Knod,1,i) - Acoef(1:Knod)*(dyuBy(1:Knod,1,i) - uBy(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,1,i), NeuRHS)

             ENDIF
           ! Internal north face commoned with the south face of its north neighbor
           ELSE
             ! North Face Stuff
             Acoef(1:Knod) = Acoef(1:Knod) / Jac(1:Knod)
             Bcoef(1:Knod) = Bcoef(1:Knod) / Jac(1:Knod)
             NeuRHS(1:Knod) = - Acoef(1:Knod)*(dyuBy(1:Knod,1,i) - uBy(1:Knod,1,i)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! South Face Stuff from the North neighbor
             ip = elemID(3,i)
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! grad at y-dir                   no grad along x-dir              x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! no grad at y-dir            grad along x-dir                     x/y-coord of geom
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) / Jac(1:Knod)
             Bcoef(1:Knod) = (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)) / Jac(1:Knod)

             NeuRHS(1:Knod) =  NeuRHS(1:Knod) + Acoef(1:Knod)*(dyuBy(1:Knod,0,ip) + uBy(1:Knod,0,ip)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = NeuMatrix(j,1:Knod) + Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! Get the common face value of the unknown
             CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,1,i), NeuRHS)
             comy(1:Knod,0,ip) = comy(1:Knod,1,i)

             ! Get the common g~ value at the face
             DO j = 1, Knod
               crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO
             dcomy(1:Knod,0,ip) = Acoef(1:Knod) * (dyuBy(1:Knod,0,ip) - (comy(1:Knod,0,ip) - uBy(1:Knod,0,ip))*gLprime) -  &
                                  Bcoef(1:Knod) * crossDuB(1:Knod)
             dcomy(1:Knod,1,i) = dcomy(1:Knod,0,ip)

           ENDIF

           ! South Face
           IF (elemID(1,i) .lt. 0) THEN
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
             Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

             ! Dirichlet BC
             IF (BC_Switch(-elemID(1,i)) .eq. DirichletBC) THEN
               comy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i))

               dyuBy(1:Knod,0,i) = dyuBy(1:Knod,0,i) - (comy(1:Knod,0,i) - uBy(1:Knod,0,i))*gLprime   ! South corrected derivative (7.17b) in paper
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO
               dcomy(1:Knod,0,i) = (Acoef(1:Knod)*dyuBy(1:Knod,0,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod) 

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(1,i)) .eq. NeumannBC) THEN
               dcomy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomy(1:Knod,0,i) - Acoef(1:Knod)*(dyuBy(1:Knod,0,i) + uBy(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,0,i), NeuRHS)

             ENDIF
           ENDIF
         ENDDO

       ENDIF

       DO i = 1, Nel

         DO jy = 1, Lnod
           DO jx = 1, Lnod
             xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
             yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
           ENDDO
         ENDDO

         ! Get first derivatives of the variable
         DO jy = 1, Knod
           DO jx = 1, Knod
             ! Get geometry stuff at all solution nodes (per element)
             sdx_dxsi(jx,jy) = 0.d0
             sdx_deta(jx,jy) = 0.d0
             sdy_dxsi(jx,jy) = 0.d0
             sdy_deta(jx,jy) = 0.d0
             DO ly = 1, Lnod
               DO lx = 1, Lnod
                 sdx_dxsi(jx,jy) = sdx_dxsi(jx,jy) + &
                            ! grad at x-dir                    no grad at y-dir               x/y-coord of geom
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_dxsi(jx,jy) = sdy_dxsi(jx,jy) + &
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
                 sdx_deta(jx,jy) = sdx_deta(jx,jy) + &
                            ! no grad at x-dir             grad at y-dir                      x/y-coord of geom
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_deta(jx,jy) = sdy_deta(jx,jy) + &
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
               ENDDO
             ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             sdx_dxsi(jx,jy) = -thetad * (radm + sps(jy)*radd) * sin(thetam + sps(jx)*thetad)
             sdy_dxsi(jx,jy) =  thetad * (radm + sps(jy)*radd) * cos(thetam + sps(jx)*thetad)
             sdx_deta(jx,jy) =  radd * cos(thetam + sps(jx)*thetad)
             sdy_deta(jx,jy) =  radd * sin(thetam + sps(jx)*thetad)
 endif

             sJac(jx,jy) = sdx_dxsi(jx,jy) * sdy_deta(jx,jy) - sdx_deta(jx,jy) * sdy_dxsi(jx,jy)

             ! Get grads of unknown uin along xsi (for eta=const) and along eta (for xsi=const
             sdu_dxsi = dot_product(uin(1:Knod,jy,i), SolnNodesGradLgrangeBasis(1:Knod,jx)) + &
                        dot_product(comx(jy,0:1,i) - uBx(jy,0:1,i), NodesGradRadau(jx,0:1))
             sdu_deta = dot_product(uin(jx,1:Knod,i), SolnNodesGradLgrangeBasis(1:Knod,jy)) + &
                        dot_product(comy(jx,0:1,i) - uBy(jx,0:1,i), NodesGradRadau(jy,0:1))

             ! Get f~ and g~ as per Huynh's paper
             f_tilda(jx,jy) = (sdu_dxsi * (sdx_deta(jx,jy)*sdx_deta(jx,jy) + sdy_deta(jx,jy)*sdy_deta(jx,jy)) - &
                               sdu_deta * (sdx_dxsi(jx,jy)*sdx_deta(jx,jy) + sdy_dxsi(jx,jy)*sdy_deta(jx,jy))) / sJac(jx,jy)
             g_tilda(jx,jy) = (sdu_deta * (sdx_dxsi(jx,jy)*sdx_dxsi(jx,jy) + sdy_dxsi(jx,jy)*sdy_dxsi(jx,jy)) - &
                               sdu_dxsi * (sdx_dxsi(jx,jy)*sdx_deta(jx,jy) + sdy_dxsi(jx,jy)*sdy_deta(jx,jy))) / sJac(jx,jy)
           ENDDO
         ENDDO

         ! Get ?_tildas at the mesh boundaries
         DO j = 1, Knod
           f_tildaB(j,0) = dot_product(f_tilda(1:Knod,j), SolnBndryLgrangeBasis(1:Knod,0))      ! Left of mesh - derivative
           f_tildaB(j,1) = dot_product(f_tilda(1:Knod,j), SolnBndryLgrangeBasis(1:Knod,1))      ! Right of mesh - derivative
           g_tildaB(j,0) = dot_product(g_tilda(j,1:Knod), SolnBndryLgrangeBasis(1:Knod,0))      ! South of mesh - derivative
           g_tildaB(j,1) = dot_product(g_tilda(j,1:Knod), SolnBndryLgrangeBasis(1:Knod,1))      ! North of mesh - derivative
         ENDDO

         ! Get potential velocity
         IF (getCurl .eq. 0) THEN
           DO jy = 1, Knod
             DO jx = 1, Knod
               tmpx = f_tilda(jx,jy) + dot_product(dcomx(jy,0:1,i) - f_tildaB(jy,0:1), NodesRadau(jx,0:1))
               tmpy = g_tilda(jx,jy) + dot_product(dcomy(jx,0:1,i) - g_tildaB(jx,0:1), NodesRadau(jy,0:1))
               Uvel(jx,jy,i) = Uvel(jx,jy,i) + (tmpx*sdx_dxsi(jx,jy) + tmpy*sdx_deta(jx,jy)) / sJac(jx,jy)
               Vvel(jx,jy,i) = Vvel(jx,jy,i) + (tmpx*sdy_dxsi(jx,jy) + tmpy*sdy_deta(jx,jy)) / sJac(jx,jy)
             ENDDO
           ENDDO
         ! Get vortical velocity
         ELSEIF (getCurl .eq. 1) THEN
           DO jy = 1, Knod
             DO jx = 1, Knod
               tmpx = f_tilda(jx,jy) + dot_product(dcomx(jy,0:1,i) - f_tildaB(jy,0:1), NodesRadau(jx,0:1))
               tmpy = g_tilda(jx,jy) + dot_product(dcomy(jx,0:1,i) - g_tildaB(jx,0:1), NodesRadau(jy,0:1))
               Uvel(jx,jy,i) =  (tmpx*sdy_dxsi(jx,jy) + tmpy*sdy_deta(jx,jy)) / sJac(jx,jy)
               Vvel(jx,jy,i) = -(tmpx*sdx_dxsi(jx,jy) + tmpy*sdx_deta(jx,jy)) / sJac(jx,jy)
             ENDDO
           ENDDO
         ELSE
           print *,'ERROR: getCurl is either 0 (for potential velocity) or 1 (for vortical velocity)!"'
           stop
         ENDIF
       ENDDO

END SUBROUTINE GetLaplacGrads


SUBROUTINE GetLaplacianP(HuynhSolver_type, uin, ddu)
       USE params
       USE variables

       implicit NONE
       integer HuynhSolver_type
       real*8, dimension(Knod, Knod, Nel) :: uin, ddu

       integer el, i, j, l, m, ij, lm, jm1, mm1

       ddu = 0.d0

       DO el = 1, Nel
         DO j = 1, Knod
           jm1 = (j-1) * Knod
           DO i = 1, Knod
             ij = jm1 + i
!!!!             ddu(i,j,el) = ddu(i,j,el) + BndrySrc(ij,el)
             DO m = 1, Knod
               mm1 = (m-1) * Knod
               DO l = 1, Knod
                 lm = mm1 + l
! The Laplacian* variables need to be modified to use matrix A  (this is for diffusion solver)
!                 ddu(i,j,el) = ddu(i,j,el) + LaplacianCenter(lm,ij,el) * uin(l,m,el)
!                 IF (elemID(2,el) > 0) ddu(i,j,el) = ddu(i,j,el) + LaplacianEast(lm,ij,el) * uin(l,m,elemID(2,el))
!                 IF (elemID(4,el) > 0) ddu(i,j,el) = ddu(i,j,el) + LaplacianWest(lm,ij,el) * uin(l,m,elemID(4,el))
!                 IF (elemID(3,el) > 0) ddu(i,j,el) = ddu(i,j,el) + LaplacianNorth(lm,ij,el) * uin(l,m,elemID(3,el))
!                 IF (elemID(1,el) > 0) ddu(i,j,el) = ddu(i,j,el) + LaplacianSouth(lm,ij,el) * uin(l,m,elemID(1,el))
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO

END SUBROUTINE GetLaplacianP


SUBROUTINE GetDiffusedFlux(HuynhSolver_type, uin, ddu)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type
       real*8, dimension(Knod, Knod, Nel) :: uin, ddu

       integer i, j, ip, jx, jy, lx, ly
       real*8, dimension(Knod, 0:1, Nel) :: uBx,  uBy,  comx,  comy    ! u at the L/R (Bx) and S/N (By) boundaries of a mesh
                                                                       ! com(mon) values at the interface along x and y
       real*8, dimension(Knod, 0:1, Nel):: dxuBx, dyuBy                ! derivatives of uBx/uBy along x (dx) and y (dy) directions
       real*8, dimension(Knod, 0:1, Nel):: dcomx, dcomy                ! derivatives of common values

       real*8, dimension(Knod) :: crossDuB                             ! cross derivatives using common values com
       real*8, dimension(Knod) :: Jac, dx_dxsi, dx_deta, dy_dxsi, dy_deta  ! coord transformation stuff
       real*8, dimension(Knod, Knod) :: sJac, f_tilda, g_tilda
       real*8, dimension(Knod, 0:1) :: f_tildaB, g_tildaB
       real*8 gLprime, sdx_dxsi, sdx_deta, sdy_dxsi, sdy_deta, sdu_dxsi, sdu_deta
       real*8 NeuMatrix(Knod, Knod), NeuRHS(Knod), Acoef(Knod), Bcoef(Knod) ! small dense matrix to obtain comx (per element) for a given Neumann BC

       real*8, dimension(Lnod, Lnod) :: xloc, yloc

       comx=0.d0
       comy=0.d0
       dcomx=0.d0
       dcomy=0.d0
       ubx=0.d0
       uby=0.d0
       dxubx=0.d0
       dyuby=0.d0

       ! Extrapolate the unknown, uin, and its derivative to the mesh boundaries using Lagrange polynomials of order Knod-1
       DO i = 1, Nel
         DO j = 1, Knod
           uBx(j,0,i) = dot_product(uin(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,0))            ! Left of mesh  - value
           uBx(j,1,i) = dot_product(uin(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,1))            ! Right of mesh - value
           uBy(j,0,i) = dot_product(uin(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,0))            ! South of mesh  - value
           uBy(j,1,i) = dot_product(uin(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,1))            ! North of mesh  - value

           dxuBx(j,0,i) = dot_product(uin(1:Knod,j,i), SolnBndryGradLgrangeBasis(1:Knod,0))      ! Left of mesh  - derivative along x
           dxuBx(j,1,i) = dot_product(uin(1:Knod,j,i), SolnBndryGradLgrangeBasis(1:Knod,1))      ! Right of mesh - derivative along x
           dyuBy(j,0,i) = dot_product(uin(j,1:Knod,i), SolnBndryGradLgrangeBasis(1:Knod,0))      ! South of mesh  - derivative along y
           dyuBy(j,1,i) = dot_product(uin(j,1:Knod,i), SolnBndryGradLgrangeBasis(1:Knod,1))      ! North of mesh - derivative along y
         ENDDO
       ENDDO

       IF (HuynhSolver_type .eq. 2) THEN

         gLprime = 0.5d0*Knod*Knod

         ! Get the common values of uin at the mesh interfaces (we use simple averaging in this case)
         ! We definitely need a more efficient strategy - right now we're saving com-s twice, and can probably save on d?uB? as well
         ! We probably can and should write a more compact method so that we don't repeat the same stuff for NEWS (which becomes messier for NEWSBT in 3D)
         DO i = 1, Nel

           DO jy = 1, Lnod
             DO jx = 1, Lnod
               xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
               yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
             ENDDO
           ENDDO

           ! Right face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(2,i) .gt. 0) THEN

             ! Right Comx: Average of uBx(i,right) + uBx(i+1,left)
             comx(1:Knod,1,i) = 0.5d0*(uBx(1:Knod,1,i) + uBx(1:Knod,0,elemID(2,i)))

             ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
             dxuBx(1:Knod,1,i) = dxuBx(1:Knod,1,i) + (comx(1:Knod,1,i) - uBx(1:Knod,1,i))*gLprime   ! Right corrected derivative (7.17a) in paper

             ! Get the corrected values of grad(u) along the mesh interface;
             DO j = 1, Knod
               crossDuB(j) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             ! Get the interface f
             dxuBx(1:Knod,1,i) = (dxuBx(1:Knod,1,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod)

           ELSEIF (elemID(2,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(2,i)) .eq. DirichletBC) THEN

               comx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i))

               ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
               dxuBx(1:Knod,1,i) = dxuBx(1:Knod,1,i) + (comx(1:Knod,1,i) - uBx(1:Knod,1,i))*gLprime   ! Right corrected derivative (7.17a) in paper

               ! Get the corrected values of grad(u) along the mesh interface;
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               ! Get the interface f
               dxuBx(1:Knod,1,i) = (dxuBx(1:Knod,1,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod)

               dcomx(1:Knod,1,i) = dxuBx(1:Knod,1,i)

             ELSEIF (BC_Switch(-elemID(2,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*BC_Values(1:Knod,-elemID(2,i)) - &
                                Acoef(1:Knod)*(dxuBx(1:Knod,1,i) - uBx(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,1,i), NeuRHS)

             ENDIF

           ENDIF

           ! Left face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(4,i) .gt. 0) THEN

             ! Saving Left Comx is unnecessary ;I'm doing it here cuz it's easier to see the approach at this development stage
             comx(1:Knod,0,i) = 0.5d0*(uBx(1:Knod,0,i) + uBx(1:Knod,1,elemID(4,i))) 

             ! Left corrected derivative (7.17b) in paper             
             dxuBx(1:Knod,0,i) = dxuBx(1:Knod,0,i) - (comx(1:Knod,0,i) - uBx(1:Knod,0,i))*gLprime

             DO j = 1, Knod
               crossDuB(j) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             dxuBx(1:Knod,0,i) = (dxuBx(1:Knod,0,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod)       

           ELSEIF (elemID(4,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(4,i)) .eq. DirichletBC) THEN

               comx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i))

               ! Left corrected derivative (7.17b) in paper             
               dxuBx(1:Knod,0,i) = dxuBx(1:Knod,0,i) - (comx(1:Knod,0,i) - uBx(1:Knod,0,i))*gLprime

               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               dxuBx(1:Knod,0,i) = (dxuBx(1:Knod,0,i) * (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod)       

               dcomx(1:Knod,0,i) = dxuBx(1:Knod,0,i)

             ELSEIF (BC_Switch(-elemID(4,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*BC_Values(1:Knod,-elemID(4,i)) - &
                                Acoef(1:Knod)*(dxuBx(1:Knod,0,i) + uBx(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,0,i), NeuRHS)

             ENDIF

           ENDIF

           ! North face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(3,i) .gt. 0) THEN

             ! North Comy: Average of uBy(i,north) + uBy(i+1,south)
             comy(1:Knod,1,i) = 0.5d0*(uBy(1:Knod,1,i) + uBy(1:Knod,0,elemID(3,i)))

             ! North corrected derivative (7.17a) in paper             
             dyuBy(1:Knod,1,i) = dyuBy(1:Knod,1,i) + (comy(1:Knod,1,i) - uBy(1:Knod,1,i))*gLprime

             DO j = 1, Knod
               crossDuB(j) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             dyuBy(1:Knod,1,i) = (dyuBy(1:Knod,1,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod)            

           ELSEIF (elemID(3,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(3,i)) .eq. DirichletBC) THEN

               comy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i))
             
               ! North corrected derivative (7.17a) in paper
               dyuBy(1:Knod,1,i) = dyuBy(1:Knod,1,i) + (comy(1:Knod,1,i) - uBy(1:Knod,1,i))*gLprime

               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               dyuBy(1:Knod,1,i) = (dyuBy(1:Knod,1,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod)            

               dcomy(1:Knod,1,i) = dyuBy(1:Knod,1,i)

             ELSEIF (BC_Switch(-elemID(3,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*BC_Values(1:Knod,-elemID(3,i)) - &
                                Acoef(1:Knod)*(dyuBy(1:Knod,1,i) - uBy(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,1,i), NeuRHS)

             ENDIF

           ENDIF

           ! South face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)

           IF (elemID(1,i) .gt. 0) THEN

             ! Saving South Comy is unnecessary
             comy(1:Knod,0,i) = 0.5d0*(uBy(1:Knod,0,i) + uBy(1:Knod,1,elemID(1,i)))

             ! South corrected derivative (7.17b) in paper
             dyuBy(1:Knod,0,i) = dyuBy(1:Knod,0,i) - (comy(1:Knod,0,i) - uBy(1:Knod,0,i))*gLprime

             DO j = 1, Knod
               crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO

             dyuBy(1:Knod,0,i) = (dyuBy(1:Knod,0,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                   crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                 ) / Jac(1:Knod) 

           ELSEIF (elemID(1,i) .lt. 0) THEN

             IF (BC_Switch(-elemID(1,i)) .eq. DirichletBC) THEN

               comy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i))

               ! South corrected derivative (7.17b) in paper
               dyuBy(1:Knod,0,i) = dyuBy(1:Knod,0,i) - (comy(1:Knod,0,i) - uBy(1:Knod,0,i))*gLprime

               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               dyuBy(1:Knod,0,i) = (dyuBy(1:Knod,0,i) * (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) -   &
                                     crossDuB(1:Knod) * (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod))     &
                                   ) / Jac(1:Knod) 

               dcomy(1:Knod,0,i) = dyuBy(1:Knod,0,i)

             ELSEIF (BC_Switch(-elemID(1,i)) .eq. NeumannBC) THEN

               Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
               Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

               dcomy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*BC_Values(1:Knod,-elemID(1,i)) - &
                                Acoef(1:Knod)*(dyuBy(1:Knod,0,i) + uBy(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,0,i), NeuRHS)

             ENDIF

           ENDIF

         ENDDO

         DO i = 1, Nel
           ! Right face
           IF (elemID(2,i) .gt. 0) dcomx(1:Knod,1,i) = 0.5d0*(dxuBx(1:Knod,1,i) + dxuBx(1:Knod,0,elemID(2,i)))   ! Right dcomx: Average of dxuBx(i,right) + dxuBx(i+1,left)
           ! Left face
           IF (elemID(4,i) .gt. 0) dcomx(1:Knod,0,i) = 0.5d0*(dxuBx(1:Knod,0,i) + dxuBx(1:Knod,1,elemID(4,i)))   ! Left dcomx
           ! North face
           IF (elemID(3,i) .gt. 0) dcomy(1:Knod,1,i) = 0.5d0*(dyuBy(1:Knod,1,i) + dyuBy(1:Knod,0,elemID(3,i)))   ! North dcomy: Average of dyuBy(i,north) + dyuBy(i+1,south)
           ! South face
           IF (elemID(1,i) .gt. 0) dcomy(1:Knod,0,i) = 0.5d0*(dyuBy(1:Knod,0,i) + dyuBy(1:Knod,1,elemID(1,i)))   ! South dcomy
         ENDDO

       ELSEIF (HuynhSolver_type .eq. 11) THEN

         gLprime = 0.5*Knod*(Knod+1)

         DO i = 1, Nel

           DO jy = 1, Lnod
             DO jx = 1, Lnod
               xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
               yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
             ENDDO
           ENDDO

           ! Right Face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
           Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
           Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

           ! Boundary face
           IF (elemID(2,i) .lt. 0) THEN
             ! Dirichlet BC
             IF (BC_Switch(-elemID(2,i)) .eq. DirichletBC) THEN
               comx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i))

               ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
               dxuBx(1:Knod,1,i) = dxuBx(1:Knod,1,i) + (comx(1:Knod,1,i) - uBx(1:Knod,1,i))*gLprime   ! Right corrected derivative (7.17a) in paper

               ! Get the corrected values of grad(u) along the mesh interface;
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               ! Get the interface f
               dcomx(1:Knod,1,i) = (Acoef(1:Knod)*dxuBx(1:Knod,1,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod)

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(2,i)) .eq. NeumannBC) THEN
               dcomx(1:Knod,1,i) = BC_Values(1:Knod,-elemID(2,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomx(1:Knod,1,i) - Acoef(1:Knod)*(dxuBx(1:Knod,1,i) - uBx(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,1,i), NeuRHS)

               VelocJump(1:Knod,-elemID(2,i)) = comx(1:Knod,1,i)

             ENDIF
           ! Internal right face commoned with the left face of its right neighbor
           ELSE
             ! Right Face Stuff
             Acoef(1:Knod) = Acoef(1:Knod) / Jac(1:Knod)
             Bcoef(1:Knod) = Bcoef(1:Knod) / Jac(1:Knod)
             NeuRHS(1:Knod) = - Acoef(1:Knod) * (dxuBx(1:Knod,1,i) - uBx(1:Knod,1,i)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! Left Face Stuff from the Right neighbor
             ip = elemID(2,i)
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! grad at x-dir                   no grad along y-dir              x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! no grad at x-dir            grad along y-dir                     x/y-coord of geom
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = (dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)) / Jac(1:Knod)
             Bcoef(1:Knod) = (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)) / Jac(1:Knod)

             NeuRHS(1:Knod) =  NeuRHS(1:Knod) + Acoef(1:Knod)*(dxuBx(1:Knod,0,ip) + uBx(1:Knod,0,ip)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = NeuMatrix(j,1:Knod) + Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! Get the common face value of the unknown
             CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,1,i), NeuRHS)
             comx(1:Knod,0,ip) = comx(1:Knod,1,i)

             ! Get the common f~ value at the face
             DO j = 1, Knod
               crossDuB(j) = dot_product(comx(1:Knod,0,ip), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO
             dcomx(1:Knod,0,ip) = Acoef(1:Knod) * (dxuBx(1:Knod,0,ip) - (comx(1:Knod,0,ip) - uBx(1:Knod,0,ip))*gLprime) -  &
                                  Bcoef(1:Knod) * crossDuB(1:Knod)
             dcomx(1:Knod,1,i) = dcomx(1:Knod,0,ip)
           ENDIF

           ! Left Face
           IF (elemID(4,i) .lt. 0) THEN
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = dx_deta(1:Knod)*dx_deta(1:Knod) + dy_deta(1:Knod)*dy_deta(1:Knod)
             Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

             ! Dirichlet BC
             IF (BC_Switch(-elemID(4,i)) .eq. DirichletBC) THEN
               comx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i))
             
               dxuBx(1:Knod,0,i) = dxuBx(1:Knod,0,i) - (comx(1:Knod,0,i) - uBx(1:Knod,0,i))*gLprime   ! Left corrected derivative (7.17b) in paper
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comx(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO
               dcomx(1:Knod,0,i) = (Acoef(1:Knod)*dxuBx(1:Knod,0,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod)       

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(4,i)) .eq. NeumannBC) THEN
               dcomx(1:Knod,0,i) = BC_Values(1:Knod,-elemID(4,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomx(1:Knod,0,i) - Acoef(1:Knod)*(dxuBx(1:Knod,0,i) + uBx(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comx(1:Knod,0,i), NeuRHS)

               VelocJump(1:Knod,-elemID(4,i)) = comx(1:Knod,0,i)

             ENDIF
           ENDIF

           ! North Face
           dx_dxsi = 0.d0
           dx_deta = 0.d0
           dy_dxsi = 0.d0
           dy_deta = 0.d0
           DO jy = 1, Lnod
             DO jx = 1, Lnod
               dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                 ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                 GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                 ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
               dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                 GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             ENDDO
           ENDDO
           Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
           Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
           Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

           ! Boundary face
           IF (elemID(3,i) .lt. 0) THEN
             ! Dirichlet BC
             IF (BC_Switch(-elemID(3,i)) .eq. DirichletBC) THEN
               comy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i))

               ! Get the corrected values of grad(u) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG'=-gRB' (=-Knod^2/2)
               dyuBy(1:Knod,1,i) = dyuBy(1:Knod,1,i) + (comy(1:Knod,1,i) - uBy(1:Knod,1,i))*gLprime   ! North corrected derivative (7.17a) in paper

               ! Get the corrected values of grad(u) along the mesh interface;
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,1,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO

               ! Get the interface g
               dcomy(1:Knod,1,i) = (Acoef(1:Knod)*dyuBy(1:Knod,1,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod)

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(3,i)) .eq. NeumannBC) THEN
               dcomy(1:Knod,1,i) = BC_Values(1:Knod,-elemID(3,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomy(1:Knod,1,i) - Acoef(1:Knod)*(dyuBy(1:Knod,1,i) - uBy(1:Knod,1,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,1,i), NeuRHS)

               VelocJump(1:Knod,-elemID(3,i)) = comy(1:Knod,1,i)

             ENDIF
           ! Internal north face commoned with the south face of its north neighbor
           ELSE
             ! North Face Stuff
             Acoef(1:Knod) = Acoef(1:Knod) / Jac(1:Knod)
             Bcoef(1:Knod) = Bcoef(1:Knod) / Jac(1:Knod)
             NeuRHS(1:Knod) = - Acoef(1:Knod)*(dyuBy(1:Knod,1,i) - uBy(1:Knod,1,i)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! South Face Stuff from the North neighbor
             ip = elemID(3,i)
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! grad at y-dir                   no grad along x-dir              x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! no grad at y-dir            grad along x-dir                     x/y-coord of geom
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod)*xcoord(nodeID(t2f(jx,jy),ip))
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod)*ycoord(nodeID(t2f(jx,jy),ip))
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = (dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)) / Jac(1:Knod)
             Bcoef(1:Knod) = (dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)) / Jac(1:Knod)

             NeuRHS(1:Knod) =  NeuRHS(1:Knod) + Acoef(1:Knod)*(dyuBy(1:Knod,0,ip) + uBy(1:Knod,0,ip)*gLprime)
             DO j = 1, Knod
               NeuMatrix(j,1:Knod) = NeuMatrix(j,1:Knod) + Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
               NeuMatrix(j,j) = NeuMatrix(j,j) + Acoef(j)*gLprime
             ENDDO

             ! Get the common face value of the unknown
             CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,1,i), NeuRHS)
             comy(1:Knod,0,ip) = comy(1:Knod,1,i)

             ! Get the common g~ value at the face
             DO j = 1, Knod
               crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
             ENDDO
             dcomy(1:Knod,0,ip) = Acoef(1:Knod) * (dyuBy(1:Knod,0,ip) - (comy(1:Knod,0,ip) - uBy(1:Knod,0,ip))*gLprime) -  &
                                  Bcoef(1:Knod) * crossDuB(1:Knod)
             dcomy(1:Knod,1,i) = dcomy(1:Knod,0,ip)

           ENDIF

           ! South Face
           IF (elemID(1,i) .lt. 0) THEN
             dx_dxsi = 0.d0
             dx_deta = 0.d0
             dy_dxsi = 0.d0
             dy_deta = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 dx_deta(1:Knod) = dx_deta(1:Knod) + &
                                   ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
                 dy_deta(1:Knod) = dy_deta(1:Knod) + &
                                   GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
                 dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                                   ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
                 dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                                   GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
               ENDDO
             ENDDO
             Jac(1:Knod) = dx_dxsi(1:Knod) * dy_deta(1:Knod) - dx_deta(1:Knod) * dy_dxsi(1:Knod)
             Acoef(1:Knod) = dx_dxsi(1:Knod)*dx_dxsi(1:Knod) + dy_dxsi(1:Knod)*dy_dxsi(1:Knod)
             Bcoef(1:Knod) = dx_dxsi(1:Knod)*dx_deta(1:Knod) + dy_dxsi(1:Knod)*dy_deta(1:Knod)

             ! Dirichlet BC
             IF (BC_Switch(-elemID(1,i)) .eq. DirichletBC) THEN
               comy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i))

               dyuBy(1:Knod,0,i) = dyuBy(1:Knod,0,i) - (comy(1:Knod,0,i) - uBy(1:Knod,0,i))*gLprime   ! South corrected derivative (7.17b) in paper
               DO j = 1, Knod
                 crossDuB(j) = dot_product(comy(1:Knod,0,i), SolnNodesGradLgrangeBasis(1:Knod,j))
               ENDDO
               dcomy(1:Knod,0,i) = (Acoef(1:Knod)*dyuBy(1:Knod,0,i) - Bcoef(1:Knod)*crossDuB(1:Knod)) / Jac(1:Knod) 

             ! Neumann BC
             ELSEIF (BC_Switch(-elemID(1,i)) .eq. NeumannBC) THEN
               dcomy(1:Knod,0,i) = BC_Values(1:Knod,-elemID(1,i)) * Sqrt(Acoef(1:Knod))

               NeuRHS(1:Knod) = Jac(1:Knod)*dcomy(1:Knod,0,i) - Acoef(1:Knod)*(dyuBy(1:Knod,0,i) + uBy(1:Knod,0,i)*gLprime)
               DO j = 1, Knod
                 NeuMatrix(j,1:Knod) = -Bcoef(j) * SolnNodesGradLgrangeBasis(1:Knod,j)
                 NeuMatrix(j,j) = NeuMatrix(j,j) - Acoef(j)*gLprime
               ENDDO

               CALL Gauss_Solver(Knod, NeuMatrix, comy(1:Knod,0,i), NeuRHS)

               VelocJump(1:Knod,-elemID(1,i)) = comy(1:Knod,0,i)

             ENDIF
           ENDIF
         ENDDO
       ENDIF

       DO i = 1, Nel

         DO jy = 1, Lnod
           DO jx = 1, Lnod
             xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
             yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
           ENDDO
         ENDDO

         ! Get first derivatives of the variable
         DO jy = 1, Knod
           DO jx = 1, Knod
             ! Get geometry stuff at all solution nodes (per element)
             sdx_dxsi = 0.d0
             sdx_deta = 0.d0
             sdy_dxsi = 0.d0
             sdy_deta = 0.d0
             DO ly = 1, Lnod
               DO lx = 1, Lnod
                 sdx_dxsi = sdx_dxsi + &
                            ! grad at x-dir                    no grad at y-dir               x/y-coord of geom
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_dxsi = sdy_dxsi + &
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
                 sdx_deta = sdx_deta + &
                            ! no grad at x-dir             grad at y-dir                      x/y-coord of geom
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_deta = sdy_deta + &
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
               ENDDO
             ENDDO
             sJac(jx,jy) = sdx_dxsi * sdy_deta - sdx_deta * sdy_dxsi

             ! Get grads of unknown uin along xsi (for eta=const) and along eta (for xsi=const
             sdu_dxsi = dot_product(uin(1:Knod,jy,i), SolnNodesGradLgrangeBasis(1:Knod,jx)) + &
                        dot_product(comx(jy,0:1,i) - uBx(jy,0:1,i), NodesGradRadau(jx,0:1))
             sdu_deta = dot_product(uin(jx,1:Knod,i), SolnNodesGradLgrangeBasis(1:Knod,jy)) + &
                        dot_product(comy(jx,0:1,i) - uBy(jx,0:1,i), NodesGradRadau(jy,0:1))

             ! Get f~ and g~ as per Huynh's paper
             f_tilda(jx,jy) = (sdu_dxsi * (sdx_deta*sdx_deta + sdy_deta*sdy_deta) - &
                               sdu_deta * (sdx_dxsi*sdx_deta + sdy_dxsi*sdy_deta)) / sJac(jx,jy)
             g_tilda(jx,jy) = (sdu_deta * (sdx_dxsi*sdx_dxsi + sdy_dxsi*sdy_dxsi) - &
                               sdu_dxsi * (sdx_dxsi*sdx_deta + sdy_dxsi*sdy_deta)) / sJac(jx,jy)
           ENDDO
         ENDDO

         ! Get ?_tildas at the mesh boundaries
         DO j = 1, Knod
           f_tildaB(j,0) = dot_product(f_tilda(1:Knod,j), SolnBndryLgrangeBasis(1:Knod,0))      ! Left of mesh - derivative
           f_tildaB(j,1) = dot_product(f_tilda(1:Knod,j), SolnBndryLgrangeBasis(1:Knod,1))      ! Right of mesh - derivative
           g_tildaB(j,0) = dot_product(g_tilda(j,1:Knod), SolnBndryLgrangeBasis(1:Knod,0))      ! South of mesh - derivative
           g_tildaB(j,1) = dot_product(g_tilda(j,1:Knod), SolnBndryLgrangeBasis(1:Knod,1))      ! North of mesh - derivative
         ENDDO

         ! Now get the Laplacian
         DO jy = 1, Knod
           DO jx = 1, Knod
             ddu(jx,jy,i) = (dot_product(f_tilda(1:Knod,jy), SolnNodesGradLgrangeBasis(1:Knod,jx))   + &
                             dot_product(dcomx(jy,0:1,i) - f_tildaB(jy,0:1), NodesGradRadau(jx,0:1)) + &
                             dot_product(g_tilda(jx,1:Knod), SolnNodesGradLgrangeBasis(1:Knod,jy))   + &
                             dot_product(dcomy(jx,0:1,i) - g_tildaB(jx,0:1), NodesGradRadau(jy,0:1))) / sJac(jx,jy)
           ENDDO
         ENDDO

       ENDDO
  
END SUBROUTINE GetDiffusedFlux


SUBROUTINE GetConvectedFlux(HuynhSolver_type, uin, du)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type
       real*8, dimension(Knod, Knod, Nel) :: uin, du

       integer i, j, ip, jx, jy, lx, ly
       real*8, dimension(Knod, Knod, Nel) :: DFluxX, DFluxY            ! nodal values of Discontinous Fluxes in X and Y dirs
       real*8, dimension(Knod, 0:1, Nel) :: uBx,  uBy                  ! u at the L/R (Bx) and S/N (By) boundaries of a mesh
                                                                       ! com(mon) values at the interface along x and y
       real*8, dimension(Knod, 0:1, Nel):: FlxUBx, FlxUBy, FlxVBx, FlxVBy      ! Flux at the interfaces along x and y directions
       real*8, dimension(Knod, 0:1, Nel):: jFlxBx, jFlxBy              ! Flux jumps at the interfaces along x and y directions

       real*8, dimension(Knod) :: Jac, dx_dxsi, dx_deta, dy_dxsi, dy_deta  ! coord transformation stuff
       real*8, dimension(Knod, Knod) :: sJac, f_tilda, g_tilda
       real*8, dimension(Knod, 0:1) :: f_tildaB, g_tildaB

       real*8, dimension(Lnod, Lnod) :: xloc, yloc

       real*8 sdx_dxsi, sdx_deta, sdy_dxsi, sdy_deta, sdu_dxsi, sdu_deta
       real*8 au, av, upwind, conv
       real*8 tmp(Knod)


       uBx = 0.d0
       uBy = 0.d0
       DFluxX = 0.d0
       DFluxY = 0.d0
       FlxUBx = 0.d0
       FlxUBy = 0.d0
       FlxVBx = 0.d0
       FlxVBy = 0.d0
       jFlxBx = 0.d0
       jFlxBy = 0.d0

       DO i = 1, Nel

         DO jy = 1, Lnod
           DO jx = 1, Lnod
             xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
             yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
           ENDDO
         ENDDO

         ! Get first derivatives of the variable
         DO jy = 1, Knod
           DO jx = 1, Knod
             ! Get geometry stuff at all solution nodes (per element)
             sdx_dxsi = 0.d0
             sdx_deta = 0.d0
             sdy_dxsi = 0.d0
             sdy_deta = 0.d0
             DO ly = 1, Lnod
               DO lx = 1, Lnod
                 sdx_dxsi = sdx_dxsi + &
                            ! grad at x-dir                    no grad at y-dir               x/y-coord of geom
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_dxsi = sdy_dxsi + &
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
                 sdx_deta = sdx_deta + &
                            ! no grad at x-dir             grad at y-dir                      x/y-coord of geom
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_deta = sdy_deta + &
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
               ENDDO
             ENDDO

             ! Get f~ and g~ as per Huynh's paper
             DFluxX(jx,jy,i) =  sdy_deta * Uvel(jx,jy,i) - sdx_deta * Vvel(jx,jy,i)
             DFluxY(jx,jy,i) = -sdy_dxsi * Uvel(jx,jy,i) + sdx_dxsi * Vvel(jx,jy,i)
           ENDDO
         ENDDO

       ENDDO

       ! Extrapolate the unknown, uin, and its derivative to the mesh boundaries using Lagrange polynomials of order Knod-1
       DO i = 1, Nel
         DO j = 1, Knod
           uBx(j,0,i) = dot_product(uin(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,0))            ! Left of mesh  - value
           uBx(j,1,i) = dot_product(uin(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,1))            ! Right of mesh - value
           uBy(j,0,i) = dot_product(uin(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,0))            ! South of mesh  - value
           uBy(j,1,i) = dot_product(uin(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,1))            ! North of mesh  - value

           FlxUBx(j,0,i) = dot_product(DFluxX(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,0))    ! Left of mesh  - Flux
           FlxUBx(j,1,i) = dot_product(DFluxX(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,1))    ! Left of mesh  - Flux
           FlxVBx(j,0,i) = dot_product(DFluxY(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,0))    ! Left of mesh  - Flux
           FlxVBx(j,1,i) = dot_product(DFluxY(1:Knod,j,i), SolnBndryLgrangeBasis(1:Knod,1))    ! Left of mesh  - Flux
           FlxUBy(j,0,i) = dot_product(DFluxX(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,0))    ! Left of mesh  - Flux
           FlxUBy(j,1,i) = dot_product(DFluxX(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,1))    ! Left of mesh  - Flux
           FlxVBy(j,0,i) = dot_product(DFluxY(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,0))    ! Left of mesh  - Flux
           FlxVBy(j,1,i) = dot_product(DFluxY(j,1:Knod,i), SolnBndryLgrangeBasis(1:Knod,1))    ! Left of mesh  - Flux

           DFluxX(j,1:Knod,i) = DFluxX(j,1:Knod,i) * uin(j,1:Knod,i)                                       ! Left-to-Right nodal fluxes
           DFluxY(j,1:Knod,i) = DFluxY(j,1:Knod,i) * uin(j,1:Knod,i)                                       ! South-to-North nodal fluxes
         ENDDO
       ENDDO


       ! Get the common values of uin at the mesh interfaces (we use simple averaging in this case)
       ! We definitely need a more efficient strategy - right now we're saving com-s twice, and can probably save on d?uB? as well
       ! We probably can and should write a more compact method so that we don't repeat the same stuff for NEWS (which becomes messier for NEWSBT in 3D)
       DO i = 1, Nel

         DO jy = 1, Lnod
           DO jx = 1, Lnod
             xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
             yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
           ENDDO
         ENDDO

         ! Right face
         dx_dxsi = 0.d0
         dx_deta = 0.d0
         dy_dxsi = 0.d0
         dy_deta = 0.d0
         DO jy = 1, Lnod
           DO jx = 1, Lnod
             dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                               ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                               GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
             dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                               GeomBndryGradLgrangeBasis(jx,1) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             dx_deta(1:Knod) = dx_deta(1:Knod) + &
                               ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                               GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
             dy_deta(1:Knod) = dy_deta(1:Knod) + &
                               GeomBndryLgrangeBasis(jx,1) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
           ENDDO
         ENDDO

         IF (elemID(2,i) .lt. 0) THEN
!           IF (BC_Switch(-elemID(2,i)) .eq. DirichletBC) THEN
!!!! **** Later on this has to be changed into local coordinates similar to Uvel and Vvel - right now, this is zero so it doesn't
!matter
             tmp(1:Knod) = BC_Values(1:Knod,-elemID(2,i))
!           ELSEIF (BC_Switch(-elemID(2,i)) .eq. NeumannBC) THEN
!           ENDIF
         ELSE
           tmp(1:Knod) = FlxUBx(1:Knod,1,i)
         ENDIF
         FlxUBx(1:Knod,1,i) = tmp(1:Knod) * uBx(1:Knod,1,i)
         FlxVBx(1:Knod,1,i) = tmp(1:Knod)


         ! Left face
         dx_dxsi = 0.d0
         dx_deta = 0.d0
         dy_dxsi = 0.d0
         dy_deta = 0.d0
         DO jy = 1, Lnod
           DO jx = 1, Lnod
             dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                               ! grad at x-dir                   no grad along y-dir                x/y-coord of geom
                               GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
             dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                               GeomBndryGradLgrangeBasis(jx,0) * GeomNodesLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
             dx_deta(1:Knod) = dx_deta(1:Knod) + &
                               ! no grad at x-dir            grad along y-dir                       x/y-coord of geom
                               GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * xloc(jx,jy)
             dy_deta(1:Knod) = dy_deta(1:Knod) + &
                               GeomBndryLgrangeBasis(jx,0) * GeomNodesGradLgrangeBasis(jy,1:Knod) * yloc(jx,jy)
           ENDDO
         ENDDO

         IF (elemID(4,i) .lt. 0) THEN
!           IF (BC_Switch(-elemID(4,i)) .eq. DirichletBC) THEN
             tmp(1:Knod) = BC_Values(1:Knod,-elemID(4,i))
!           ELSEIF (BC_Switch(-elemID(4,i)) .eq. NeumannBC) THEN
!           ENDIF
         ELSE
           tmp(1:Knod) = FlxUBx(1:Knod,0,i)
         ENDIF
         FlxUBx(1:Knod,0,i) = tmp(1:Knod) * uBx(1:Knod,0,i)
         FlxVBx(1:Knod,0,i) = tmp(1:Knod)


         ! North face
         dx_dxsi = 0.d0
         dx_deta = 0.d0
         dy_dxsi = 0.d0
         dy_deta = 0.d0
         DO jy = 1, Lnod
           DO jx = 1, Lnod
             dx_deta(1:Knod) = dx_deta(1:Knod) + &
                               ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                               GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
             dy_deta(1:Knod) = dy_deta(1:Knod) + &
                               GeomBndryGradLgrangeBasis(jy,1) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                               ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                               GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
             dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                               GeomBndryLgrangeBasis(jy,1) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
           ENDDO
         ENDDO

         IF (elemID(3,i) .lt. 0) THEN
!           IF (BC_Switch(-elemID(3,i)) .eq. DirichletBC) THEN
             tmp(1:Knod) = BC_Values(1:Knod,-elemID(3,i))
!           ELSEIF (BC_Switch(-elemID(3,i)) .eq. NeumannBC) THEN
!           ENDIF
         ELSE
           tmp(1:Knod) = FlxVBy(1:Knod,1,i) 
         ENDIF
         FlxVBy(1:Knod,1,i) = tmp(1:Knod) * uBy(1:Knod,1,i)
         FlxUBy(1:Knod,1,i) = tmp(1:Knod)

         ! South face
         dx_dxsi = 0.d0
         dx_deta = 0.d0
         dy_dxsi = 0.d0
         dy_deta = 0.d0
         DO jy = 1, Lnod
           DO jx = 1, Lnod
             dx_deta(1:Knod) = dx_deta(1:Knod) + &
                               ! grad at y-dir                   no grad along x-dir                x/y-coord of geom
                               GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
             dy_deta(1:Knod) = dy_deta(1:Knod) + &
                               GeomBndryGradLgrangeBasis(jy,0) * GeomNodesLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
             dx_dxsi(1:Knod) = dx_dxsi(1:Knod) + &
                               ! no grad at y-dir            grad along x-dir                       x/y-coord of geom
                               GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * xloc(jx,jy)
             dy_dxsi(1:Knod) = dy_dxsi(1:Knod) + &
                               GeomBndryLgrangeBasis(jy,0) * GeomNodesGradLgrangeBasis(jx,1:Knod) * yloc(jx,jy)
           ENDDO
         ENDDO

         IF (elemID(1,i) .lt. 0) THEN
!           IF (BC_Switch(-elemID(1,i)) .eq. DirichletBC) THEN
             tmp(1:Knod) = BC_Values(1:Knod,-elemID(1,i))
!           ELSEIF (BC_Switch(-elemID(1,i)) .eq. NeumannBC) THEN
!           ENDIF
         ELSE
           tmp(1:Knod) = FlxVBy(1:Knod,0,i) 
         ENDIF
         FlxVBy(1:Knod,0,i) = tmp(1:Knod) * uBy(1:Knod,0,i)
         FlxUBy(1:Knod,0,i) = tmp(1:Knod)

       ENDDO

       DO i = 1, Nel
         ! Right face
         IF (elemID(2,i) .gt. 0) THEN
           DO j = 1, Knod
             IF (Abs(uBx(j,1,i) - uBx(j,0,elemID(2,i))) .gt. 1.0d-10) THEN
               conv = (FlxUBx(j,1,i) - FlxUBx(j,0,elemID(2,i))) / (uBx(j,1,i) - uBx(j,0,elemID(2,i)))
             ELSE
               conv = FlxVBx(j,1,i)
             ENDIF
             upwind = 0.5d0 * ( FlxUBx(j,1,i) + FlxUBx(j,0,elemID(2,i)) + Abs(conv) * ( uBx(j,1,i) - uBx(j,0,elemID(2,i)) ) )
             jFlxBx(j,1,i) = upwind !- jFlxBx(j,1,i)
           ENDDO
         ELSEIF (elemID(2,i) .lt. 0) THEN
           ! Dirichlet BC
!           IF (BC_Switch(-elemID(2,i)) .eq. DirichletBC) THEN
             jFlxBx(1:Knod,1,i) = FlxUBx(1:Knod,1,i) !- jFlxBx(j,1,i)
           ! Neumann BC
!           ELSEIF (BC_Switch(-elemID(2,i)) .eq. NeumannBC) THEN
!           ENDIF
         ENDIF

         ! Left face
         IF (elemID(4,i) .gt. 0) THEN
           DO j = 1, Knod
             IF (Abs(uBx(j,0,i) - uBx(j,1,elemID(4,i))) .gt. 1.0d-10) THEN
               conv = (FlxUBx(j,0,i) - FlxUBx(j,1,elemID(4,i))) / (uBx(j,0,i) - uBx(j,1,elemID(4,i)))
             ELSE
               conv = FlxVBx(j,0,i)
             ENDIF
             upwind = 0.5d0 * ( FlxUBx(j,0,i) + FlxUBx(j,1,elemID(4,i)) - Abs(conv) * ( uBx(j,0,i) - uBx(j,1,elemID(4,i)) ) )
             jFlxBx(j,0,i) = upwind !- jFlxBx(j,0,i)
           ENDDO
         ELSEIF (elemID(4,i) .lt. 0) THEN
           ! Dirichlet BC
!           IF (BC_Switch(-elemID(4,i)) .eq. DirichletBC) THEN
             jFlxBx(1:Knod,0,i) = FlxUBx(1:Knod,0,i) !- jFlxBx(j,0,i)
           ! Neumann BC
!           ELSEIF (BC_Switch(-elemID(4,i)) .eq. NeumannBC) THEN
!           ENDIF
         ENDIF

         ! North face
         IF (elemID(3,i) .gt. 0) THEN
           DO j = 1, Knod
             IF (Abs(uBy(j,1,i) - uBy(j,0,elemID(3,i))) .gt. 1.0d-10) THEN
               conv = (FlxVBy(j,1,i) - FlxVBy(j,0,elemID(3,i))) / (uBy(j,1,i) - uBy(j,0,elemID(3,i)))
             ELSE
               conv = FlxUBy(j,1,i)
             ENDIF
             upwind = 0.5d0 * ( FlxVBy(j,1,i) + FlxVBy(j,0,elemID(3,i)) + Abs(conv) * ( uBy(j,1,i) - uBy(j,0,elemID(3,i)) ) )
             jFlxBy(j,1,i) = upwind !- jFlxBy(j,1,i)
           ENDDO
         ELSEIF (elemID(3,i) .lt. 0) THEN
           ! Dirichlet BC
!           IF (BC_Switch(-elemID(3,i)) .eq. DirichletBC) THEN
             jFlxBy(1:Knod,1,i) = FlxVBy(1:Knod,1,i) !- jFlxBy(j,1,i)
           ! Neumann BC
!           ELSEIF (BC_Switch(-elemID(3,i)) .eq. NeumannBC) THEN
!           ENDIF
         ENDIF

         ! South face
         IF (elemID(1,i) .gt. 0) THEN
           DO j = 1, Knod
             IF (Abs(uBy(j,0,i) - uBy(j,1,elemID(1,i))) .gt. 1.0d-10) THEN
               conv = (FlxVBy(j,0,i) - FlxVBy(j,1,elemID(1,i))) / (uBy(j,0,i) - uBy(j,1,elemID(1,i)))
             ELSE
               conv = FlxUBy(j,0,i)
             ENDIF
             upwind = 0.5d0 * ( FlxVBy(j,0,i) + FlxVBy(j,1,elemID(1,i)) - Abs(conv) * ( uBy(j,0,i) - uBy(j,1,elemID(1,i)) ) )
             jFlxBy(j,0,i) = upwind !- jFlxBy(j,0,i)
           ENDDO
         ELSEIF (elemID(1,i) .lt. 0) THEN
           ! Dirichlet BC
!           IF (BC_Switch(-elemID(1,i)) .eq. DirichletBC) THEN
             jFlxBy(1:Knod,0,i) = FlxVBy(1:Knod,0,i) !- jFlxBy(j,0,i)
           ! Neumann BC
!           ELSEIF (BC_Switch(-elemID(1,i)) .eq. NeumannBC) THEN
!           ENDIF
         ENDIF

       ENDDO


       DO i = 1, Nel

         DO jy = 1, Lnod
           DO jx = 1, Lnod
             xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
             yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
           ENDDO
         ENDDO

         ! Get first derivatives of the variable
         DO jy = 1, Knod
           DO jx = 1, Knod
             ! Get geometry stuff at all solution nodes (per element)
             sdx_dxsi = 0.d0
             sdx_deta = 0.d0
             sdy_dxsi = 0.d0
             sdy_deta = 0.d0
             DO ly = 1, Lnod
               DO lx = 1, Lnod
                 sdx_dxsi = sdx_dxsi + &
                            ! grad at x-dir                    no grad at y-dir               x/y-coord of geom
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_dxsi = sdy_dxsi + &
                            GeomNodesGradLgrangeBasis(lx,jx) * GeomNodesLgrangeBasis(ly,jy) * yloc(lx,ly)
                 sdx_deta = sdx_deta + &
                            ! no grad at x-dir             grad at y-dir                      x/y-coord of geom
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * xloc(lx,ly)
                 sdy_deta = sdy_deta + &
                            GeomNodesLgrangeBasis(lx,jx) * GeomNodesGradLgrangeBasis(ly,jy) * yloc(lx,ly)
               ENDDO
             ENDDO
             sJac(jx,jy) = sdx_dxsi * sdy_deta - sdx_deta * sdy_dxsi

             ! Get f~ and g~ as per Huynh's paper
             f_tilda(jx,jy) = DFluxX(jx,jy,i)
             g_tilda(jx,jy) = DFluxY(jx,jy,i)
           ENDDO
         ENDDO

         ! Get ?_tildas at the mesh boundaries
         DO j = 1, Knod
           f_tildaB(j,0) = dot_product(f_tilda(1:Knod,j), SolnBndryLgrangeBasis(1:Knod,0))      ! Left of mesh - derivative
           f_tildaB(j,1) = dot_product(f_tilda(1:Knod,j), SolnBndryLgrangeBasis(1:Knod,1))      ! Right of mesh - derivative
           g_tildaB(j,0) = dot_product(g_tilda(j,1:Knod), SolnBndryLgrangeBasis(1:Knod,0))      ! South of mesh - derivative
           g_tildaB(j,1) = dot_product(g_tilda(j,1:Knod), SolnBndryLgrangeBasis(1:Knod,1))      ! North of mesh - derivative
         ENDDO

         ! Now get the convection
         DO jy = 1, Knod
           DO jx = 1, Knod
!             du(jx,jy,i) = du(jx,jy,i) -  &
             du(jx,jy,i) = - &
                             (dot_product(f_tilda(1:Knod,jy), SolnNodesGradLgrangeBasis(1:Knod,jx))   + &
                              dot_product(jFlxBx(jy,0:1,i) - f_tildaB(jy,0:1), NodesGradRadau(jx,0:1)) + &
                              dot_product(g_tilda(jx,1:Knod), SolnNodesGradLgrangeBasis(1:Knod,jy))   + &
                              dot_product(jFlxBy(jx,0:1,i) - g_tildaB(jx,0:1), NodesGradRadau(jy,0:1))) / sJac(jx,jy)
           ENDDO
         ENDDO

       ENDDO
  
END SUBROUTINE GetConvectedFlux


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


SUBROUTINE dumpResult(numStep, Reyn, dt, HuynhSolver_type, tIntegrator_type, prob_type)
       USE params
       USE variables
       USE CnvrtTensor2FemIndices

       implicit NONE
       integer HuynhSolver_type, tIntegrator_type, numStep, prob_type
       real*8 Reyn, dt

       integer i, jx, jy, kx, ky, ntime
       real*8 x, y, time, timep, Vort_ex, Uvel_ex, Linf_Norm, L1_Norm, L2_Norm
       character(132) filename

       real*8 theta2,theta1,thetam,thetad,rad2,rad1,radm,radd
       real*8, dimension(Lnod,Lnod) :: xloc, yloc
       real*8, dimension(Lnod) :: xloc1D, yloc1D

   INTERFACE
     FUNCTION Vort_xct(x, y, t, Re, prob_type)
       implicit NONE
       real*8 x, y, t, Re, Vort_xct
       integer prob_type
     END FUNCTION Vort_xct
     FUNCTION Uvel_xct(x, y, t, Re, prob_type)
       implicit NONE
       real*8 x, y, t, Re, Uvel_xct
       integer prob_type
     END FUNCTION Uvel_xct
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
       DO i = 1, Nel
         DO ky = 1, Knod
           DO kx = 1, Knod

             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
                 yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
               ENDDO
             ENDDO

             ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
             x = 0.d0
             y = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                         ! along x-dir                  along y-dir                    x-coord of geom
                 x = x + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * xloc(jx,jy)
                 y = y + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * yloc(jx,jy)
               ENDDO
             ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             y = (radm + sps(ky)*radd) * sin(thetam + sps(kx)*thetad)
             x = (radm + sps(ky)*radd) * cos(thetam + sps(kx)*thetad)
 endif

             Vort_ex = Vort_xct(x, y, timep, Reyn, prob_type)
             Linf_Norm = Max(Linf_Norm, Abs(Vort_ex-Vort(kx,ky,i)))
             L1_Norm = L1_Norm + Abs(Vort_ex-Vort(kx,ky,i))
             L2_Norm = L2_Norm + (Vort_ex-Vort(kx,ky,i))**2
!             write(*,*) x, y, Vort(i,kx,ky), Vort_ex, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
             write(9,*) x, y, Vort(kx,ky,i), Vort_ex, abs(Vort_ex-Vort(kx,ky,i)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
!             write(9,*) x, y, Vort(kx,ky,i), Uvel(kx,ky,i), Vvel(kx,ky,i)
           ENDDO
           write(9,*) ' '
         ENDDO
       ENDDO

       close(unit=9)

       open(unit=9, file=trim(filename)//'.stat',status='unknown')

       IF (.TRUE.) THEN

       L1_Norm = L1_Norm / (Nel*Knod**2)
       L2_Norm = Sqrt(L2_Norm / (Nel*Knod**2))
       write(9,*) trim('dt, 1.d0/Nel, Nel, Nel*Knod^2, Linf_Norm, L1_Norm, L2_Norm')
       write(9,*) dt, 1.d0/Nel, Nel, Nel*Knod**2, Linf_Norm, L1_Norm, L2_Norm
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod**2, Linf_Norm, L1_Norm, L2_Norm
       close(unit=9)

       open(unit=9, file=trim(filename)//'_du.csv',status='unknown')

       Linf_Norm = -1.d0
       L1_Norm = 0.d0
       L2_Norm = 0.d0
       DO i = 1, Nel
         DO ky = 1, Knod
           DO kx = 1, Knod

             DO jy = 1, Lnod
               DO jx = 1, Lnod
                 xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
                 yloc(jx,jy) = ycoord(nodeID(t2f(jx,jy),i))
               ENDDO
             ENDDO

             ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
             x = 0.d0
             y = 0.d0
             DO jy = 1, Lnod
               DO jx = 1, Lnod
                         ! along x-dir                  along y-dir                    x-coord of geom
                 x = x + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * xloc(jx,jy)
                 y = y + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * yloc(jx,jy)
               ENDDO
             ENDDO

 if (exact) then
             if (abs(yloc(Lnod,1)) < 1.d-12) yloc(Lnod,1) = 0.d0
             if (abs(xloc(Lnod,1)) < 1.d-12) xloc(Lnod,1) = 0.d0
             if (abs(yloc(1,1)) < 1.d-12) yloc(1,1) = 0.d0
             if (abs(xloc(1,1)) < 1.d-12) xloc(1,1) = 0.d0
             theta2 = atan2(yloc(Lnod,1),xloc(Lnod,1))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc(1,1),xloc(1,1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

             y = (radm + sps(ky)*radd) * sin(thetam + sps(kx)*thetad)
             x = (radm + sps(ky)*radd) * cos(thetam + sps(kx)*thetad)
 endif

             Uvel_ex = Uvel_xct(x, y, timep, Reyn, prob_type)
             Linf_Norm = Max(Linf_Norm, Abs(Uvel_ex-Uvel(kx,ky,i)))
             L1_Norm = L1_Norm + Abs(Uvel_ex-Uvel(kx,ky,i))
             L2_Norm = L2_Norm + (Uvel_ex-Uvel(kx,ky,i))**2
!             Linf_Norm = Max(Linf_Norm, Abs(Uvel_ex-Vvel(kx,ky,i)))
!             L1_Norm = L1_Norm + Abs(Uvel_ex-Vvel(kx,ky,i))
!             L2_Norm = L2_Norm + (Uvel_ex-Vvel(kx,ky,i))**2
!             write(*,*) x, y, Vort(i,kx,ky), Vort_ex, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
             write(9,*) x, y, Vvel(kx,ky,i), Uvel_ex, abs(Uvel_ex-Uvel(kx,ky,i)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
!             write(9,*) x, y, Vvel(kx,ky,i), Uvel_ex, abs(Uvel_ex-Vvel(kx,ky,i)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
!             write(9,*) x, y, Vort(kx,ky,i), Uvel(kx,ky,i), Vvel(kx,ky,i)
           ENDDO
           write(9,*) ' '
         ENDDO
       ENDDO

       close(unit=9)

       open(unit=9, file=trim(filename)//'_du.stat',status='unknown')

       ENDIF

       L1_Norm = L1_Norm / (Nel*Knod**2)
       L2_Norm = Sqrt(L2_Norm / (Nel*Knod**2))
       write(9,*) trim('dt, 1.d0/Nel, Nel, Nel*Knod^2, Linf_Norm, L1_Norm, L2_Norm')
       write(9,*) dt, 1.d0/Nel, Nel, Nel*Knod**2, Linf_Norm, L1_Norm, L2_Norm
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod**2, Linf_Norm, L1_Norm, L2_Norm
       close(unit=9)

       open(unit=9, file=trim(filename)//'_boundary.csv',status='unknown')

       Linf_Norm = -1.d0
       L1_Norm = 0.d0
       L2_Norm = 0.d0
       DO i = 1, Nelb/2

         DO jx = 1, Lnod
           xloc1D(jx) = xcoord(boundarynodeID(t2f(jx),i))
           yloc1D(jx) = ycoord(boundarynodeID(t2f(jx),i))
         ENDDO

         DO kx = 1, Knod
           ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
           x = 0.d0
           y = 0.d0
           DO jx = 1, Lnod
                     ! along x-dir                  x-coord of geom
             x = x + GeomNodesLgrangeBasis(jx,kx) * xloc1D(jx)
             y = y + GeomNodesLgrangeBasis(jx,kx) * yloc1D(jx)
           ENDDO

 if (exact) then
             if (abs(yloc1D(Lnod)) < 1.d-12) yloc1D(Lnod) = 0.d0
             if (abs(xloc1D(Lnod)) < 1.d-12) xloc1D(Lnod) = 0.d0
             if (abs(yloc1D(1)) < 1.d-12) yloc1D(1) = 0.d0
             if (abs(xloc1D(1)) < 1.d-12) xloc1D(1) = 0.d0
             theta2 = atan2(yloc1D(Lnod),xloc1D(Lnod))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc1D(1),xloc1D(1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

!             y = (radm + sps(ky)*radd) * sin(thetam + sps(kx)*thetad)
!             x = (radm + sps(ky)*radd) * cos(thetam + sps(kx)*thetad)
             y = 0.5d0 * sin(thetam + sps(kx)*thetad)
             x = 0.5d0 * cos(thetam + sps(kx)*thetad)
 endif

             Vort_ex = Vort_xct(x, y, timep, Reyn, prob_type)
             Linf_Norm = Max(Linf_Norm, Abs(Vort_ex-VelocJump(kx,i)))
             L1_Norm = L1_Norm + Abs(Vort_ex-VelocJump(kx,i))
             L2_Norm = L2_Norm + (Vort_ex-VelocJump(kx,i))**2
!             write(*,*) x, y, Vort(i,kx,ky), Vort_ex, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
             write(9,*) x, y, VelocJump(kx,i), Vort_ex, abs(Vort_ex-VelocJump(kx,i)) !, 100*(Vort_ex-Vort(i,kx,ky))/Vort_ex
!             write(9,*) x, y, Vort(kx,ky,i), Uvel(kx,ky,i), Vvel(kx,ky,i)
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
       DO i = 1, Nelb

         DO jx = 1, Lnod
           xloc1D(jx) = xcoord(boundarynodeID(t2f(jx),i))
           yloc1D(jx) = ycoord(boundarynodeID(t2f(jx),i))
         ENDDO

         DO kx = 1, Knod
           ! interpolate x/y-coord of sol pt at (jx,jy) using nodal coordinates of element i
           x = 0.d0
           y = 0.d0
           DO jx = 1, Lnod
                     ! along x-dir                  x-coord of geom
             x = x + GeomNodesLgrangeBasis(jx,kx) * xloc1D(jx)
             y = y + GeomNodesLgrangeBasis(jx,kx) * yloc1D(jx)
           ENDDO

 if (exact) then
             if (abs(yloc1D(Lnod)) < 1.d-12) yloc1D(Lnod) = 0.d0
             if (abs(xloc1D(Lnod)) < 1.d-12) xloc1D(Lnod) = 0.d0
             if (abs(yloc1D(1)) < 1.d-12) yloc1D(1) = 0.d0
             if (abs(xloc1D(1)) < 1.d-12) xloc1D(1) = 0.d0
             theta2 = atan2(yloc1D(Lnod),xloc1D(Lnod))
             if (theta2 < 0.d0) theta2 = theta2 + 2.d0*pi
             theta1 = atan2(yloc1D(1),xloc1D(1))
             if (theta1 <= 0.d0) theta1 = theta1 + 2.d0*pi
             thetam = 0.5d0 * (theta2 + theta1)
             thetad = 0.5d0 * (theta2 - theta1)
             rad2 = sqrt(xloc(1,Lnod)**2 + yloc(1,Lnod)**2)
             rad1 = sqrt(xloc(1,1)**2 + yloc(1,1)**2)
             radm = 0.5d0 * (rad2 + rad1)
             radd = 0.5d0 * (rad2 - rad1)

!             y = (radm + sps(ky)*radd) * sin(thetam + sps(kx)*thetad)
!             x = (radm + sps(ky)*radd) * cos(thetam + sps(kx)*thetad)
             y = 0.5d0 * sin(thetam + sps(kx)*thetad)
             x = 0.5d0 * cos(thetam + sps(kx)*thetad)
 endif

             Uvel_ex = Uvel_xct(x, y, timep, Reyn, prob_type)
             Linf_Norm = Max(Linf_Norm, Abs(Uvel_ex-VelocJump(kx,i)))
             L1_Norm = L1_Norm + Abs(Uvel_ex-VelocJump(kx,i))
             L2_Norm = L2_Norm + (Uvel_ex-VelocJump(kx,i))**2
             write(9,*) x, y, VelocJump(kx,i), Uvel_ex, abs(Uvel_ex-VelocJump(kx,i))
           ENDDO
           write(9,*) ' '
         ENDDO
       L1_Norm = L1_Norm / (Nel*Knod)
       L2_Norm = Sqrt(L2_Norm / (Nel*Knod))
       write(*,*) dt, 1.d0/Nel, Nel, Nel*Knod, Linf_Norm, L1_Norm, L2_Norm

       close(unit=9)

CONTAINS

FUNCTION itoa(i) RESULT(res)

       character(11) res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
       write(tmp,'(i0)') i
!       res = trim(tmp)
       res = tmp

END FUNCTION itoa

END SUBROUTINE dumpResult


FUNCTION Vort_xct(x, y, t, Re, prob_type)
       implicit NONE
       real*8, parameter :: pi = 3.1415926535897931d0
       real*8 x, y, t, Re, Vort_xct
       integer prob_type

       integer i
       real*8 error, tmp, r, xm, ym

       IF (prob_type .eq. 1) THEN
         Vort_xct = Sin(pi*y)  / pi**2

       ELSEIF (prob_type .eq. 3) THEN
         Vort_xct = Sin(pi*y)  / pi**2

       ELSEIF (prob_type .eq. 6) THEN
!         Vort_xct = x * (1.d0 + Cos(pi*x)) * Sin(pi*y)
         Vort_xct = y * (1.d0 + Cos(pi*y)) * Sin(pi*x)

       ELSEIF (prob_type .eq. 4) THEN
!         xm = 1.d0 - x
!         Vort_xct = xm * (1.d0 + Cos(pi*xm)) * Sin(pi*y)
         ym = 1.d0 - y
         Vort_xct = ym * (1.d0 + Cos(pi*ym)) * Sin(pi*x)

       ELSEIF (prob_type .eq. 7) THEN
!         Vort_xct = 1.d0 + 0.5d0 * Log(x*x + y*y) / Log(2.d0)
         Vort_xct = - 0.5d0 * Log(x*x + y*y) / Log(2.d0)

       ELSEIF (prob_type .eq. 8) THEN
         r = sqrt(x*x + y*y)
         Vort_xct = - (r**3 - Log(r**3)/8.d0 -1.d0) / 9.d0

       ENDIF
      
END FUNCTION Vort_xct


FUNCTION Uvel_xct(x, y, t, Re, prob_type)
       implicit NONE
       real*8, parameter :: pi = 3.1415926535897931d0
       real*8 x, y, t, Re, Uvel_xct
       integer prob_type

       integer i
       real*8 error, tmp, r, xm, ym

       IF (prob_type .eq. 1) THEN
         Uvel_xct = Cos(pi*y) / pi
       ELSEIF (prob_type .eq. 3) THEN
         Uvel_xct = Cos(pi*y) / pi

       ELSEIF (prob_type .eq. 6) THEN
!         Uvel_xct = (1.d0 + Cos(pi*x) - pi*x*Sin(pi*x)) * Sin(pi*y)
!         Uvel_xct = (1.d0 + Cos(pi*y) - pi*y*Sin(pi*y)) * Sin(pi*x)
         Uvel_xct = (1.d0 + Cos(pi*y) - pi*y*Sin(pi*y)) * Sin(pi*x)
!         Uvel_xct = - pi * (pi*y*(1.d0 + 2.d0*Cos(pi*y)) + 2.d0*Sin(pi*y)) * Sin(pi*x)
       ELSEIF (prob_type .eq. 4) THEN
!         xm = 1.d0 - x
!         Uvel_xct = -(1.d0 + Cos(pi*xm) - pi*xm*Sin(pi*xm)) * Sin(pi*y)
         ym = 1.d0 - y
         Uvel_xct = -(1.d0 + Cos(pi*ym) - pi*ym*Sin(pi*ym)) * Sin(pi*x)
!         Uvel_xct = ( (1.d0 - 5.d0*x**4)/20.d0 + (-pi*Sin(pi*x) + 2.d0)/pi**2 ) * &
!                    ( (y - y**5)/20.d0 + (Cos(pi*y) + 2.d0*y - 1.d0)/pi**2 )
       ELSEIF (prob_type .eq. 7) THEN
         r = x*x + y*y
!         Uvel_xct = y / (r * Log(2.d0))
         Uvel_xct = - y / (r * Log(2.d0))
!!         Uvel_xct = 0.d0
!!         Uvel_xct = x / (r * Log(2.d0))
!         Uvel_xct = 1.d0 / (Sqrt(r) * Log(2.d0))
       ELSEIF (prob_type .eq. 8) THEN
         r = sqrt(x*x + y*y)
         Uvel_xct = y / (3.d0 * r) * (1.d0 / (8.d0 * r) - r*r)
       ENDIF
      
END FUNCTION Uvel_xct


