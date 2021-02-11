#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include "preprocess.h"
#include <Eigen/Eigenvalues> 
#include<Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/OrderingMethods>
#include<Eigen/IterativeLinearSolvers>
typedef Eigen::Triplet<double> Trplet;

#define DirichletBC 0
#define NeumannBC 1

struct LR_boundary { // for the shape functions and Radua where everything reverts back to 1D and on the boundaries
	double left=1., right=1.;
};

// the main class to hold the primitve and fundamental variables and methods
class HO_2D {

private:
	Mesh mesh;
	unsigned int Knod; //Number of sps in each direction,i.e. number of gauss point in each cell in each direction
	double g_prime[2]; //derivative of the correction function g on the [0]: left and [1]: right boundaries
	unsigned char* BC_switch; //specify if a boundary edge is dirichlet(0 or neumann(1) BC
	double** BC_values; // BC values for the Knod points of NelB boundary elements (in terms of velocity or gradient)
	double** BC_parl_vel, ** BC_normal_vel; //the parallel and normal contravariant flux velocity on all global boundary elements
	double** BC_psi; // wall normal velocity (in terms of psi)
	bool* BC_no_slip; //if the boundary conditions are no slip wall. if yes, then it is true
	double*** vorticity; //the vorticity field: first index goes for element, second and third go for j (eta dir) and i (csi dir) sps
	double*** initial_vorticity; //The initial (t=0) vorticity field
	double*** stream_function; //The stream function field
	Cmpnts2*** velocity_cart; //The velocity field
	double Reyn_inv; // the inverse of the Reynolds number
	int HuynhSolver_type; //based on Huynh's scheme type in Table 6.1 of his diffusion paper; types: 2 = stndard DG; 11 = higher order
	int time_integration_type; //time integration method; 1 = Euler; 2 = Mod.Euler; 4 = RK4
	unsigned int problem_type; //Solve different problems (types 1 and 2). prob_type=10 reads mesh from file
	unsigned int num_time_steps; //total number of time steps to march in time
	double dt; //time step size
	unsigned int dump_frequency; //the frequency of time saving
	bool fast; //0 for original, 1 for compact formulation (use gLB, gRB)
	unsigned int N_gps; //number of geometry nodes in the domain
	unsigned int N_el_boundary; //number of edges on the boundaries
	double *sps_local_coor, *sps_weight, *gps_local_coor; //the arrays to store the Gauss-Legendre and their weight, geometry points
	double** sps_boundary_basis; //The value Lagrange shape function of Knod solution points on the left[0](csi=-1) and right[1] (csi=+1) boundaries
	double** sps_boundary_grad_basis; //The derivative of Lagrange shape function of Knod solution points on the left[0](csi=-1) and right[1] (csi=+1) boundaries
	double** sps_sps_grad_basis;   //derivative of Knod solution points on the Knod solution points. The Knod sps is the first index(say i) and other nodes are second index(say j) : sps_sps_grad_basis[i][j]
	double** gps_boundary_basis; //The value of Lagrange shape function of Lnod geometry nodes on the left[0](csi=-1) and right[1] (csi=+1) boundaries
	double** gps_sps_basis;		//value of shape functions of Lnod gps on the Knod sps. size is [Lnod][Knod]
	double** gps_boundary_grad_basis; //derivative of shape functions of Lnod gps on the left and right boundaries, so the size is[Lnod][2]
	double** gps_sps_grad_basis;		//derivative of shape functions of Lnod gps on the Knod sps. size is [Lnod][Knod]
	double** sps_radau; //The left and right Radau function value on the Knod sps points. Left radau=.left, ...
	double** sps_grad_radau; //The derivative of the left and right Radau function on the Knod sps points. Left radau=.left, ...
	Cmpnts2*** vol_Dx_Dxsi; //the derivative of x wrt to xsi_s(s is 0, 1 to make dx / dxsi, dx / deta) at the sps(jx, jy) in tensor form on element el.This forms Vol_Dx_iDxsi_j(el, jy, jx).x, .y
	Cmpnts2*** vol_Dy_Dxsi; //the derivative of y wrt to xsi_s(s is 0, 1 to make dy / dxsi, dy / deta) at the sps(jx, jy) in tensor form on element el.This forms Vol_Dy_iDxsi_j(el, jy, jx).x, .y
	double*** vol_jac;
	double*** face_Acoef, *** face_Bcoef, *** face_Anorm, *** face_jac; // all have the size(N, 4, Knod)
	double***** G; // G[el][j][i]][d1=0,1][d2=0,1] is (g_d1 dot g_d2)/J at j,i of element el
	double***** GB; // GB[el][d=0,1][t=0,1]][j][d2=0,1] is (g_d dot g_d2)/J at j flux point on the face of element el that is in the d direction and t side.
	double*** RHS_advective; //the right hand side of the vorticity eq for the advective term: -div(vw) = -d_dxsi_j(V^jw)/G
	double*** RHS_diffusive;   //stores 1/Re * Laplace(w) in RHS_diffusive[el][ky][kx]
	double** velocity_jump;
	double*** boundary_source; //Poisson Solver's RHS term to be dotted by BC_Values of the edge that belongs to the element (which has this edge on the global boundary); size is [N_edges_boundary][Knod*Knod][Knod]
	int LHS_type; //1 is eigen, 2 is hypre
	int nnz; //number of non-zeros in the LHS matrix of poisson equation
	Eigen::SparseMatrix<double> LHS_Eigen; //to store the poisson LHS in Eigen format
	Eigen::VectorXd RHS_Eigen; //the right hand side in the poisson equation discretization
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double, Eigen::Lower | Eigen::Upper, Eigen::AMDOrdering<int> > > cg_Eigen;  //for SPD only
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg_Eigen;  //for SPD only
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper, Eigen::IncompleteCholesky<double> > cg_Eigen;  //for SPD only
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> cg_Eigen;  //for SPD only
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> bicg_Eigen;  //BICGSTAB without preconditioner
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicg_Eigen;  //BICGSTAB with ILU preconditioner
	//Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> LU_Eigen;  //sparseLU method 


public:
	HO_2D() //default constructor
	{
		Knod = 1;
		dump_frequency = 1000;
		dt = 0.001;
		num_time_steps = 10000;
		problem_type = 10;
		time_integration_type = 1;
		HuynhSolver_type = 2;
		Reyn_inv = 0.001;
		Knod = 2;
	}; 
	HO_2D(const HO_2D& HO) :
		vorticity(HO.vorticity), stream_function(HO.stream_function) {}

	HO_2D& operator=(const HO_2D& HO)
	{
		if (this != &HO)
		{
			vorticity = HO.vorticity;
			stream_function = HO.stream_function;
		}
		return *this;
	}

	~HO_2D()
	{
		release_memory();
	}; // desctructor
	void save_output(int time);
	void release_memory();
	int read_input_file(const std::string filename);
	char allocate_arrays();
	char setup_mesh(); //based on the problem_type reads/creates the mesh
	void setup_sps_gps(); //setup the sps and their weights, gps
	void form_bases(); //form the sps lagrangian shepe function (&derivatives), gps shape functions (&derivatives), form Radau bases
	void setup_IC_BC_SRC(); //setup initial condition, boundary condition and the source/sink terms
	void form_metrics();
	unsigned int tensor2FEM(unsigned int i); // converts the tensor index to the FEM node ordering for 1D case
	unsigned int tensor2FEM(unsigned int i, unsigned int j); // converts the tensor index to the FEM node ordering for 2D case
	char solve_vorticity_streamfunction(); //solves the two set of equations: 1) advection diffucion of vorticity + 2) Poisson eq for streamfunction
	char solve_advection_diffusion(); //solves the advection diffusion equation for the vorticity field
	char Euler_time_integrate(); //use the explicit Euler method to integrate in time the vorticity evolution equation
	char solve_Poisson(); //solves the Poisson's equation for the streamfunction field
	char calc_RHS_advection(); //calculates the advective term on the RHS = -div(vw), v is Cartesian velocity vector, w is vorticity. stores in RHS_advective[el][Knod][Knod]
	char calc_RHS_diffusion();  //stores Laplace(w) in RHS_diffusive[el][Knod][Knod]
	void calc_internal_comm_vals_meth2(int el, int ijp, int ijpm, int ibnd, double* B_vort, double* B_vortm, double* B_G_vort, double* com_vort);
	void calc_boundary_comm_vals_meth2(int el, int el_b, int ijp, int ibnd, double* vort, double* Dvort, double* com_vort, double* com_grad_vort);
	void form_Laplace_operator_matrix(); //forms the left hand side matrix derived form the Laplace discretization. The matrix is sparse and in Eigen format
	void Poisson_solver_Hypre_setup(double*** laplacian_center, double**** laplacian_neighbor);  //setup and fill the LHS and RHS for Poisson solver via the Hypre 
	void Create_Hypre_Matrix();
	void Poisson_solver_Eigen_setup(double*** laplacian_center, double**** laplacian_neighbor);  //setup and fill the LHS and RHS for Poisson solver via the Eigen
};
