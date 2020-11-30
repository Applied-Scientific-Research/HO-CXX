#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "preprocess.h"

#define DirichletBC 0
#define NeumannBC 1

struct LR_boundary { // for the shape functions and Radua where everything reverts back to 1D and on the boundaries
	double left=1., right=1.;
};

cell_sides i2f[] = {west,east,south,north}; //(i2f[ijp] is the local side index (south, east, ... CCW) of the current cell: south =0:ibnd=0,idir=1; east=1:ibnd=1, idir=0, ...
unsigned int nbr[] = {east, south, west, north}; //ijpm = nbr[ijp]: neighbor of ijp face, returned in the same ijp ordering
double sgn[] = { -1., 1., -1., 1. }; //to calculate the upwind flux fronm the two neighboring fluxes at the common face

// the main class to hold the primitve and fundamental variables and methods
class HO_2D{

private:
	Mesh mesh;
	unsigned int Knod; //Number of sps in each direction,i.e. number of gauss point in each cell in each direction
	double g_prime[2]; //derivative of the correction function g on the [0]: left and [1]: right boundaries
	unsigned char* BC_switch; //specify if a boundary edge is dirichlet(0 or neumann(1) BC
	double** BC_values; // BC values for the Knod points of NelB boundary elements (in terms of velocity or gradient)
	double** BC_parl_vel, ** BC_normal_vel; //the parallel and normal contravariant flux velocity on all global boundary elements
	bool* BC_no_slip; //if the boundary conditions are no slip wall. if yes, then it is true
	double*** vorticity; //the vorticity field: first index goes for element, second and third go for j (eta dir) and i (csi dir) sps
	double*** initial_vorticity; //The initial (t=0) vorticity field
	double*** stream_function; //The stream function field
	Cmpnts2*** velocity_cart; //The velocity field
	double Reyn_inv; // the inverse of the Reynolds number
	unsigned char HuynhSolver_type; //based on Huynh's scheme type in Table 6.1 of his diffusion paper; types: 2 = stndard DG; 11 = higher order
	unsigned char time_integration_type; //time integration method; 1 = Euler; 2 = Mod.Euler; 4 = RK4
	unsigned int problem_type; //Solve different problems (types 1 and 2). prob_type=10 reads mesh from file
	unsigned int num_time_steps; //total number of time steps to march in time
	double dt; //time step size
	unsigned int dump_frequency; //the frequency of time saving
	bool fast; //0 for original, 1 for compact formulation (use gLB, gRB)
	unsigned int N_boundary; //number of boundary conditions (number of different boundaries)
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
	double*** RHS_advective; //the right hand side of the vorticity eq for the advective term: -div(vw) = -d_dxsi_j(V^jw)/G
	double*** RHS_diffusive;   //stores Laplace(w) in RHS_diffusive[el][Knod][Knod]
	double** velocity_jump;
	double*** Lap_vorticity; //Laplacian of the vorticity at [el][ky][kx]
	


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
		N_boundary = 4;
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

	void release_memory();
	int read_input_file(const std::string const filename);
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
	void calc_internal_comm_vals_meth2(int el, int ijp, int ibnd, double* B_vort, double* B_vortm, double* B_G_vort, double* com_vort);
	void calc_boundary_comm_vals_meth2(int el, int el_b, int ijp, int ibnd, double* vort, double* Dvort, double* com_vort, double* com_grad_vort);
};
