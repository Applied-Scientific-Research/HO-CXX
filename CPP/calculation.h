#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "preprocess.h"

struct LR_boundary { // for the shape functions and Radua where everything reverts back to 1D and on the boundaries
	double left, right;
};

// the main class to hold the primitve and fundamental variables and methods
class HO_2D{

private:
	Mesh mesh;
	unsigned int Knod; //Number of sps in each direction,i.e. number of gauss point in each cell in each direction
	double*** vorticity; //the vorticity field: first index goes for element, second and third go for j (eta dir) and i (csi dir) sps
	double*** initial_vorticity; //The initial (t=0) vorticity field
	double*** stream_function; //The stream function field
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
	LR_boundary* sps_boundary_basis; //The value Lagrange shape function of Knod solution points on the left[0](csi=-1) and right[1] (csi=+1) boundaries
	LR_boundary* sps_boundary_grad_basis; //The derivative of Lagrange shape function of Knod solution points on the left[0](csi=-1) and right[1] (csi=+1) boundaries
	


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



};
