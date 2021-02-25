#include "calculation.h"
#include "misc.hpp"
#include <cmath>

void HO_2D::release_memory() { //release the memory as destructor
	//delete[] vert_coor; delete[] vert_coor_initial; delete[] vert_coor_m1; delete[] Hex_verts_comps; delete[] projected_hex_verts_comps;
	//delete[] zone_id; delete[] group_regions; delete[] facet_groups; delete[] cell_N_verts; delete[] face_N_verts; delete[] N_connectivity;
	for (int k = 0; k < mesh.N_el; ++k) {
		for (int j = 0; j < Knod; ++j) {
			delete[] vorticity[k][j];
			delete[] stream_function[k][j];
			delete[] initial_vorticity[k][j];
		}
		delete[] vorticity[k];
		delete[] stream_function[k];
		delete[] initial_vorticity[k];
	}
	delete[] vorticity;
	delete[] stream_function;
	delete[] initial_vorticity;
}

int HO_2D::read_input_file(const std::string const filename) {
	//reads the essential solver / problem settings form the file filenme
	int retval=0;
	std::string temp_string;
	char temp_char;
	std::cout << "     Reading problem/solver settings file   ***** " << filename << " *****   opened for reading ..." << std::endl << std::endl;

	std::ifstream file_handle(filename);
	if (file_handle.fail()) {
		std::cout << "Input file opening failed.\n";
		exit(1);
	}

	try {
		//     retrieve number of quadrilateral elements in case of structured mesh (prespecified problem)
		getline(file_handle, temp_string);
		file_handle >> mesh.N_el_i >>temp_char>> mesh.N_el_j; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');   // ignore until next line
		// retrieve Knod: number of solution points in each element in each direction
		getline(file_handle, temp_string); 
		file_handle >> Knod; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read on Lnod_in: degree of the polynomial of the geometry edges: number of nodes on the elements edges -1
		getline(file_handle, temp_string);
		file_handle >> mesh.Lnod_in; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		mesh.Lnod = mesh.Lnod_in + 1; //number ofgeometry nodes on each edge
		// read the name of the mesh file and store it in the mesh.input_msh_file
		getline(file_handle, temp_string);
		getline(file_handle,temp_string); file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		mesh.input_msh_file = temp_string.c_str();
		//ignore the next 2 lines as they correspond to the exact geometry for concentric cylinders
		getline(file_handle, temp_string);
		// read in the Reynolds number and store the inverse of Reynolds number in Reyn_inv
		getline(file_handle, temp_string);
		file_handle >> Reyn_inv; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		Reyn_inv = Reyn_inv > 1.e-3 ? 1. / Reyn_inv : -1.; //if negative then ignore the vicous term
		// read inmultiplier to expand the grid geometrically by factor dx_ratio
		getline(file_handle, temp_string);
		file_handle >> mesh.dx_ratio; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in multiplier to make grids randomly non-uniform; uniform: fac=0.
		getline(file_handle, temp_string);
		file_handle >> mesh.fac; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in HuynhSolver_type
		getline(file_handle, temp_string);
		file_handle >> HuynhSolver_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in time integration method, 1 = Euler; 2 = Mod.Euler; 4 = RK4
		getline(file_handle, temp_string);
		file_handle >> time_integration_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in problem type
		getline(file_handle, temp_string);
		file_handle >> problem_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in total number of time steps
		getline(file_handle, temp_string);
		file_handle >> num_time_steps; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in time step size
		getline(file_handle, temp_string);
		file_handle >> dt; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in file saving frequency
		getline(file_handle, temp_string);
		file_handle >> dump_frequency; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in if it uses large stencil or compact stencil. fast-true means compact
		getline(file_handle, temp_string);
		file_handle >> fast; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// read in the type of matrix formation for poisson equation (1=Eigen, 2=Hypre)
		getline(file_handle, temp_string);
		file_handle >> LHS_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		// Later add to read in the type of solver and preconditioner. for now dont worry about it

		std::cout << "         Done" << std::endl;
		file_handle.close();
	}

	catch (...) {
		//		error.printError();
		retval = 2;
	}

	return retval;

}

char HO_2D::setup_mesh() {
	//based on the problem type (problem_type variable) creates the proper mesh. if problem_type is 10 then read from msh file
	char success;
	success = problem_type == 10 ? mesh.read_msh_file() : mesh.setup_mesh_problem(problem_type);
	return success;
}

void HO_2D::setup_sps_gps() {  
	//setup the local coordinate of the solution points (Gauss-Legendre) and their weights for integration, local coor of geometry points
	if (Knod == 1) {
		sps_local_coor[0] = 0.;
		sps_weight[0] = 2.;   //Gaussian weights
	}
	else if (Knod == 2) {
		sps_local_coor[1] = 1. / std::sqrt(3.);    //collocation nodes per elem are{ -C1, C1 }
		sps_local_coor[0] = -sps_local_coor[1];
		sps_weight[0] = sps_weight[1] = 1.0;  //Gaussian weights
	}
	else if (Knod == 3) {
		sps_local_coor[2] = std::sqrt(0.6);	//collocation nodes per elem are{ -C1, 0, C1 }
		sps_local_coor[1] = 0.;
		sps_local_coor[0] = -sps_local_coor[2];
		sps_weight[2] = sps_weight[0] = 5.0 / 9.0;
		sps_weight[1] = 8.0 / 9.0;
	}
	else if (Knod == 4) {
		sps_local_coor[3] = std::sqrt((15. + 2. * std::sqrt(30.)) / 35.);	//collocation nodes per elem are{ -C2, -C1, C1, C2 }
		sps_local_coor[2] = std::sqrt((15. - 2. * std::sqrt(30.)) / 35.);
		sps_local_coor[1] = -sps_local_coor[2];
		sps_local_coor[0] = -sps_local_coor[3];
		sps_weight[3] = sps_weight[0] = (18. - std::sqrt(30.)) / 36.;
		sps_weight[2] = sps_weight[1] = (18. + std::sqrt(30.)) / 36.;
	}
	else if (Knod == 5) {
		sps_local_coor[4] = std::sqrt(5. + 2. * std::sqrt(70.) / 7.) / 3.;     //collocation nodes per elem are{ -C2, -C1, 0, C1, C2 }
		sps_local_coor[3] = std::sqrt(5. - 2. * std::sqrt(70.) / 7.) / 3.;
		sps_local_coor[2] = 0.;
		sps_local_coor[1] = -sps_local_coor[3];
		sps_local_coor[0] = -sps_local_coor[4];
		sps_weight[0] = sps_weight[4] = (322. - 13. * std::sqrt(70.)) / 900.;
		sps_weight[1] = sps_weight[3] = (322. + 13. * std::sqrt(70.)) / 900.;
		sps_weight[2] = 128. / 225.;
	}
	else {
		std::cout << "The highest Knod allowed is 5!" << std :: endl;
		exit(2);
	}

	if (mesh.Lnod == 2) {//they are all equally spaced on the interval[-1, 1] including the end points - 1, 1
		gps_local_coor[1] = 1.;
		gps_local_coor[0] = -1.;
	}
	else if (mesh.Lnod == 3) {
		gps_local_coor[2] = 1.0;
		gps_local_coor[1] = 0.0;
		gps_local_coor[0] = -1.;
	}
	else if (mesh.Lnod == 4) {
		gps_local_coor[3] = 1.0;
		gps_local_coor[2] = 1./3.;
		gps_local_coor[1] = -1. / 3.;
		gps_local_coor[0] = -1.;
	}
	else {
		std::cout << "Only up to cubic elements (4 nodes on each edge) is allowd!" << std::endl;
		exit(3);
	}

	if (HuynhSolver_type == 2) { //g = g_DG
		g_prime[0] = -0.5 * Knod * Knod;
		g_prime[1] = 0.5 * Knod * Knod; //0.5 * minus_one_to_power(Knod) * Knod;
	}

}

char HO_2D::allocate_arrays() {
	//allocate the memory for the arrays and resize them.
	int N_el = mesh.N_el;
	int Lnod = mesh.Lnod;
	initial_vorticity = new double** [N_el];
	vorticity = new double** [N_el];
	stream_function = new double** [N_el];
	velocity_cart = new Cmpnts2** [N_el];
	for (int j = 0; j < N_el; ++j) {
		initial_vorticity[j] = new double* [Knod];
		vorticity[j] = new double* [Knod];
		stream_function[j] = new double* [Knod];
		velocity_cart[j] = new Cmpnts2* [Knod];
		for (int i = 0; i < Knod; ++i) {
			initial_vorticity[j][i] = new double[Knod];
			vorticity[j][i] = new double[Knod];
			stream_function[j][i] = new double[Knod];
			velocity_cart[j][i] = new Cmpnts2[Knod];
		}
	}

	sps_local_coor = new double[Knod];
	sps_weight = new double[Knod];
	gps_local_coor = new double[Lnod];
	sps_boundary_basis = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_boundary_basis[i] = new double[2]; //for left and right
	sps_boundary_grad_basis = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_boundary_grad_basis[i] = new double[2]; //for left and right
	sps_sps_grad_basis = new double* [Knod];
	for (int i = 0; i < Knod; ++i) sps_sps_grad_basis[i] = new double[Knod];

	gps_boundary_basis = new double* [Lnod]; for (int i = 0; i < Lnod; ++i) gps_boundary_basis[i] = new double[2];
	gps_boundary_grad_basis = new double* [Lnod]; for (int i = 0; i < Lnod; ++i) gps_boundary_grad_basis[i] = new double[2];
	gps_sps_basis = new double* [Lnod];
	gps_sps_grad_basis = new double* [Lnod];
	for (int i = 0; i < Lnod; ++i) gps_sps_basis[i] = new double[Knod];
	for (int i = 0; i < Lnod; ++i) gps_sps_grad_basis[i] = new double[Knod];

	sps_radau = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_radau[i] = new double[2]; //left[0] and right [1] radau functions on sps
	sps_grad_radau = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_grad_radau[i] = new double[2];

	vol_Dx_Dxsi = new Cmpnts2** [N_el];
	vol_Dy_Dxsi = new Cmpnts2** [N_el];
	G = new double**** [N_el];
	GB = new double**** [N_el];
	vol_jac = new double** [N_el];
	for (int i = 0; i < N_el; ++i) {
		vol_Dx_Dxsi[i] = new Cmpnts2* [Knod];
		vol_Dy_Dxsi[i] = new Cmpnts2* [Knod];
		vol_jac[i] = new double* [Knod];
		G[i] = new double*** [Knod];
		GB[i] = new double*** [2]; //for the direction d
		for (int j = 0; j < Knod; ++j) {
			vol_Dx_Dxsi[i][j] = new Cmpnts2 [Knod];
			vol_Dy_Dxsi[i][j] = new Cmpnts2 [Knod];
			vol_jac[i][j] = new double[Knod];
			G[i][j] = new double** [Knod];
			for (int r = 0; r < Knod; ++r) {
				G[i][j][r] = new double* [2];
				for (int s = 0; s < 2; ++s) G[i][j][r][s] = new double[2];
			}
		}
		for (int r = 0; r < 2; ++r) {
			GB[i][r] = new double** [2];
			for (int s = 0; s < 2; ++s) {
				GB[i][r][s] = new double* [Knod];
				for (int j = 0; j < Knod; ++j)
					GB[i][r][s][j] = new double[2];
			}
		}
	}

	face_Acoef = new double** [N_el];
	face_Bcoef = new double** [N_el];
	face_jac = new double** [N_el];
	face_Anorm = new double** [N_el];
	for (int i = 0; i < N_el; ++i) {
		face_Acoef [i] = new double* [4];
		face_Bcoef [i] = new double* [4];
		face_jac[i] = new double* [4];
		face_Anorm[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			face_Acoef[i][j] = new double [Knod];
			face_Bcoef[i][j] = new double[Knod];
			face_jac[i][j] = new double[Knod];
			face_Anorm[i][j] = new double[Knod];
		}

	}

	RHS_advective = new double** [mesh.N_el];
	for (int e = 0; e < mesh.N_el; ++e) {
		RHS_advective[e] = new double* [Knod];
		for (int i = 0; i < Knod; ++i)
			RHS_advective[e][i] = new double[Knod];
	}

	RHS_diffusive = new double** [N_el];
	for (int el = 0; el < N_el; ++el) {
		RHS_diffusive[el] = new double* [Knod];
		for (int k = 0; k < Knod; ++k)
			RHS_diffusive[el][k] = new double[Knod];
	}
	
	BC_no_slip = new bool[mesh.N_Gboundary]; //Set NoSlip to all walls (as default)
	boundary_source = new double** [mesh.N_edges_boundary];
	BC_switch = new unsigned char[mesh.N_edges_boundary];
	BC_psi = new double* [mesh.N_edges_boundary];
	BC_values = new double* [mesh.N_edges_boundary];
	velocity_jump = new double* [mesh.N_edges_boundary];
	BC_parl_vel = new double* [mesh.N_edges_boundary];
	BC_normal_vel = new double* [mesh.N_edges_boundary];
	for (int i = 0; i < mesh.N_edges_boundary; ++i) {
		boundary_source[i] = new double* [Knod*Knod];
		BC_psi[i] = new double[Knod];
		BC_values[i] = new double[Knod];
		velocity_jump[i] = new double[Knod];
		BC_parl_vel[i] = new double[Knod];
		BC_normal_vel[i] = new double[Knod];
		for (int j = 0; j < Knod*Knod; ++j) boundary_source[i][j] = new double[Knod]; //RHS for streamfunction Poisson's eq.
	}

	return 0;
}

void HO_2D::form_bases() {
	//This method froms the Lagrange shape functions: form the sps lagrangian shape function (&derivatives), gps shape functions (&derivatives), form Radau bases
	//--------------------------------------------------------------------------------------------------------------
	// *************************************** in This sub section ***************************************
	// calculate the value of shape functions of Knod solution points on the local left (-1) and right (+1) boundaries(sps_boundary_basis).
	// calculates the derivative of the Lagrange shape functions of Knod solution points on the local left (-1) and right (+1) boundaries(sps_boundary_grad_basis).
	// calculates the derivative of Knod solution points on the Knod solution points(sps_sps_grad_basis). the Knod sps is the first index(say i) and other nodes are second index(say j) : sps_sps_grad_basis(i, j)

	unsigned int Lnod = mesh.Lnod;
	double denominator;
	double grad_numerator[2];
	double* sps_grad_numerator = new double [Knod];

	for (int k = 0; k < Knod; ++k) {
		sps_boundary_basis[k][0] = sps_boundary_basis[k][1] = 1.;
		sps_boundary_grad_basis[k][0] = sps_boundary_grad_basis[k][1] = 0.0;
		for (int m = 0; m < Knod; ++m) sps_sps_grad_basis[k][m] = 0.;
		denominator = 1.;

		for (int j = 0; j < Knod; ++j) {
			if (j == k) continue;
			denominator *= sps_local_coor[k] - sps_local_coor[j];
			sps_boundary_basis[k][0] *= -1.0 - sps_local_coor[j]; //local left  boundary, i.e.xsi = -1
			sps_boundary_basis[k][1] *= 1.0 - sps_local_coor[j]; //local right  boundary, i.e.xsi = +1
			grad_numerator[0] = grad_numerator[1] = 1.;  //1 sentence in the numerator. The derivative has Knod-1 sentences. it has 2 indices 0:left boundary and 1 :right boundary
			for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] = 1.;

			for (int i = 0; i < Knod; ++i) {
				if (i == k || i == j) continue;
				grad_numerator[0] *= -1.0 - sps_local_coor[i];
				grad_numerator[1] *= 1.0 - sps_local_coor[i];
				for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] *= sps_local_coor[m] - sps_local_coor[i];
			}
			for (int LR=0; LR<2; ++LR) sps_boundary_grad_basis[k][LR] += grad_numerator[LR];

			for (int m = 0; m < Knod; ++m) sps_sps_grad_basis[k][m] += sps_grad_numerator[m];
		}
		
		for (int LR = 0; LR < 2; ++LR) {
			sps_boundary_basis[k][LR] /= denominator;
			sps_boundary_grad_basis[k][LR] /= denominator;
		}
		
		for (int m = 0; m < Knod; ++m) sps_sps_grad_basis[k][m] /= denominator;
	}

	// *************************************** in This sub section ***************************************
	// calculates the value of shape function of Lnod gps on the 2 boundaries(left and right)(gps_boundary_basis).
	// calculates the shape function of Lnod gps on the Knod sps so the size is(Lnod, Knod)(gps_sps_basis).
	// calculates the derivative of shape functions of Lnod gps on the left and right boundaries, so the size is(Lnod, 2) (gps_boundary_grad_basis)
	// calculates the derivative of shape function of Lnod gps on the Knod sps so the size is [Lnod][Knod] (gps_sps_grad_basis).

	for (int k = 0; k <Lnod; ++k) {
		gps_boundary_basis[k][0] = gps_boundary_basis[k][1] = 1.;
		gps_boundary_grad_basis[k][0] = gps_boundary_grad_basis[k][1] = 0.;
		for (int m = 0; m < Knod; ++m) gps_sps_grad_basis[k][m] = 0.0;
		for (int m = 0; m < Knod; ++m) gps_sps_basis[k][m] = 1;

		denominator = 1.;
		for (int j = 0; j < Lnod; ++j) {
			if (j == k) continue;
			denominator *= gps_local_coor[k] - gps_local_coor[j];
			gps_boundary_basis[k][1] *= 1. - gps_local_coor[j]; //right boundary
			gps_boundary_basis[k][0] *= -1. - gps_local_coor[j];  //left boundary
			for (int m = 0; m < Knod; ++m) gps_sps_basis[k][m] *= sps_local_coor[m] - gps_local_coor[j];
			grad_numerator[0] = grad_numerator[1] = 1.;
			for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] = 1.;

			for (int i = 0; i < Lnod; ++i) {
				if (i == k || i == j) continue;
				grad_numerator[0] *= -1.0 - gps_local_coor[i];
				grad_numerator[1] *= 1.0 - gps_local_coor[i];
				for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] *= sps_local_coor[m] - gps_local_coor[i];
			}
			for (int LR=0; LR<2; ++LR) gps_boundary_grad_basis[k][LR] += grad_numerator[LR];

			for (int m = 0; m < Knod; ++m) gps_sps_grad_basis[k][m] += sps_grad_numerator[m];
		}

		for (int LR = 0; LR < 2; ++LR) gps_boundary_basis[k][LR] /= denominator;
		for (int LR = 0; LR < 2; ++LR) gps_boundary_grad_basis[k][LR] /= denominator;
		for (int m = 0; m < Knod; ++m) gps_sps_grad_basis[k][m] /= denominator;
		for (int m = 0; m < Knod; ++m) gps_sps_basis[k][m] /= denominator;
	}

	// *************************************** In this sub section ***************************************
	// calculates and stores the value of right(sps_radau[k][0]) and left(sps_radau[k][1]) radau functions on the Knod sps and stores in sps_radau[k][0:1]
	// calculates and stores the derivative of right(sps_grad_radau[k][0]) and left(sps_grad_radau[k][1]) radau functions on the Knod sps and stores in sps_grad_radau[k][0:1]

	double coef, coefd;

		//		R_k(x) = (-1) ^ k / 2 * (P(x)_k - P(x)_(k - 1)), where P(x)_k is Legendre polynomial of order k, 3.15 in Hyun diffusion paper
		//		derivative :D[R_k(x), x] = (k / 2) * (-1) ^ k * (P(x)_k + P(x)_(k - 1)) / (x + 1)
	if (!Knod) {
		std::cout << "ERROR: Order 0 Radau is Undefined!" << std::endl;
		exit(1);
	}
	else if (Knod == 1) {
		sps_radau[0][0] = sps_radau[0][1] = 0.5;
		sps_grad_radau[0][1] = 0.5;
		sps_grad_radau[0][0] = -0.5;
	}
	else {
		coef =  0.5 * minus_one_to_power(Knod);
		coefd = 0.5 * Knod * minus_one_to_power(Knod);
		for (int i = 0; i < Knod; ++i) {
			//	remember right radau(xsi) = left radau(-xsi) : R_k(x) | Right = R_k(-x) | Left
			sps_radau[i][0] = coef * (Legendre(Knod, sps_local_coor[i]) - Legendre(Knod-1, sps_local_coor[i])); //value of right Radau function(g_LB) on the sps
			sps_radau[i][1]  = coef * (Legendre(Knod, -sps_local_coor[i]) - Legendre(Knod - 1, -sps_local_coor[i])); //value of left Radau function(g_LB) on the sps; 

			// D[R_k(x), x] | Right = -D[R_k(-x), x] | Left
			sps_grad_radau[i][0] = coefd * (Legendre(Knod, sps_local_coor[i]) + Legendre(Knod - 1, sps_local_coor[i])) / (1. + sps_local_coor[i]);
			sps_grad_radau[i][1] = -coefd * (Legendre(Knod, -sps_local_coor[i]) + Legendre(Knod - 1, -sps_local_coor[i])) / (1. - sps_local_coor[i]);
		}		
	}
	delete[] sps_grad_numerator;
}

void HO_2D::setup_IC_BC_SRC() {
	//setup initial condition, boundary condition and the source/sink terms	
	for (int el = 0; el < mesh.N_el; ++el)
		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i) {
				initial_vorticity[el][j][i] = 0.;
				vorticity[el][j][i] = initial_vorticity[el][j][i];
			}

	for (int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b)
		for (int m = 0; m < Knod; ++m) {
			BC_normal_vel[el_b][m] = 0.;
			BC_values[el_b][m] = 0.;
			BC_parl_vel[el_b][m] = 0.;
		}
			
	std::fill(BC_switch, BC_switch + mesh.N_edges_boundary, NeumannBC);
	std::fill(BC_no_slip, BC_no_slip + mesh.N_Gboundary, true);
	for (unsigned int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b) std::fill(BC_psi[el_b], BC_psi[el_b] + Knod, 0.);

	if (problem_type == 1) {
		// Periodic moving wall problem; bottom velocity is 1; top is 0
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) BC_switch[el_b] = DirichletBC;
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) {
			if (std::fabs(mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y) < 1.e-6 && std::fabs(mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y) < 1.e-6) //boundary on the bottom
				for (int m = 0; m < Knod; ++m) BC_parl_vel[el_b][m] = 1.0;    //Wall parallel velocity
		}
	}

	else if (problem_type == 3) {
		//Square Cavity Problem; top velocity is 1
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) BC_switch[el_b] = DirichletBC;
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) {
			if (std::fabs(mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y - 1.) < 1.e-6 && std::fabs(mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y - 1.) < 1.e-6) //top boundary
				for (int m = 0; m < Knod; ++m) BC_parl_vel[el_b][m] = 1.0; //top wall parallel velocity
		}
	}

	else if (problem_type == 10) {
		/*
		Need to make BC_psi generic to account for width of inlet automatically instead of the current manual setup
        We want a uniform inlet velocity of 1, so psi = y / (2*HalfHeight)
        But we want the values at the wall solution points, so y = 0.5 * ( (y_2 + y_1) + sps * (y_2 - y_1) )
        where y_2 and y_1 are the coords of the end points of a boundary element and sps is the solution point position in local coords
		*/

		std::fill(BC_switch, BC_switch + mesh.N_edges_boundary, DirichletBC);
		double y2=0., y1=0.;
		for (int Gboundary = 0; Gboundary < mesh.N_Gboundary; ++Gboundary) {
			if (mesh.boundaries[Gboundary].name == "top") {
				int edge_index = mesh.boundaries[Gboundary].edges[0]; //the index of the first edge forming the top boundary
				int node_index = mesh.edges[edge_index].nodes[0];
				y2 = mesh.nodes[node_index].coor.y; // the y coordinate of the first node on the first edge on the top boundary
			}
			else if (mesh.boundaries[Gboundary].name == "bottom") {
				int edge_index = mesh.boundaries[Gboundary].edges[0]; //the index of the first edge forming the bottom boundary
				int node_index = mesh.edges[edge_index].nodes[0];
				y1 = mesh.nodes[node_index].coor.y; // the y coordinate of the first node on the first edge on the bottom boundary
			}
		}
		double half_height = 0.5 * std::fabs(y2 - y1);

		for (int Gboundary = 0; Gboundary < mesh.N_Gboundary; ++Gboundary) {
			if (mesh.boundaries[Gboundary].name == "top") {
				// Set slip BC to not allow diffusion
				BC_no_slip[Gboundary] = false;
				for (unsigned int el_b : mesh.boundaries[Gboundary].edges)
					std::fill(BC_psi[el_b], BC_psi[el_b] + Knod, half_height); //all the edges on the top boundary are a streamline matching the streamline at the top side of the inlet
			}

			else if (mesh.boundaries[Gboundary].name == "bottom") {
				BC_no_slip[Gboundary] = false;
				for (unsigned int el_b : mesh.boundaries[Gboundary].edges)
					std::fill(BC_psi[el_b], BC_psi[el_b] + Knod, -half_height); //all the edges on the top boundary are a streamline matching the streamline at the top side of the inlet
			}

			else if (mesh.boundaries[Gboundary].name == "inlet") {
				// Set slip BC to not allow diffusion
				BC_no_slip[Gboundary] = false;
				for (unsigned int el_b : mesh.boundaries[Gboundary].edges) {
					for (int ky = 0; ky < Knod; ++ky) {
						//interpolate y - coord of sol pt at(ky) using nodal coordinates of element j
						double y = 0.;
						for (int j = 0; j < mesh.Lnod; ++j) {// on edge j is 0 1 2 3 4
							int node_index = mesh.edges[el_b].nodes[tensor2FEM(j)]; // retrives 0 2 3 1 on edge (local numbering)
							y += gps_sps_basis[j][ky] * mesh.nodes[node_index].coor.y; // the y coordinate of the sps ky
						}
						BC_psi[el_b][ky] = y;
						BC_normal_vel[el_b][ky] = 1.0;
					} //for ky
				} //for edges on the "inlet" boundary
			}

			else if (mesh.boundaries[Gboundary].name == "outlet") {
				BC_no_slip[Gboundary] = false;
				for (unsigned int el_b : mesh.boundaries[Gboundary].edges) {
					for (int ky = 0; ky < Knod; ++ky) {
						double y = 0.;
						for (int j = 0; j < mesh.Lnod; ++j) {// on edge j is 0 1 2 3 4
							int node_index = mesh.edges[el_b].nodes[tensor2FEM(j)]; // retrives 0 2 3 1 on edge (local numbering)
							y += gps_sps_basis[j][ky] * mesh.nodes[node_index].coor.y; // the y coordinate of the sps ky
						}
						BC_psi[el_b][ky] = y;
						BC_normal_vel[el_b][ky] = 1.0;
					} //for ky
				} //for edges on the "inlet" boundary
			}
			else {
				std::cout << "The type of Boundary condition is not specified in specific forms, check the BC" << std::endl;
				exit(1);
			}
		} //for Gboundary
	}  //if problem_tpe==10
}

void HO_2D::form_metrics() {
	/* forms
	vol_Dx_Dxsi; //the derivative of x wrt to xsi_s(s is 0, 1 to make dx / dxsi, dx / deta) at the sps(jx, jy) in tensor form on element el.This forms Vol_Dx_iDxsi_j(el, jy, jx).x, .y
	vol_Dy_Dxsi; //the derivative of y wrt to xsi_s(s is 0, 1 to make dy / dxsi, dy / deta) at the sps(jx, jy) in tensor form on element el.This forms Vol_Dy_iDxsi_j(el, jy, jx).x, .y		
	vol_jac[el][jy][jx] : is the cross product of vectors g_1 = (dx / dcsi, dy / dcsi) x g_2 = (dx / deta, dy / deta) at the sps(jx, jy) of element el.g_1 and g_2 are based on Fotis notes(covariant bases).So, the cross product is G = det[dx / dcsi, dy / dxsi; dx / deta, dy / deta] = ratio of volume of original element / volume of transformed element(del_csi(2) * del_eta(2) = 4)
		!Haji : Face_Jac(i(1:Knod), r(0:3), el(1:Nel)) is the G = dx / dxsi * dy / deta - dx / deta * dy / dxsi on the boundary face at the r = 0(left) or r = 1 (right) or r = 2(south) or r = 3(north)face on the i sps of the element el
	Face_Acoef[el][ijp][Knod] is the Case 1) squared length of g_2 = (dx / deta, dy / deta), i.e.g_2(dot)g_2, on the leftand right boundaries(ijp=0,1); Case 2) squared length of g_1 = (dx / dcsi, dy / dcsi), i.e.g_1(dot)g_1, on the bottom and top boundaries(ijp= 2, 3)
		!Haji: Face_Bcoef(i(1:Knod), r(0:3), el(1:Nel)) is the dot product of g_1 with g_2 = g_1(dot)g_2 = dx / dxsi * dx / deta + dy / dxsi * dy / deta) on the 4 boundaries(left, right, bottom, top(r = 0, 1, 2, 3)
	Face_ANorm[el][ijp][Knod] is the norm of Face_Acoef, i.e.Face_Norm[el][0:3][k] = Sqrt(Face_Acoef[el][0:3][k]
			!Haji : The Face_Acoef and Face_Bcoef, forms the covariant metric tensor g_ij = [g11 g12; g12, g22]] (symmetric matrix))\

	G[el][j][i]][d1][d2] is (g_d1 dot g_d2)/J at j,i of element el, so the size is[N_el][Knod][[Knod][2][2]. here g_0=(-partial(x)/partial(eta),partial(y)/partial(eta)), g_1=(partial(x)/partial(csi),-partial(y)/partial(csi))
	GB[el][d][t]][j][d2] is (g_d dot g_d2)/J at the jth flux point on the d direction , t side of the element el, so the size is[N_el][2][[2][Knod][2]. here g_0=(-partial(x)/partial(eta),partial(y)/partial(eta)), g_1=(partial(x)/partial(csi),-partial(y)/partial(csi)) on the face

	*/


	unsigned int Lnod = mesh.Lnod;
	double dx_dxsi, dy_dxsi, dx_deta, dy_deta;
	Cmpnts2 g[2]; //to store g_0=(-partial(x)/partial(eta),partial(y)/partial(eta)) and g_1=(partial(x)/partial(csi),-partial(y)/partial(csi))
	Cmpnts2** local_coor = new Cmpnts2 * [Lnod]; //local array tp store the coor of the gps in an element, dont really need it just for convenience
	for (int i = 0; i < Lnod; ++i) local_coor[i] = new Cmpnts2 [Lnod];

	for (int el = 0; el < mesh.N_el; ++el) {
		for (int j = 0; j < Lnod; ++j)
			for (int i = 0; i < Lnod; ++i)
				local_coor[j][i] = mesh.nodes[mesh.elements[el].nodes[tensor2FEM(i, j)]].coor;

		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i) {
				for (int n = 0; n < Lnod; ++n)
					for (int m = 0; m < Lnod; ++m) {
						vol_Dx_Dxsi[el][j][i].x += gps_sps_grad_basis[m][i] * gps_sps_basis[n][j] * local_coor[n][m].x;	//grad at x - dir  * no grad at y - dir  * x / y - coord of geom
						vol_Dy_Dxsi[el][j][i].x += gps_sps_grad_basis[m][i] * gps_sps_basis[n][j] * local_coor[n][m].y;	//grad at x - dir  * no grad at y - dir  * x / y - coord of geom

						vol_Dx_Dxsi[el][j][i].y += gps_sps_grad_basis[n][j] * gps_sps_basis[m][i] * local_coor[n][m].x;	//grad at y - dir  * no grad at x - dir  * x / y - coord of geom
						vol_Dy_Dxsi[el][j][i].y += gps_sps_grad_basis[n][j] * gps_sps_basis[m][i] * local_coor[n][m].y;	//grad at y - dir  * no grad at x - dir  * x / y - coord of geom
					}
				vol_jac[el][j][i] = vol_Dx_Dxsi[el][j][i].x * vol_Dy_Dxsi[el][j][i].y - vol_Dx_Dxsi[el][j][i].y * vol_Dy_Dxsi[el][j][i].x;
				g[0].set_coor(-vol_Dx_Dxsi[el][j][i].y, vol_Dy_Dxsi[el][j][i].y);
				g[1].set_coor(vol_Dx_Dxsi[el][j][i].x, -vol_Dy_Dxsi[el][j][i].x);
				for (int r = 0; r < 2; ++r)
					for (int s = 0; s < 2; ++s) 
						G[el][j][i][r][s] = DOT(g[r], g[s])/vol_jac[el][j][i];
			}


		for (int k = 0; k < Knod; ++k) { //loop on all sps (flux points) on the faces to calculate metrics
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			// ****** sps on the left (west) boundary (xsi=-1) *********
			for (int j=0; j<Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) { 
					dx_dxsi += gps_boundary_grad_basis[i][0] * gps_sps_basis[j][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_grad_basis[i][0] * gps_sps_basis[j][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_basis[i][0] * gps_sps_grad_basis[j][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_basis[i][0] * gps_sps_grad_basis[j][k] * local_coor[j][i].y;
				}
			// note: On lines of constant xsi, it is proven in page 6 of my curvilinear notes that dS = g_2 d_eta, so ||dS|| = ||g_2|| * d_eta, dS is element length
			face_jac[el][0][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][0][k] = dx_deta * dx_deta + dy_deta * dy_deta; //g_2 dot g_2 in fotis notes = g_{22}
			face_Bcoef[el][0][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_Anorm[el][0][k] = std::sqrt(face_Acoef[el][0][k]); //||g_2||
			g[0].set_coor(-dx_deta, dy_deta);
			g[1].set_coor(dx_dxsi, -dy_dxsi);
			for (int r = 0; r < 2; ++r)
				GB[el][0][0][k][r] = DOT(g[0], g[r]) / face_jac[el][0][k];

			// **********************************************************
			// ****** sps on the right (east) boundary (xsi=+1) *********
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			for (int j = 0; j < Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) {
					dx_dxsi += gps_boundary_grad_basis[i][1] * gps_sps_basis[j][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_grad_basis[i][1] * gps_sps_basis[j][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_basis[i][1] * gps_sps_grad_basis[j][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_basis[i][1] * gps_sps_grad_basis[j][k] * local_coor[j][i].y;
				}
			face_jac[el][1][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][1][k] = dx_deta * dx_deta + dy_deta * dy_deta; //g_2 dot g_2 in fotis notes = g_{22}
			face_Bcoef[el][1][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_Anorm[el][1][k] = std::sqrt(face_Acoef[el][1][k]); //||g_2||
			g[0].set_coor(-dx_deta, dy_deta);
			g[1].set_coor(dx_dxsi, -dy_dxsi);
			for (int r = 0; r < 2; ++r)
				GB[el][0][1][k][r] = DOT(g[0], g[r]) / face_jac[el][1][k];
			// ************************************************************
			// ****** sps on the bottom (south) boundary (eta=-1) *********
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			for (int j = 0; j < Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) {
					dx_dxsi += gps_boundary_basis[j][0] * gps_sps_grad_basis[i][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_basis[j][0] * gps_sps_grad_basis[i][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_grad_basis[j][0] * gps_sps_basis[i][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_grad_basis[j][0] * gps_sps_basis[i][k] * local_coor[j][i].y;
				}
			face_jac[el][2][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][2][k] = dx_dxsi * dx_dxsi + dy_dxsi * dy_dxsi; //g_1 dot g_1 in fotis notes = g_{11}
			face_Bcoef[el][2][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_Anorm[el][2][k] = std::sqrt(face_Acoef[el][2][k]); //||g_1||
			g[0].set_coor(-dx_deta, dy_deta);
			g[1].set_coor(dx_dxsi, -dy_dxsi);
			for (int r = 0; r < 2; ++r)
				GB[el][1][0][k][r] = DOT(g[1], g[r]) / face_jac[el][2][k];
			// ************************************************************
			// ****** sps on the top (north) boundary (eta=+1) *********
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			for (int j = 0; j < Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) {
					dx_dxsi += gps_boundary_basis[j][1] * gps_sps_grad_basis[i][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_basis[j][1] * gps_sps_grad_basis[i][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_grad_basis[j][1] * gps_sps_basis[i][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_grad_basis[j][1] * gps_sps_basis[i][k] * local_coor[j][i].y;
				}
			face_jac[el][3][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][3][k] = dx_dxsi * dx_dxsi + dy_dxsi * dy_dxsi; //g_1 dot g_1 in fotis notes = g_{11}
			face_Bcoef[el][3][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_Anorm[el][3][k] = std::sqrt(face_Acoef[el][3][k]); //||g_1||
			g[0].set_coor(-dx_deta, dy_deta);
			g[1].set_coor(dx_dxsi, -dy_dxsi);
			for (int r = 0; r < 2; ++r)
				GB[el][1][1][k][r] = DOT(g[1], g[r]) / face_jac[el][3][k];
		}  //for k=0; k<Knod
	} //for el
	
	for (int i = 0; i < Lnod; ++i)
		delete[] local_coor[i];
	delete[] local_coor;
	
}

unsigned int HO_2D::tensor2FEM(unsigned int i) {
	// converts the tensor index to the FEM node ordering for 1D case
	if (!i) return 0;
	else if (i == mesh.Lnod_in) return 1;
	else return (i + 1);
}

unsigned int HO_2D::tensor2FEM(unsigned int i, unsigned int j) {
	// converts the tensor index to the FEM node ordering for 1D case
	unsigned  t2f;
	unsigned int Lnod_in = mesh.Lnod_in;

	if (Lnod_in == 1) t2f = i+3*j-2*i*j;   //linear element
		
	else if (Lnod_in == 2) {
		if (!j) { //first row
			if (!i) t2f = 0;
			else if (i == 1) t2f = 4;
			else t2f = 1;
		}
		else if (j == 1) {
			if (!i) t2f = 7;
			else if (i == 1) t2f = 8;
			else t2f = 5;
		}
		else {
			if (!i) t2f = 3;
			else if (i == 1) t2f = 6;
			else t2f = 2;
		}
	}

	else if (Lnod_in == 3) {
		if (!j) { //first row
			if (!i) t2f = 0;
			else if (i == 1) t2f = 4;
			else if (i == 2) t2f = 5;
			else t2f = 1;
		}
		else if (j == 1) {
			if (!i) t2f = 11;
			else if (i == 1) t2f = 12;
			else if (i == 2) t2f = 13;
			else t2f = 6;
		}

		else if (j == 2) {
			if (!i) t2f = 10;
			else if (i == 1) t2f = 15;
			else if (i == 2) t2f = 14;
			else t2f = 7;
		}

		else {
			if (!i) t2f = 3;
			else if (i == 1) t2f = 9;
			else if (i == 2) t2f = 8;
			else t2f = 2;
		}
	}
	else {
		std::cout << "Only up to cubic finite elements (Lnod = 4) are supported!" << std::endl;
		exit(2);
	}
	return t2f;
}

char HO_2D::solve_vorticity_streamfunction() {
	//solves the two set of equations: 1) advection diffucion of vorticity + 2) Poisson eq for streamfunction
	for (int el = 0; el<mesh.N_el; ++el)
		for (int j = 0; j < Knod; ++j) 
			for (int i = 0; i < Knod; ++i)
				vorticity[el][j][i] = initial_vorticity[el][j][i];
		

	save_output(0);
	form_Laplace_operator_matrix(); //form the LHS matrix in Laplace discretization in Eigen format. Done only once

	for (unsigned int ti = 1; ti <= num_time_steps; ++ti) {
		std::cout << "timestep  " << ti << std::endl;
		
		//solve_advection_diffusion();

		//temporarily set the right hand side (vorticity) of the poisson equation
		for (int el = 0; el < mesh.N_el; ++el) {
			for (int j = 0; j < Knod; ++j) {
				for (int i = 0; i < Knod; ++i) {
					Cmpnts2 coor(0., 0.); //coordinate of sps[j][i]
					for (int ny = 0; ny < mesh.Lnod; ++ny)
						for (int nx = 0; nx < mesh.Lnod; ++nx) {
							int node_index = mesh.elements[el].nodes[tensor2FEM(nx, ny)];
							coor.plus(gps_sps_basis[ny][j] * gps_sps_basis[nx][i], mesh.nodes[node_index].coor);
						}
					double ym = M_PI * (1. - coor.y);
					//vorticity[el][j][i] = M_PI * std::sin(M_PI * coor.x) * (ym * (1. + 2. * std::cos(ym)) + 2. * std::sin(ym));
					//vorticity[el][j][i] = std::sin(M_PI * coor.y);
					//vorticity[el][j][i] = M_PI*std::sin(M_PI*coor.x)*(2.*std::sin(M_PI*coor.y) + 2.*M_PI*coor.y*std::cos(M_PI*coor.y) + M_PI*coor.y);
					vorticity[el][j][i] = -4. * (coor.x * coor.x + coor.y * coor.y + 1.) * std::exp(coor.x * coor.x + coor.y * coor.y);
					//vorticity[el][j][i] = -4.;
				}
			}
		}
		solve_Poisson();

		if (!(ti % dump_frequency) || ti==1)
			save_output(ti);
 	}

	return 0;
}

char HO_2D::solve_advection_diffusion() {
	//solves the advection diffusion eq. for the vorticity field
	if (time_integration_type == 1) { // first order Euler time integration method
		Euler_time_integrate(); //calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
	}
	return 0;
}

void HO_2D:: form_Laplace_operator_matrix() {
	// This subroutine forms the left hand side matrix derived form the Laplace discretization. The matrix is sparse and in Eigen format
	int N_el = mesh.N_el;
	int a, b, q, Mq; //temporary indices as in my notes
	int eln, ijp, ijpm; //neighbor element; element side, neighbor element side
	int dn, tn; // the direction and boundary side of the neighboring element eln, adjacent to d direction and t side of the element el
	int rs, ij, i1, j1; // temporary indices
	double coeff, tmp2, tmp3;
	unsigned int Ksq = Knod * Knod, Km1 = Knod - 1, K4 = Ksq * Ksq;
	double*** laplacian_center = new double** [N_el]; //The coefficients in the left hand side matrix that has contribution from the element itself. it has the coefficients for element el, sps ij=j*Knod+i and rs=r*Knod+s, so laplacian_center[el][ij][rs]
	double**** laplacian_neighbor = new double*** [N_el]; //The coefficients in the left hand side matrix that has contribution from the 4 neighboring element (if are not located on the global boundary). it has the coefficients for element el per cell side, sps ij=j*Knod+i and rs=r*Knod+s of the neighbor element, so laplacian_center[el][4][ij][rs]
	//double** NeuMatrix_Orig = new double* [Knod], **NeuMatrix = new double* [Knod]; //small dense matrix to obtain comx(per element) for a given Neumann BC
	Eigen::MatrixXd NeuMatrix_Orig(Knod, Knod), NeuMatrix(Knod, Knod); //I am using Eigen here to save energy and time

	std::fill(BC_switch, BC_switch + mesh.N_edges_boundary, /*NeumannBC*/ DirichletBC); //to test the poisson solver
	for (int el_b = 0; el_b<mesh.N_edges_boundary; ++el_b) {
		double x1 = mesh.nodes[mesh.edges[el_b].nodes[0]].coor.x;
		double x2 = mesh.nodes[mesh.edges[el_b].nodes[1]].coor.x;
		double y1 = mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y;
		double y2 = mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y;

		if (std::sqrt(x1*x1+y1*y1)<0.5101 && std::sqrt(x2 * x2 + y2 * y2) < 0.501) BC_switch[el_b] = NeumannBC;

		/*
		if (mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y < 1e-6 && mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y < 1e-6) BC_switch[el_b] = NeumannBC;
		else if (mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y > 9.99e-1 && mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y > 9.99e-1) BC_switch[el_b] = DirichletBC;
		else if (mesh.nodes[mesh.edges[el_b].nodes[0]].coor.x < 1e-6 && mesh.nodes[mesh.edges[el_b].nodes[1]].coor.x < 1e-6) BC_switch[el_b] = NeumannBC;
		else BC_switch[el_b] = NeumannBC;
		*/
	}


	for (int i = 0; i < N_el; ++i) {
		laplacian_center[i] = new double* [Ksq];
		for (int j = 0; j < Ksq; ++j)
			laplacian_center[i][j] = new double[Ksq];
	}
	for (int i = 0; i < N_el; ++i) {
		laplacian_neighbor[i] = new double** [4];
		for (int j = 0; j < 4; ++j) {
			laplacian_neighbor[i][j] = new double* [Ksq];
			for (int k = 0; k < Ksq; ++k) laplacian_neighbor[i][j][k] = new double [Ksq];
		}
	}

	for (int el = 0; el < N_el; ++el) {
		for (int ij = 0; ij < Ksq; ++ij) {
			for (int rs = 0; rs < Ksq; ++rs) {
				laplacian_center[el][ij][rs] = 0.;
				for (int side = 0; side < 4; ++side) laplacian_neighbor[el][side][ij][rs] = 0.;
			}
		}
	}

	for (int el_b=0; el_b<mesh.N_edges_boundary; ++el_b)
		for (int ij = 0; ij < Ksq; ++ij)
			for (int i = 0; i < Knod; ++i)
				boundary_source[el_b][ij][i] = 0.;

	nnz = K4 * N_el;
	for (int el = 0; el < N_el; ++el) {
		// **************** effect of the all-cases term (unrelated to the element boundaries)*********************
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				ij = j * Knod + i; //local cumulative sps index in element k

				for (int m = 0; m < Knod; ++m) {
					for (int r = 0; r < 2; r++) {
						a = m * (1 - r) + i * r;
						b = m * r + j * (1 - r);
						tmp3 = 0.;
						for (int s = 0; s < 2; s++) tmp3 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s]; //from the M,N,P,Q of ADrins note
						double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
						for (int d = 0; d < 2; ++d) {
							double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
							int tmp1 = a * (1 - d) + b * d; //the second index for SNGLB
							for (int p = 0; p < Knod; ++p) {
								j1 = p * d + b * (1 - d);
								i1 = p * (1 - d) + a * d;
								rs = j1 * Knod + i1;
								tmp2 = 0.;
								for (int t = 0; t < 2; ++t) tmp2 += 0.5 * sps_boundary_basis[p][t] * sps_grad_radau[tmp1][t];
								laplacian_center[el][ij][rs] += A * SNGLB1 * (sps_sps_grad_basis[p][tmp1] - tmp2);  //from the L term on page 7 of Adrin's note; page 1 of my notes, 3 feb 2021
								laplacian_center[el][ij][rs] +=  A * tmp3 * (-sps_sps_grad_basis[p][tmp1] + tmp2); //from the M,N,P,Q of Adrin's note: page 2 of my notes 3 feb 2021
							}
						}
					}
				}
			}
		}
		// **********************************************************************************************************

		// effect of the elements faces
		for (int d = 0; d < 2; ++d) {
			for (int t = 0; t < 2; ++t) {
				ijp = 2 * d + t;
				cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir

				if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) { // if this face is on the global boundary:
					int edge_index = mesh.elem_neighbor[el].boundary_index[elem_side];  // the index of the edge that is on global boundary, which is on the d direction and t side of element el
					if (BC_switch[edge_index] == DirichletBC) { // Dirichlet boundary condition
						//******************************* The dirichlet BC case: *********************************
						for (int j = 0; j < Knod; ++j) {
							for (int i = 0; i < Knod; ++i) {
								ij = j * Knod + i; //local cumulative sps index in element k
								for (int m = 0; m < Knod; ++m) {
									for (int r = 0; r < 2; r++) {
										a = m * (1 - r) + i * r;
										b = m * r + j * (1 - r);
										q = a * d + b * (1 - d);
										double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
										double NGR1 = sps_grad_radau[a * (1 - d) + b * d][t];
										double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
										tmp2 = 0.;
										for (int s = 0; s < 2; s++) tmp2 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s];
										for (int p = 0; p < Knod; ++p) {
											j1 = p * d + b * (1 - d);
											i1 = p * (1 - d) + a * d;
											rs = j1 * Knod + i1; //cumulative index
											laplacian_center[el][ij][rs] -= 0.5 * SNGLB1 * A * sps_boundary_basis[p][t] * NGR1;  //page 1 of my notes 3 Feb 2021
											laplacian_center[el][ij][rs] += 0.5 * tmp2 * A * sps_boundary_basis[p][t] * NGR1;  //page 2 of my notes 3 Feb 2021
										}
										boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * q - 1) + q] += SNGLB1 * A * NGR1; //page 1 on my notes 3 Feb 2021, effect of DBC, according to the edge[edge_index] direction
										boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * q - 1) + q] -= tmp2 * A * NGR1; //page 2 on my notes 3 Feb 2021, effect of DBC, according to the edge[edge_index] direction

									}
								}
							}
						}

						for (int j = 0; j < Knod; ++j) {
							for (int i = 0; i < Knod; ++i) {
								ij = j * Knod + i; //local cumulative sps index in element k
								int alpha = i * d + j * (1 - d);
								double A = GB[el][d][t][alpha][d];  //[(g_d . g_d)/J]|k,d,t,alpha
								double NGR1 = sps_grad_radau[j * d + i * (1 - d)][t];
								for (int m = 0; m < Knod; ++m) {
									j1 = m * d + j * (1 - d);
									i1 = m * (1 - d) + i * d;
									rs = j1 * Knod + i1;
									laplacian_center[el][ij][rs] += A * (sps_boundary_grad_basis[m][t] - g_prime[t] * sps_boundary_basis[m][t]) * NGR1; //page 3 of my notes 3 Feb 2021
								}
								boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * alpha - 1) + alpha] += A * g_prime[t] * NGR1; //page 3 of my notes 3 Feb 2021
								double B = GB[el][d][t][alpha][1 - d];
								for (int m = 0; m < Knod; ++m) {
									boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * m - 1) + m] += B * sps_sps_grad_basis[m][alpha] * NGR1; //page 3 of my notes 3 Feb 2021
								}
							}
						}
					} // if the edge is located on the Dirichlet-type global boundary

					else if (BC_switch[edge_index] == NeumannBC) { //if the edge is located on Neumann type boundary
						//******************************* The Neumann BC case: *********************************
						// find NeuMatrix = inverse(NeuMatrix_Orig)
						for (int i = 0; i < Knod; ++i) { // row index
							for (int j = 0; j < Knod; ++j) //column index of the matrix
								NeuMatrix_Orig(i,j) = -GB[el][d][t][i][1 - d] * sps_sps_grad_basis[j][i];
							NeuMatrix_Orig(i,i) += -GB[el][d][t][i][d] * g_prime[t];
						}

						NeuMatrix = NeuMatrix_Orig.inverse();

						double sign = t ? 1. : -1.;

						for (int j = 0; j < Knod; ++j) {
							for (int i = 0; i < Knod; ++i) {
								ij = j * Knod + i; //local cumulative sps index in element k
								
								int alpha = i * d + j * (1 - d);
								boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * alpha - 1) + alpha] += sign * sps_grad_radau[j * d + i * (1 - d)][t] 
									 * std::sqrt(GB[el][d][t][alpha][d]*face_jac[el][ijp][alpha]); //page 9 on my notes 17 Feb 2021, effect of NBC, according to the edge[edge_index] direction

								for (int m = 0; m < Knod; ++m) {
									for (int r = 0; r < 2; r++) {
										a = m * (1 - r) + i * r;
										b = m * r + j * (1 - r);
										q = a * d + b * (1 - d);
										double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
										double NGR1 = sps_grad_radau[a * (1 - d) + b * d][t];
										double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
										tmp2 = 0.;
										for (int s = 0; s < 2; s++) tmp2 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s];
										for (int p = 0; p < Knod; ++p) {
											j1 = p * d + b * (1 - d);
											i1 = p * (1 - d) + a * d;
											rs = j1 * Knod + i1; //cumulative index
											laplacian_center[el][ij][rs] -= 0.5 * SNGLB1 * A * sps_boundary_basis[p][t] * NGR1;  //last term of page 8 of my notes 17 Feb 2021
											laplacian_center[el][ij][rs] += 0.5 * tmp2 * A * sps_boundary_basis[p][t] * NGR1;  //last term of page 9 of my notes 17 Feb 2021

											boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * p - 1) + p] += A * NGR1 * sign * NeuMatrix(q, p) *
												std::sqrt(GB[el][d][t][p][d] * face_jac[el][ijp][p]) * (tmp2 - SNGLB1); //NBC terms on page 8,9 of mu notes 17 Feb 2021

											double tmp4 = sps_boundary_grad_basis[p][t] - g_prime[t] * sps_boundary_basis[p][t]; //SNGLB(p,t)-gLp(t)*SBLB(p,t)
											for (int c = 0; c < Knod; ++c) {
												j1 = p * d + c * (1 - d);
												i1 = p * (1 - d) + c * d;
												rs = j1 * Knod + i1; //cumulative index
												laplacian_center[el][ij][rs] += SNGLB1 * A * NGR1 * tmp4 * NeuMatrix(q,c) * GB[el][d][t][c][d]; //page 8 of my notes 17 Feb 2021
												laplacian_center[el][ij][rs] -= tmp2 * A * NGR1 * tmp4 * NeuMatrix(q, c) * GB[el][d][t][c][d];  //page 9 of my notes 17 Feb 2021
											}
										}
									}
								}
							}
						}
					} //if the edge is located on the Neumann - type global boundary


				} //if face is located on the global boundary

				else { // if this face is not on the global boundary: fill laplacian_neighbor[el][elem_side][][]
					nnz += K4;
					eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
					if (problem_type == 10) ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
					else ijpm = nbr[ijp];  //local ijp of the neighboring element (eln)
					dn = ijpm / 2;
					tn = ijpm % 2;

					for (int j = 0; j < Knod; ++j) {
						for (int i = 0; i < Knod; ++i) {
							int ij = j * Knod + i; //local cumulative sps index in element k
							for (int m = 0; m < Knod; ++m) {
								for (int r = 0; r < 2; r++) {
									a = m * (1 - r) + i * r;
									b = m * r + j * (1 - r);
									q = a * d + b * (1 - d);
									Mq = ((d + t + dn + tn + 1) % 2) * (Knod - 2 * q - 1) + q; //index of the flux point in the neighbor element
									double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
									double NGR1 = sps_grad_radau[a * (1 - d) + b * d][t];
									double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
									tmp2 = 0.;
									for (int s = 0; s < 2; s++) tmp2 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s];
									for (int p = 0; p < Knod; ++p) {
										j1 = p * dn + Mq * (1 - dn); //j of the neighbor element
										i1 = p * (1 - dn) + Mq * dn; //i of the neighbor element
										rs = j1 * Knod + i1; //cumulative index in the neighbor element
										laplacian_neighbor[el][elem_side][ij][rs] += 0.5 * SNGLB1 * A * sps_boundary_basis[p][tn] * NGR1;  //contribution to the L term, page 1 of my notes, 3 Feb 2021
										laplacian_neighbor[el][elem_side][ij][rs] -= 0.5 * tmp2   * A * sps_boundary_basis[p][tn] * NGR1;  //contribution to the M,N,P,Q term, page 2 of my notes, 3 Feb 2021
									}
								}
							}
						}
					}
					//page 3 of my notes, 3 Feb 2021
					for (int j = 0; j < Knod; ++j) {
						for (int i = 0; i < Knod; ++i) {
							ij = j * Knod + i; //local cumulative sps index in element k
							int alpha = i * d + j * (1 - d);
							int Malpha = ((d + t + dn + tn + 1) % 2) * (Knod - 2 * alpha - 1) + alpha;
							double A = GB[el][d][t][alpha][d];  //[(g_d . g_d)/J]|k,d,t,alpha
							double B = GB[el][d][t][alpha][1 - d]; //[(g_d . g_1-d)/J]|k,d,t,alpha
							double An = GB[eln][dn][tn][Malpha][dn]; //[(g_dn . g_dn)/J]|kn,dn,tn,Malpha
							double Bn = GB[eln][dn][tn][Malpha][1 - dn]; //[(g_dn . g_1-dn)/J]|kn,dn,tn,Malpha
							double NGR1 = sps_grad_radau[j * d + i * (1 - d)][t];
							double sign = (t == tn) ? -1. : 1.;

							for (int m = 0; m < Knod; ++m) {
								int Mm = ((d + t + dn + tn + 1) % 2) * (Knod - 2 * m - 1) + m;
								j1 = m * d + j * (1 - d);
								i1 = m * (1 - d) + i * d;
								rs = j1 * Knod + i1;
								laplacian_center[el][ij][rs] += (0.5 * A * sps_boundary_grad_basis[m][t] - 0.25 * g_prime[t] * A * sps_boundary_basis[m][t] +
									0.25 * sign * g_prime[tn] * An * sps_boundary_basis[m][t]) * NGR1;

								j1 = m * dn + Malpha * (1 - dn);
								i1 = m * (1 - dn) + Malpha * dn;
								rs = j1 * Knod + i1;
								laplacian_neighbor[el][elem_side][ij][rs] += (0.5 * sign * An * sps_boundary_grad_basis[m][tn] - 0.25 * sign * g_prime[tn] * An * sps_boundary_basis[m][tn] +
									0.25 * g_prime[t] * A * sps_boundary_basis[m][tn]) * NGR1;

								for (int p = 0; p < Knod; ++p) {
									j1 = p * d + m * (1 - d); 
									i1 = p * (1 - d) + m * d;
									rs = j1 * Knod + i1; //cumulative index in the neighbor element
									coeff = 0.25 * (B * sps_sps_grad_basis[m][alpha] + sign * Bn * sps_sps_grad_basis[Mm][Malpha]) * NGR1;
									laplacian_center[el][ij][rs] += coeff * sps_boundary_basis[p][t];
									j1 = p * dn + Mm * (1 - dn); //j of the neighbor element
									i1 = p * (1 - dn) + Mm * dn; //i of the neighbor element
									rs = j1 * Knod + i1; //cumulative index in the neighbor element
									laplacian_neighbor[el][elem_side][ij][rs] += coeff * sps_boundary_basis[p][tn];
								}
							}
						}
					}
				}  //else if the face is not on global boundary
			} //for t
		} //for d
	}  //for el

	
	if (LHS_type==1) Poisson_solver_Eigen_setup(laplacian_center, laplacian_neighbor);
	else if (LHS_type==2) Poisson_solver_Hypre_setup(laplacian_center, laplacian_neighbor);

	
	//  **************** free unnecessary memory ****************
	for (int i = 0; i < N_el; ++i) {
		for (int j = 0; j < Ksq; ++j)
			delete[] laplacian_center[i][j];
		delete[] laplacian_center[i];
	}
	delete[] laplacian_center;
	
	for (int i = 0; i < N_el; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < Ksq; ++k) 
				delete[] laplacian_neighbor[i][j][k];
			delete[] laplacian_neighbor[i][j];
		}
		delete[] laplacian_neighbor[i];
	}
	delete[] laplacian_neighbor;

}

char HO_2D::Euler_time_integrate() {
	//uses the Euler explicit time integration to solve the advection-diffusion equation for vorticity
	calc_RHS_advection(); //stores -div(vw) in RHS_advective[el][Knod][Knod]
	calc_RHS_diffusion();  //stores Laplace(w) in RHS_diffusive[el][Knod][Knod]
	
	//update the vorticity field now
	for (int el = 0; el < mesh.N_el; ++el)
		for (int ky = 0; ky < Knod; ++ky)
			for (int kx = 0; kx < Knod; ++kx)
				vorticity[el][ky][kx] += dt * (RHS_advective[el][ky][kx] + RHS_diffusive[el][ky][kx]);
	return 0;
}

char HO_2D::solve_Poisson() {
	//This subroutine solves the Poisson's equation for the stream function field
	int Ksq = Knod * Knod;

	int N_el = mesh.N_el;
	double* RHS = new double[N_el * Ksq];
	for (int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b)
		for (int m = 0; m < Knod; ++m) {
			//set the proper normal_vel[el_b][m] and BC_parl_vel[el_b][m] to construct proper Neuman BC for sai
			Cmpnts2 coor(0., 0.); //coordinate of sps[j][i]
			for (int n = 0; n < mesh.Lnod; ++n) {
				int node_index = mesh.edges[el_b].nodes[tensor2FEM(n)];
				coor.plus(gps_sps_basis[n][m], mesh.nodes[node_index].coor);
			}

			//BC_values[el_b][m] = coor.x + coor.y; //set the Dirichlet value for the sai values at the flux points on the boundary  
			BC_values[el_b][m] = 0.; //set the Dirichlet value for the sai values at the flux points on the boundary  
			double r = std::sqrt(coor.x * coor.x + coor.y * coor.y);
			if (r < 0.5101) BC_values[el_b][m] = -2. * r * std::exp(r * r); //cylinder surface
			//else if (mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y < -1.999 && mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y < -1.999) BC_values[el_b][m] = std::exp(coor.x*coor.x+4.); //set the Dirichlet value for the sai values at the flux points on the boundary  
			//else if (mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y > 1.999 && mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y > 1.999) BC_values[el_b][m] = std::exp(coor.x * coor.x + 4.); //set the Dirichlet value for the sai values at the flux points on the boundary  
			//else if (mesh.nodes[mesh.edges[el_b].nodes[0]].coor.x < -1.999 && mesh.nodes[mesh.edges[el_b].nodes[1]].coor.x < -1.999) BC_values[el_b][m] = std::exp(coor.y * coor.y + 4.);
			else BC_values[el_b][m] = std::exp(coor.y * coor.y + coor.x*coor.x);
			//BC_values[el_b][m] = std::exp(r * r);
			//BC_values[el_b][m] = coor.x * coor.x + coor.y * coor.y;
			//if (r < 0.501) BC_values[el_b][m] = -2. * r;
			//double DelY = mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y - mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y;
			//double DelX = mesh.nodes[mesh.edges[el_b].nodes[1]].coor.x - mesh.nodes[mesh.edges[el_b].nodes[0]].coor.x;
			//double nx = DelY / std::sqrt(DelY*DelY + DelX*DelX);
			//double ny = - DelX / std::sqrt(DelY * DelY + DelX * DelX);
			//if (r < 0.501) BC_values[el_b][m] = 2. * std::exp(r * r) * (coor.x * nx + coor.y * ny);//2.*(coor.x*nx + coor.y*ny);
		}
	
	for (int el = 0; el < N_el; ++el) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				int ij = j * Knod + i;
				RHS[el * Ksq + ij] = -vorticity[el][j][i] * vol_jac[el][j][i];
				//RHS[el * Ksq + ij] = 0.; //temporary to just test the code for poisson
			}
		}
	}

	for (int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				int ij = j * Knod + i;
				int el = mesh.boundary_elem_ID[el_b].element_index;
				for (int alpha = 0; alpha < Knod; ++alpha)
					RHS[el * Ksq + ij] -= boundary_source[el_b][ij][alpha] * BC_values[el_b][alpha];
			}
		}
	}


	if (LHS_type == 1) {  //Eigen then construct the RHS by Eigen too
		//RHS_Eigen = Eigen::VectorXd::Constant(N_el*Knod*Knod, 0.);
		//RHS_Eigen = Eigen::VectorXd::Zero(N_el*Knod*Knod);// // Initiate the RHS to zero
		for (int el = 0; el < N_el; ++el)
			for (int ij = 0; ij < Ksq; ++ij) RHS_Eigen(el * Ksq + ij) = RHS[el * Ksq + ij];

		/*
		Eigen::VectorXd poisson_sol = cg_Eigen.solve(RHS_Eigen); // Conjugate gradient method
		std::cout << "#iterations:     " << cg_Eigen.iterations() << std::endl;
		std::cout << "estimated error: " << cg_Eigen.error() << std::endl;
		*/

		
		Eigen::VectorXd poisson_sol = bicg_Eigen.solve(RHS_Eigen);  //biCGstab method
		std::cout << "#iterations:     " << bicg_Eigen.iterations() << std::endl;
		std::cout << "estimated error: " << bicg_Eigen.error() << std::endl;
		

		/*
		Eigen::VectorXd poisson_sol = LU_Eigen.solve(RHS_Eigen); //sparseLU method
		*/

		for (int el = 0; el < N_el; ++el)
			for (int j = 0; j < Knod; ++j)
				for (int i = 0; i < Knod; ++i) {
					int index = el * Ksq + j * Knod + i;
					stream_function[el][j][i] = poisson_sol(index);
				}



	}  //if LHS_type==1 (Eigen)





	delete[] RHS;

	return 0;
}

char HO_2D::calc_RHS_advection() {
	// ********* This subroutine calculates the - div(VW), V is Cartesian velocity vector, W is vorticity *********
	int ijp, ijpm, eln, neighbor_side, kk;
	double bndr_vel;
	double** local_vort = new double* [Knod]; //local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	for (int i = 0; i < Knod; ++i) local_vort[i] = new double[2];
	double** local_vel  = new double* [Knod]; //local array to hold the contravariant flux on a row of csi [0] (xsi dir), and row of eta [1] direction (eta dir)
	for (int i = 0; i < Knod; ++i) local_vel[i] = new double[2];

	double ***bndr_vort, ***bndr_flx, ***upwnd_flx; //quantities per element, per ijp face per sps index
	bndr_vort = new double** [mesh.N_el]; //interpolated vorticity from a row/column of K nodes on the ijp faces
	bndr_flx = new double** [mesh.N_el]; //interpolated flux of uw and vw in ijp faces, i.e. w*f^tilda on ijp=0,1 faces and w*g^tilda on ijp=2,3 faces
	upwnd_flx = new double** [mesh.N_el]; //flux values at el element, in idir direction, at row_col row or column at k sps, hence [el][idir][row_col][k]= f^tilda[k] w[k] if idir=0 AND g^tilda[k] w[k] if idir=1

	for (int i = 0; i < mesh.N_el; ++i) {
		bndr_vort[i] = new double* [4];
		bndr_flx[i] = new double* [4];
		upwnd_flx[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			bndr_vort[i][j] = new double [Knod];
			bndr_flx[i][j] = new double[Knod];
			upwnd_flx[i][j] = new double[Knod];
		}
	}

	double**** disc_flx = new double*** [mesh.N_el];
	for (int i = 0; i < mesh.N_el; ++i) {
		disc_flx[i] = new double** [2]; //0 for xsi and 1: eta directions, so 0 means f^tilda and 1 means g^tilda
		for (int j = 0; j < 2; j++) {
			disc_flx[i][j] = new double* [Knod];
			for (int k = 0; k < Knod; k++) disc_flx[i][j][k] = new double[Knod];
		}
	}

	for (int i = 0; i < mesh.N_el; ++i) 
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < Knod; k++) 
				for (int m = 0; m < Knod; m++) disc_flx[i][j][k][m] = 0.;

	double** bndr_disc_flx = new double* [4];  //4 side of a cells
	for (int i = 0; i < 4; ++i) bndr_disc_flx[i] = new double[Knod];
	// ****************************************************************


	//Extrapolate the unknown, Phi and the Flux to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) { //form flux (U^1w/U^2w on xsi/eta dirs, respectively) on the 4 boundaries (bndr_flx) and all interior sps of all elements (disc_flx)
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (ijp = 0; ijp < 4; ++ijp) bndr_vort[el][ijp][row_col] = bndr_flx[el][ijp][row_col] = upwnd_flx[el][ijp][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_vort[i][0] = vorticity[el][row_col][i];  // xsi direction
				local_vort[i][1] = vorticity[el][i][row_col];  //eta direction
				local_vel[i][0] =  velocity_cart[el][row_col][i].x * vol_Dy_Dxsi[el][row_col][i].y - velocity_cart[el][row_col][i].y * vol_Dx_Dxsi[el][row_col][i].y;  //f^tilda = contravariant xsi FLUX component along xsi direction, based on eq. 2.11 in fotis class notes
				local_vel[i][1] = -velocity_cart[el][i][row_col].x * vol_Dy_Dxsi[el][i][row_col].x + velocity_cart[el][i][row_col].y * vol_Dx_Dxsi[el][i][row_col].x;  //g^tilda = contravariant eta FLUX component along eta direction,
			}

			ijp = 0; //ijp=0 (ibnd=idir=0)
			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				//discontinuous flux values (f^tilda, g^tilda) on (row_col,1:Knod) solution points for xsi direction and on(row_col, 1:Knod) solution points for eta direction
				for (int i = 0; i < Knod; ++i) disc_flx[el][idir][row_col][i] = local_vel[i][idir] * local_vort[i][idir];

				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					for (int m = 0; m < Knod; m++) bndr_vort[el][ijp][row_col] += local_vort[m][idir] * sps_boundary_basis[m][ibnd];
					int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
					if (mesh.elem_neighbor[el].is_on_boundary[elem_side])
						bndr_flx[el][ijp][row_col] = bndr_vort[el][ijp][row_col] * BC_normal_vel[mesh.elem_neighbor[el].boundary_index[elem_side]][row_col];
					else {
						//value of u* w; ijp = 0: on left boundary at the row_col^th row of sps, ijp = 1 : on right boundary at the row_col^th row of sps, ijp = 2 : on bottom boundary at the row_col^th column of sps, ijp = 3, top at row_col^th column
						for (int m = 0; m < Knod; ++m) bndr_flx[el][ijp][row_col] += local_vel[m][idir] * sps_boundary_basis[m][ibnd]; //flux at the ibnd boundary in idir direction
						bndr_flx[el][ijp][row_col] *= bndr_vort[el][ijp][row_col];
					}
					ijp++;
				}
				
			}
		}
	}

	for (int el = 0; el < mesh.N_el; el++) {
		for (ijp = 0; ijp < 4; ++ijp) {
			int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
			if (mesh.elem_neighbor[el].is_on_boundary[elem_side])
				for (int j = 0; j < Knod; ++j)
					upwnd_flx[el][ijp][j] = bndr_flx[el][ijp][j];
			else {
				eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
				if (problem_type == 10) ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
				else ijpm = nbr[ijp];  //local ijp of the neighboring element (eln) : face to the left/south of left/south face or right/north of right/north face

				bool revert_tangent = (ijp == ijpm) || (ijp + ijpm == 3);
				bool revert_normal = (ijp == ijpm) || (std::fabs(ijp - ijpm) == 2);

				double sign = revert_normal? -1. : 1.;
				for (int j = 0; j < Knod; ++j) {
					kk = (revert_tangent) ? Knod - 1 - j : j;
											
					upwnd_flx[el][ijp][j] = 0.5 * (bndr_flx[el][ijp][j] + sign * bndr_flx[eln][ijpm][kk]);
					double vort_diff = bndr_vort[el][ijp][j] - bndr_vort[eln][ijpm][kk];
					if (std::fabs(vort_diff) > 1.e-6) {
						bndr_vel = (bndr_flx[el][ijp][j] - sign * bndr_flx[eln][ijpm][kk]) / vort_diff; //a_tilda
						upwnd_flx[el][ijp][j] += 0.5 * sgn[ijp] * std::fabs(bndr_vel) * vort_diff;
					}
				}
			}
		}
	}

	for (int el = 0; el < mesh.N_el; el++) {
		ijp = 0;
		for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
			for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
				for (int j = 0; j < Knod; ++j) {
					bndr_disc_flx[ijp][j] = 0.;
					for (int m = 0; m < Knod; ++m) bndr_disc_flx[ijp][j] -= disc_flx[el][idir][j][m] * sps_boundary_basis[m][ibnd]; // - f(1 or -1))
					bndr_disc_flx[ijp][j] += upwnd_flx[el][ijp][j];  //f_upw - f(1 or -1))
				}
				ijp++;
			}
		}

		//Now get the convection
		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i) {
				RHS_advective[el][j][i] = 0.;
				for (int dummy = 0; dummy < Knod; ++dummy) {
					RHS_advective[el][j][i] -= disc_flx[el][0][j][dummy] * sps_sps_grad_basis[dummy][i];  //xsi direction: -d(U^1w)/dxsi at (j,i)
					RHS_advective[el][j][i] -= disc_flx[el][1][i][dummy] * sps_sps_grad_basis[dummy][j];	//eta direction: -d(U^2w)/deta at (j,i)
				}
				RHS_advective[el][j][i] -= bndr_disc_flx[0][j] * sps_grad_radau[i][0] + bndr_disc_flx[1][j] * sps_grad_radau[i][1];
				RHS_advective[el][j][i] -= bndr_disc_flx[2][i] * sps_grad_radau[j][0] + bndr_disc_flx[3][i] * sps_grad_radau[j][1];
				RHS_advective[el][j][i] /= vol_jac[el][j][i]; //confirmed
			}
	}
	



	// ********************** free memory on the heap*****************
	for (int i = 0; i < 4; ++i) 
		delete[] bndr_disc_flx[i];
	delete[] bndr_disc_flx;

	for (int i = 0; i < mesh.N_el; ++i) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < Knod; k++)
				delete[] disc_flx[i][j][k];
			delete[] disc_flx[i][j];
		}
		delete[] disc_flx[i];
	}
	delete[] disc_flx;

	for (int i = 0; i < mesh.N_el; ++i) {
		for (int j = 0; j < 4; ++j) {
			delete[] bndr_vort[i][j];
			delete[] bndr_flx[i][j];
			delete[] upwnd_flx[i][j];
		}
		delete[] bndr_vort[i];
		delete[] bndr_flx[i];
		delete[] upwnd_flx[i];
	}
	delete[] bndr_vort;
	delete[] bndr_flx;
	delete[] upwnd_flx;

	for (int i = 0; i < Knod; ++i) {
		delete[] local_vort[i];
		delete[] local_vel[i];
	}
	delete[] local_vort;
	delete[] local_vel;
	return 0;
}

char HO_2D::calc_RHS_diffusion() {
	//this method calculates the 1/Re * Laplacian of vorticity at all sps of all elements stored in Lap_vorticity[el][ky][kx]
	
	//   update the BC based on below later
	for (int Gboundary = 0; Gboundary < mesh.N_Gboundary; ++Gboundary) { //usually 4 in case of square
		if (BC_no_slip[Gboundary]) {  //no slip condition on the global boundary element Gboundary, so if no slip wall then
			for (int edge_index : mesh.boundaries[Gboundary].edges) {  //loop over edges on the Gboundary global boundary
				for (int m = 0; m < Knod; ++m) velocity_jump[edge_index][m] -= BC_parl_vel[edge_index][m];
				for (int m = 0; m < Knod; ++m) BC_values[edge_index][m] = velocity_jump[edge_index][m] / (dt * Reyn_inv);
			}
		}
		else {
			for (int edge_index : mesh.boundaries[Gboundary].edges)  //loop over edges on the Gboundary global boundary
				for (int m = 0; m < Knod; ++m) BC_values[edge_index][m] = 0.; //clear this with Adrin when I get the chance (why the normal derivative is zero? for simple boundary the concavity is not zero)
		}
	}

	for (int edge_index = 0; edge_index < mesh.N_edges_boundary; ++edge_index) BC_switch[edge_index] = NeumannBC;

	// ****************** definition, memory allocation and initialization ****************
	double du_dxsi, du_deta;
	int eln;
	unsigned char ijp, ijpm; //left boundary(xsi): ijp=0, right_boundary (xsi): ijp=1; south boundary (eta): ijp=2; north boundary (eta): ijp = 3
	double** local_vort = new double* [Knod]; //local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	for (int i = 0; i < Knod; ++i) local_vort[i] = new double[2];

	double*** bndr_vort, *** bndr_grad_vort; //vorticity and grad of vorticity at sps of the L/R (0:1) and S/N (2:3) boundaries of all elements
	double ***comm_vort, *** comm_grad_vort; //common vorticity and common grad vorticity per element, per face (west, east, south, north) per sps index
	bndr_vort = new double** [mesh.N_el];
	bndr_grad_vort = new double** [mesh.N_el];
	comm_vort = new double** [mesh.N_el];
	comm_grad_vort = new double** [mesh.N_el];

	for (int i = 0; i < mesh.N_el; ++i) {
		bndr_vort[i] = new double* [4];
		bndr_grad_vort[i] = new double* [4];
		comm_vort[i] = new double* [4];
		comm_grad_vort[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			bndr_vort[i][j] = new double[Knod];
			bndr_grad_vort[i][j] = new double[Knod];
			comm_vort[i][j] = new double[Knod];
			comm_grad_vort[i][j] = new double[Knod];
		}
	}

	double** f_tilda, ** g_tilda;
	f_tilda = new double* [Knod]; //f_tilda is the derivative of the continuous vorticity function at all sps, i.e. d/dxsi(w^C)
	g_tilda = new double* [Knod];  //g_tilda is the derivative of the continuous vorticity function at all sps, i.e. d/deta(w^C)
	for (int i = 0; i < Knod; ++i) {
		f_tilda[i] = new double[Knod];
		g_tilda[i] = new double[Knod];
	}

	double **f_tilda_B = new double* [2], **g_tilda_B = new double* [2];
	for (int i = 0; i < 2; ++i) {
		f_tilda_B[i] = new double[Knod];
		g_tilda_B[i] = new double[Knod];
	}


	//**********************************************************************

	//Extrapolate the unknown vorticity and its derivative to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) { //form vorticity xsi/eta dirs, respectively) on the 4 boundaries (bndr_flx) and all interior sps of all elements (disc_flx)
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (int side = 0; side < 4; ++side) bndr_vort[el][side][row_col] = bndr_grad_vort[el][side][row_col] = comm_vort[el][side][row_col] = comm_grad_vort[el][side][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_vort[i][0] = vorticity[el][row_col][i];  // xsi direction, sps at row i
				local_vort[i][1] = vorticity[el][i][row_col];  //eta direction, sps at column i
			}

			ijp = 0;

			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					for (int m = 0; m < Knod; m++) bndr_vort[el][ijp][row_col] += local_vort[m][idir] * sps_boundary_basis[m][ibnd]; //vorticity on the ijp side on row/col row_col
					for (int m = 0; m < Knod; ++m) bndr_grad_vort[el][ijp][row_col] += local_vort[m][idir] * sps_boundary_grad_basis[m][ibnd]; //derivative of vorticity in xsi/eta dirs on the ijp side on row / col row_col: west/east:xsi derivative; north/south: eta derivative
					ijp++;
				}
			}
		}
	}

	if (HuynhSolver_type == 2) { //Use g = g_DG, so presumably DG method
		for (int el = 0; el < mesh.N_el; el++) {
			ijp = 0;
			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
					if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) {
						int el_B = mesh.elem_neighbor[el].boundary_index[elem_side];
						calc_boundary_comm_vals_meth2(el, el_B, ijp, ibnd, bndr_vort[el][ijp], bndr_grad_vort[el][ijp], comm_vort[el][ijp], comm_grad_vort[el][ijp]);
						for (int k = 0; k < Knod; ++k) velocity_jump[el_B][k] = comm_vort[el][ijp][k];
					}
					else { //if the element boundary ijp is internal
						eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
						if (problem_type==10) ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
						else ijpm = nbr[ijp];  //local ijp of the neighboring element (eln)
						
						calc_internal_comm_vals_meth2(el, ijp, ijpm, ibnd, bndr_vort[el][ijp], bndr_vort[eln][ijpm], bndr_grad_vort[el][ijp], comm_vort[el][ijp]);
					}
					ijp++;
				}
			}
		}

		for (int el = 0; el < mesh.N_el; el++) {
			for (ijp = 0; ijp < 4; ++ijp) {
				cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to ijp
				if (!(mesh.elem_neighbor[el].is_on_boundary[elem_side])) { //if it is internal local boundary
					eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
					if (problem_type==10) ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
					else ijpm = nbr[ijp];
					//comm_grad_vort[el][ijp][Knod] is Average of bndr_grad_vort[el][ijp][Knod] + [eln][ijpm][Knod]; same works for Left, South, and North
					bool revert_tangent = (ijp == ijpm) || (ijp + ijpm == 3);
					bool revert_normal = (ijp == ijpm) || (std::fabs(ijp - ijpm) == 2);
					if (revert_tangent)
						if (revert_normal)
							for (int k = 0; k < Knod; ++k) comm_grad_vort[el][ijp][k] = 0.5 * (bndr_grad_vort[el][ijp][k] - bndr_grad_vort[eln][ijpm][Knod - 1 - k]);
						else 
							for (int k = 0; k < Knod; ++k) comm_grad_vort[el][ijp][k] = 0.5 * (bndr_grad_vort[el][ijp][k] + bndr_grad_vort[eln][ijpm][Knod - 1 - k]);
					else
						if (revert_normal)
							for (int k = 0; k < Knod; ++k) comm_grad_vort[el][ijp][k] = 0.5 * (bndr_grad_vort[el][ijp][k] - bndr_grad_vort[eln][ijpm][k]);
						else
							for (int k = 0; k < Knod; ++k) comm_grad_vort[el][ijp][k] = 0.5 * (bndr_grad_vort[el][ijp][k] + bndr_grad_vort[eln][ijpm][k]);
				}
			}
		}
	}

	double dw_dxsi, dw_deta;
	for (int el = 0; el < mesh.N_el; ++el) {
		//calculate first derivative of the variable
		for (int ky = 0; ky < Knod; ++ky) {
			for (int kx = 0; kx < Knod; ++kx) {
				//calc grads of unknown uin along xsi(for eta = const) and along eta(for xsi = const); du_dxsi and du_deta are the du / dxsi and du / deta at sps[ky][kx]
				dw_dxsi = dw_deta = 0.;
				for (int dummy = 0; dummy < Knod; ++dummy) dw_dxsi += vorticity[el][ky][dummy] * sps_sps_grad_basis[dummy][kx];
				dw_dxsi += (comm_vort[el][0][ky] - bndr_vort[el][0][ky]) * sps_grad_radau[kx][0];
				dw_dxsi += (comm_vort[el][1][ky] - bndr_vort[el][1][ky]) * sps_grad_radau[kx][1];

				for (int dummy = 0; dummy < Knod; ++dummy) dw_deta += vorticity[el][dummy][kx] * sps_sps_grad_basis[dummy][ky];
				dw_deta += (comm_vort[el][2][kx] - bndr_vort[el][2][kx]) * sps_grad_radau[ky][0];
				dw_deta += (comm_vort[el][3][kx] - bndr_vort[el][3][kx]) * sps_grad_radau[ky][1];

				// now calculate f^tilda g^tilda as per Huynh's paper, i.e. eq 7.16 and a line below in Hyun diff paper
				f_tilda[ky][kx] = (dw_dxsi * (vol_Dx_Dxsi[el][ky][kx].y * vol_Dx_Dxsi[el][ky][kx].y + vol_Dy_Dxsi[el][ky][kx].y * vol_Dy_Dxsi[el][ky][kx].y)
					- dw_deta * (vol_Dx_Dxsi[el][ky][kx].x * vol_Dx_Dxsi[el][ky][kx].y + vol_Dy_Dxsi[el][ky][kx].x * vol_Dy_Dxsi[el][ky][kx].y)) / vol_jac[el][ky][kx];

				g_tilda[ky][kx] = (dw_deta * (vol_Dx_Dxsi[el][ky][kx].x * vol_Dx_Dxsi[el][ky][kx].x + vol_Dy_Dxsi[el][ky][kx].x * vol_Dy_Dxsi[el][ky][kx].x)
					- dw_dxsi * (vol_Dx_Dxsi[el][ky][kx].x * vol_Dx_Dxsi[el][ky][kx].y + vol_Dy_Dxsi[el][ky][kx].x * vol_Dy_Dxsi[el][ky][kx].y)) / vol_jac[el][ky][kx];
			}
		}
		// calculate tailda values at the mesh boundaries
		for (int ibnd = 0; ibnd < 2; ++ibnd) {
			for (int k = 0; k < Knod; k++) {
				f_tilda_B[ibnd][k] = comm_grad_vort[el][ibnd][k];
				for (int dummy = 0; dummy < Knod; ++dummy) f_tilda_B[ibnd][k] -= f_tilda[k][dummy] * sps_boundary_basis[dummy][ibnd]; //common derivative - extrapolated derivative of continuous function at left/right of mesh
				g_tilda_B[ibnd][k] = comm_grad_vort[el][ibnd + 2][k];
				for (int dummy = 0; dummy < Knod; ++dummy) g_tilda_B[ibnd][k] -= g_tilda[dummy][k] * sps_boundary_basis[dummy][ibnd]; //common derivative - extrapolated derivative of continuous function at left/right of mesh
			}
		}
	
		//now calculate the Laplacian
		for (int ky = 0; ky < Knod; ++ky) {
			for (int kx = 0; kx < Knod; ++kx) {
				RHS_diffusive[el][ky][kx] = 0.;
				for (int dummy = 0; dummy < Knod; ++dummy) RHS_diffusive[el][ky][kx] += f_tilda[ky][dummy] * sps_sps_grad_basis[dummy][kx];
				RHS_diffusive[el][ky][kx] += f_tilda_B[0][ky] * sps_grad_radau[kx][0] + f_tilda_B[1][ky] * sps_grad_radau[kx][1]; //contribution from xsi dir

				for (int dummy = 0; dummy < Knod; ++dummy) RHS_diffusive[el][ky][kx] += g_tilda[dummy][kx] * sps_sps_grad_basis[dummy][ky];
				RHS_diffusive[el][ky][kx] += g_tilda_B[0][kx] * sps_grad_radau[ky][0] + g_tilda_B[1][kx] * sps_grad_radau[ky][1];
				RHS_diffusive[el][ky][kx] *= Reyn_inv/vol_jac[el][ky][kx];

			}
		}
	}
	
	// ******************* free memory on heap ***************
	for (int i = 0; i < Knod; ++i)
		delete[] local_vort[i];
	delete[] local_vort;

	for (int j = 0; j < mesh.N_el; ++j) {
		for (int i = 0; i < 4; ++i) {
			delete[] bndr_vort[j][i];
			delete[] bndr_grad_vort[j][i];
			delete[] comm_vort[j][i];
			delete[] comm_grad_vort[j][i];
		}
		delete[] bndr_vort[j];
		delete[] bndr_grad_vort[j];
		delete[] comm_vort[j];
		delete[] comm_grad_vort[j];
	}
	delete[] bndr_vort;
	delete[] bndr_grad_vort;
	delete[] comm_vort;
	delete[] comm_grad_vort;

	for (int i = 0; i < Knod; ++i) {
		delete[] f_tilda[i];
		delete[] g_tilda[i];
	}
	delete[] f_tilda;
	delete[] g_tilda;

	for (int i = 0; i < 2; ++i) {
		delete[] f_tilda_B[i];
		delete[] g_tilda_B[i];
	}
	delete[] f_tilda_B;
	delete[] g_tilda_B;

	return 0;
}

void HO_2D::calc_internal_comm_vals_meth2(int el, int ijp, int ijpm, int ibnd, double* vort, double* vort_nbr, double* Dvort, double* com_vort) {
	/*
	vort[0:Knod-1]: in: Knod vorticity values on the ijp internal boundary of a cell which stores extrapolated vorticity from the internal sps nodes on xsi or eta constant sps
	vort_nbr[0:Knod-1]: in: Knod vorticity values on the same local boundary as in vort, but extrapolated from the vorticity of the neighboring sps. so vort and vort_nbr correspond to the same flux points on the boundary, but extrapolated from two separate sides
	Dvort[0:Knod-1]: in: gradient of Knod vorticity at the flux points of ijp side of element el. These gradients are calculated based extrapolation of sps vorticity*derivative of shape functions; out: normal flux to the global boundary on the ijp side
	com_vort[0:Knod-1]: out: The common vorticity (average of left and right (top+bottom) at all flux points on the ijp side of element el
	*/

	double* cross_Dvort = new double[Knod];   //cross derivatives using common values com_vort
	// if the orientation on the two sides of the interface are different, revert the interface bndr_vort of the neighboring element
	bool revert = (ijp == ijpm) || (ijp + ijpm == 3);
	if (revert) 
		for (int k = 0; k < Knod; ++k) com_vort[k] = 0.50 * (vort[k] + vort_nbr[Knod - k-1]); //The common value of Phi at all flux points on one edge of element(eq . 3.1 in hyun diffusion paper))
	else 
		for (int k = 0; k < Knod; ++k) com_vort[k] = 0.50 * (vort[k] + vort_nbr[k]); //The common value of Phi at all flux points on one edge of element(eq . 3.1 in hyun diffusion paper))

	//calc the corrected values of grad(Phi) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG' = -gRB' (=-Knod^2/2)

	for (int k = 0; k < Knod; ++k) Dvort[k] += (com_vort[k] - vort[k]) * g_prime[ibnd]; //derivative of eq 3.7 to use in 4.3 (common derivative))

	//calc the values of grad(Phi) along the mesh interface;

	for (int k = 0; k < Knod; ++k) {
		cross_Dvort[k] = 0.;
		for (int dummy = 0; dummy < Knod; ++dummy)
			cross_Dvort[k] += com_vort[dummy] * sps_sps_grad_basis[dummy][k]; //derivative of common value(at the element edge) in the other direction, i.e.derivative w.r.t eta if ijp=0,1 and derivative w.r.t xsi if ijp=2,3
	}

	//calculate the common flux
	for (int k = 0; k < Knod; ++k)
		Dvort[k] = (Dvort[k] * face_Acoef[el][ijp][k] - cross_Dvort[k] * face_Bcoef[el][ijp][k]) / face_jac[el][ijp][k]; //the common flux normal to the global boundary

	delete[] cross_Dvort;
}

void HO_2D::calc_boundary_comm_vals_meth2(int el, int el_b, int ijp, int ibnd, double* vort, double* Dvort, double* com_vort, double* com_grad_vort) {
	/*
	el_b: global boundary index
	vort = bndr_vort[el][ijp][0:Knod-1]: extrapolated vorticity at the ijp side of an element
	Dvort = bndr_grad_vort[el][ijp][0:Knod-1]; in: the extrapolated grad of vorticity on the boundary, out: the corrected grad of vorticity on the boundary
	com_vort = comm_vort[el][ijp][0:Knod-1; out: common vorticity value at the ijp side (which is on global boundary)
	com_grad_vort = comm_grad_vort[el][ijp]; out: the common flux on the global boundary
	*/

	double* cross_Dvort = new double[Knod];   //cross derivatives using common values com_vort
	double** NeuMatrix = new double* [Knod], * NeuRHS = new double[Knod];
	for (int m = 0; m < Knod; ++m) NeuMatrix[m] = new double[Knod];

	if (BC_switch[el_b] == DirichletBC) {  //the vorticity (velocity) is provided on the global boundary as BC

		for (int m = 0; m < Knod; ++m) com_vort[m] = BC_values[el_b][m]; //the vorticity at the global boundary, this is taken as the common value

		//Get the corrected values of grad(vort) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG' = -gRB' (=-Knod^2/2)
		for (int k = 0; k < Knod; ++k)
			Dvort[k] += (com_vort[k] - vort[k]) * g_prime[ibnd]; //corrected grad vorticity at the ibnd global boundary

		//calc the corrected values of grad(Phi) along the mesh interface (parallel to the global boundary);
		for (int k = 0; k < Knod; ++k) {
			cross_Dvort[k] = 0.;
			for (int m = 0; m < Knod; ++m)
				cross_Dvort[k] += com_vort[m] * sps_sps_grad_basis[m][k];  // parallel to the boundary component of grad(vorticity)
		}

		//calculate the common flux
		for (int k = 0; k < Knod; ++k)
			com_grad_vort[k] = (Dvort[k] * face_Acoef[el][ijp][k] - cross_Dvort[k] * face_Bcoef[el][ijp][k]) / face_jac[el][ijp][k]; // the common flux normal to the global boundary

	}
	else if (BC_switch[el_b] == NeumannBC) { //the normal gradient of vorticity is provided on the global boundary as BC: grad(w).n^hat
		for (int k = 0; k < Knod; ++k)
			com_grad_vort[k] = BC_values[el_b][k] * face_Anorm[el][ijp][k];  //flux created from the vorticity normal gradient: serves as common flux

		//convert comon flux into common value for vorticity on the boundary el_b
		for (int k = 0; k < Knod; ++k) NeuRHS[k] = face_jac[el][ijp][k] * com_grad_vort[k] - face_Acoef[el][ijp][k] * (Dvort[k] - vort[k] * g_prime[ibnd]); //right hand side vector to solve for common values of vorticity on the boundary
		
		for (int k = 0; k < Knod; ++k) {
			for (int i = 0; i < Knod; ++i)
				NeuMatrix[k][i] = -face_Bcoef[el][ijp][k] * sps_sps_grad_basis[i][k];
			NeuMatrix[k][k] += face_Acoef[el][ijp][k] * g_prime[ibnd];
			
		}			

		Gauss_solver(Knod, NeuMatrix, NeuRHS, com_vort);  //solves the system [NeuMatrix]*{com_vort} = {NeuRHS} to find com_vort
	}


	delete[] cross_Dvort;
	delete[] NeuRHS;
	for (int i = 0; i < Knod; ++i) delete[] NeuMatrix[i];
	delete[] NeuMatrix;
}

void HO_2D::Poisson_solver_Hypre_setup(double*** laplacian_center, double**** laplacian_neighbor) {
	// sets up the LHS of the Poisson equation via the HYPRE library
	//Create_Hypre_Matrix(laplacian_center, laplacian_neighbor);
	//Create_Hypre_Vector();

}

void HO_2D::Create_Hypre_Matrix () {
	//int localsize_p;
	//VecGetLocalSize(user->Phi2, &localsize_p);

	//int p_lower = user->p_global_begin;
	//int p_upper = p_lower + localsize_p - 1;

	//std::vector<int> nz_p(localsize_p);
	//std::fill(nz_p.begin(), nz_p.end(), 19);

	//PetscPrintf(TE_comm, "\nbegin HYPRE_IJMatrixCreate gap(%d)\n", user->gap_index);
	//HYPRE_IJMatrixCreate(TE_comm, p_lower, p_upper, p_lower, p_upper, &Ap);
	//HYPRE_IJMatrixSetObjectType(Ap, HYPRE_PARCSR);
	//HYPRE_IJMatrixSetRowSizes(Ap, &nz_p[0]);
	////HYPRE_IJMatrixSetMaxOffProcElmts (Ap, 10);
	//HYPRE_IJMatrixInitialize(Ap);
	//PetscPrintf(TE_comm, "end HYPRE_IJMatrixCreate for gap (%d)\n\n", user->gap_index);
}

void HO_2D::Poisson_solver_Eigen_setup(double*** laplacian_center, double**** laplacian_neighbor) {
	// sets up the LHS of the Poisson equation via the Eigen library
	int N_el = mesh.N_el;
	int Ksq = Knod * Knod, K4 = Ksq * Ksq;;
	std::vector<Trplet> coefficientsVec;            // list of non-zeros coefficients in the LHS matrix
	coefficientsVec.reserve(nnz);  //maximum number of triplets
	nnz = 0;
	for (int el = 0; el < N_el; ++el)
		for (int ij = 0; ij < Ksq; ++ij)
			for (int rs = 0; rs < Ksq; ++rs)
				if (std::fabs(laplacian_center[el][ij][rs]) > 1.e-14) {
					coefficientsVec.push_back(Trplet(el * Ksq + ij, el * Ksq + rs, laplacian_center[el][ij][rs]));
					nnz++;
				}

	for (int el = 0; el < N_el; ++el)
		for (int elem_side = south; elem_side <= west; elem_side++)
			if (!(mesh.elem_neighbor[el].is_on_boundary[elem_side])) {
				int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
				for (int ij = 0; ij < Ksq; ++ij)
					for (int rs = 0; rs < Ksq; ++rs)
						if (std::fabs(laplacian_neighbor[el][elem_side][ij][rs]) > 1.e-14) {
							coefficientsVec.push_back(Trplet(el * Ksq + ij, eln * Ksq + rs, laplacian_neighbor[el][elem_side][ij][rs]));
							nnz++;
						}




			}

	LHS_Eigen.resize(N_el * Ksq, N_el * Ksq);
	LHS_Eigen.data().squeeze();
	LHS_Eigen.setFromTriplets(coefficientsVec.begin(), coefficientsVec.end());
	LHS_Eigen.makeCompressed();

	// ************************ The Conjugate gradient method ******************
	/*
	cg_Eigen.analyzePattern(LHS_Eigen); 
	cg_Eigen.factorize(LHS_Eigen);
	*/
	 
	//****************** the BICG method ****************
	
	bicg_Eigen.compute(LHS_Eigen);
	
	// ****************** Sparse LU **********************
	/*
	LU_Eigen.analyzePattern(LHS_Eigen);
	LU_Eigen.factorize(LHS_Eigen);
	*/
	//***************************************************

	RHS_Eigen.resize(N_el * Ksq);
}

void HO_2D::save_output(int n) {
	std::string file_name = "problem";
	file_name.append(std::to_string(problem_type));
	file_name.append("_timestep");
	file_name.append(std::to_string(n));
	file_name.append(".dat");

	std::string file_name_norm = "Norms_timesteps_";
	file_name_norm.append(std::to_string(n));
	file_name_norm.append(".dat");
	std::cout << "     Writing the results after " << n << "  timesteps into the file:  " << file_name << std::endl;


	double time = n * dt;
	std::ofstream file_handle(file_name);
	std::ofstream file_handle_norm(file_name_norm);
	double u_exact = 0.;
	double Linf_Norm = -1.0, L1_Norm = 0.0, L2_Norm = 0.0;
	for (int el = 0; el < mesh.N_el; ++el) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				Cmpnts2 coor(0., 0.); //coordinate of sps[j][i]
				for (int ny = 0; ny < mesh.Lnod; ++ny)
					for (int nx = 0; nx < mesh.Lnod; ++nx) {
						int node_index = mesh.elements[el].nodes[tensor2FEM(nx, ny)];
						coor.plus(gps_sps_basis[ny][j] * gps_sps_basis[nx][i], mesh.nodes[node_index].coor);
					}

				if (problem_type == 3)
					u_exact = 0.25 * (3. * std::sin(M_PI * coor.x) * std::exp(-M_PI * M_PI * time) - std::sin(3. * M_PI * coor.x) * std::exp(-9.0 * M_PI * M_PI * time))
					* (3. * std::sin(M_PI * coor.y) * std::exp(-M_PI * M_PI * time) - std::sin(3. * M_PI * coor.y) * std::exp(-9.0 * M_PI * M_PI * time));

				else if (problem_type == 2) {
					u_exact = 0.0;
					double error = 1.0;
					for (int iter = 1; iter <= 200; ++iter) {
						double tmp = 0.50 * (2. * iter - 1.) * M_PI;
						u_exact += 2.0 * minus_one_to_power(iter + 1) * std::cos(tmp * coor.x) * exp(-tmp * tmp * time) / tmp;
						if (std::fabs(error) < 1.E-10) break;
					}
				}

				else if (problem_type == 8) {//heat conduction on cylinder (periodic BC)
					double r = coor.norm2();
					double beta[] = { 1., 2.,3. }; //temporary fixation

					double tmp = 0.;
					for (int iter = 0; iter < 45; ++iter) {
						tmp = tmp + exp(-beta[iter] * beta[iter] * time) * _j0(0.5 * beta[iter]) *
							(_j0(beta[iter] * r) * _y0(beta[iter]) - _j0(beta[iter]) * _y0(beta[iter] * r)) /
							(_j0(0.5 * beta[iter]) + _j0(beta[iter]));
					}
					u_exact = M_PI * tmp;
				}

				else if (problem_type == 10) {
					//u_exact = coor.x + coor.y;
					double ym = 1. - coor.y;
					u_exact = ym * (1. + std::cos(M_PI * ym)) * sin(M_PI * coor.x);
					u_exact = std::sin(M_PI * coor.y) / (M_PI*M_PI);
					u_exact = 2. * coor.x * coor.x - coor.y * coor.y - coor.x + coor.y + coor.x * coor.y;
					u_exact = 3.*coor.x + coor.y-1;
					u_exact = 3. * coor.x * coor.x * coor.y - 2. * coor.x * coor.x * coor.x + coor.x * coor.y * coor.y + coor.y * coor.y * coor.y + 1.0;
					u_exact = coor.y * (1.+std::cos(M_PI*coor.y)) * std::sin(M_PI*coor.x);
					u_exact = std::exp(coor.x * coor.x + coor.y * coor.y);
					//u_exact = coor.x * coor.x + coor.y * coor.y;
				}

				Linf_Norm = std::max(Linf_Norm, fabs(u_exact - /*vorticity*/stream_function[el][j][i]));
				L1_Norm = L1_Norm + fabs(u_exact - /*vorticity*/stream_function[el][j][i]);
				L2_Norm = L2_Norm + std::pow(u_exact - /*vorticity*/stream_function[el][j][i], 2);

				file_handle << el << " \t " << j << " \t " << i << " \t " << coor.x << " \t " << coor.y << " \t " << /*vorticity*/stream_function[el][j][i] << " \t " << u_exact << " \t " << 100. * (u_exact - /*vorticity*/stream_function[el][j][i]) / u_exact << std::endl;
			}
		}
	}
	L1_Norm /= ((double)mesh.N_el * Knod * Knod);
	L2_Norm = std::sqrt(L2_Norm / ((double)mesh.N_el * Knod * Knod));
	file_handle_norm << L1_Norm << " \t " << L2_Norm << " \t " << Linf_Norm;
	std::cout << "Done writing to file" << std::endl;
	file_handle.close();
	file_handle_norm.close();
}