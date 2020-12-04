#include "calculation.h"
#include "misc.hpp"

void HO_2D::release_memory() { //release the memory as destructor
	//delete[] vert_coor; delete[] vert_coor_initial; delete[] vert_coor_m1; delete[] Hex_verts_comps; delete[] projected_hex_verts_comps;
	//delete[] zone_id; delete[] group_regions; delete[] facet_groups; delete[] cell_N_verts; delete[] face_N_verts; delete[] N_connectivity;
	for (int k = 0; k < 10/*mesh.N_elements*/; ++k) {
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
		sps_local_coor[1] = 1. / sqrt(3.);    //collocation nodes per elem are{ -C1, C1 }
		sps_local_coor[0] = -sps_local_coor[1];
		sps_weight[0] = sps_weight[1] = 1.0;  //Gaussian weights
	}
	else if (Knod == 3) {
		sps_local_coor[2] = sqrt(0.6);	//collocation nodes per elem are{ -C1, 0, C1 }
		sps_local_coor[1] = 0.;
		sps_local_coor[0] = -sps_local_coor[2];
		sps_weight[2] = sps_weight[0] = 5.0 / 9.0;
		sps_weight[1] = 8.0 / 9.0;
	}
	else if (Knod == 4) {
		sps_local_coor[3] = sqrt((15. + 2. * sqrt(30.)) / 35.);	//collocation nodes per elem are{ -C2, -C1, C1, C2 }
		sps_local_coor[2] = sqrt((15. - 2. * sqrt(30.)) / 35.);
		sps_local_coor[1] = -sps_local_coor[2];
		sps_local_coor[0] = -sps_local_coor[3];
		sps_weight[3] = sps_weight[0] = (18. - sqrt(30.)) / 36.;
		sps_weight[2] = sps_weight[1] = (18. + sqrt(30.)) / 36.;
	}
	else if (Knod == 5) {
		sps_local_coor[4] = sqrt(5. + 2. * sqrt(70.) / 7.) / 3.;     //collocation nodes per elem are{ -C2, -C1, 0, C1, C2 }
		sps_local_coor[3] = sqrt(5. - 2. * sqrt(70.) / 7.) / 3.;
		sps_local_coor[2] = 0.;
		sps_local_coor[1] = -sps_local_coor[3];
		sps_local_coor[0] = -sps_local_coor[4];
		sps_weight[0] = sps_weight[4] = (322. - 13. * sqrt(70.)) / 900.;
		sps_weight[1] = sps_weight[3] = (322. + 13. * sqrt(70.)) / 900.;
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
		g_prime[1] = 0.5 * minus_one_to_power(Knod) * Knod;
	}

}

char HO_2D::allocate_arrays() {
	//allocate the memory for the arrays and resize them.

	initial_vorticity = new double** [mesh.N_el];
	vorticity = new double** [mesh.N_el];
	velocity_cart = new Cmpnts2** [mesh.N_el];
	for (int j = 0; j < mesh.N_el; ++j) {
		initial_vorticity[j] = new double* [Knod];
		vorticity[j] = new double* [Knod];
		velocity_cart[j] = new Cmpnts2* [Knod];
		for (int i = 0; i < Knod; ++i) {
			initial_vorticity[j][i] = new double[Knod];
			vorticity[j][i] = new double[Knod];
			velocity_cart[j][i] = new Cmpnts2[Knod];
		}
	}

	int N_el = mesh.N_el;
	sps_local_coor = new double[Knod];
	sps_weight = new double[Knod];
	gps_local_coor = new double[mesh.Lnod];
	sps_boundary_basis = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_boundary_basis[i] = new double[2]; //for left and right
	sps_boundary_grad_basis = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_boundary_grad_basis[i] = new double[2]; //for left and right
	sps_sps_grad_basis = new double* [Knod];
	for (int i = 0; i < Knod; ++i) sps_sps_grad_basis[i] = new double[Knod];

	gps_boundary_basis = new double* [mesh.Lnod]; for (int i = 0; i < mesh.Lnod; ++i) gps_boundary_basis[i] = new double[2];
	gps_boundary_grad_basis = new double* [mesh.Lnod]; for (int i = 0; i < mesh.Lnod; ++i) gps_boundary_grad_basis[i] = new double[2];
	gps_sps_basis = new double* [mesh.Lnod];
	gps_sps_grad_basis = new double* [mesh.Lnod];
	for (int i = 0; i < mesh.Lnod; ++i) gps_sps_basis[i] = new double[Knod];
	for (int i = 0; i < mesh.Lnod; ++i) gps_sps_grad_basis[i] = new double[Knod];

	sps_radau = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_radau[i] = new double[2]; //left[0] and right [1] radau functions on sps
	sps_grad_radau = new double* [Knod]; for (int i = 0; i < Knod; ++i) sps_grad_radau[i] = new double[2];

	vol_Dx_Dxsi = new Cmpnts2** [N_el];
	vol_Dy_Dxsi = new Cmpnts2** [N_el];
	vol_jac = new double** [N_el];
	for (int i = 0; i < N_el; ++i) {
		vol_Dx_Dxsi[i] = new Cmpnts2* [Knod];
		vol_Dy_Dxsi[i] = new Cmpnts2* [Knod];
		vol_jac[i] = new double* [Knod];
		for (int j = 0; j < Knod; ++j) {
			vol_Dx_Dxsi[i][j] = new Cmpnts2 [Knod];
			vol_Dy_Dxsi[i][j] = new Cmpnts2 [Knod];
			vol_jac[i][j] = new double[Knod];
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
	mesh.boundaries.resize(mesh.N_Gboundary);

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

			for (int i = 0; i < Knod; ++i) {
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
	BC_switch = new unsigned char [mesh.N_edges_boundary];
	
	BC_values = new double* [mesh.N_edges_boundary];
	velocity_jump = new double* [mesh.N_edges_boundary];
	BC_parl_vel = new double* [mesh.N_edges_boundary];
	BC_normal_vel = new double* [mesh.N_edges_boundary];
	for (int i = 0; i < mesh.N_edges_boundary; ++i) {
		BC_values[i] = new double[Knod];
		velocity_jump[i] = new double[Knod];
		BC_parl_vel[i] = new double[Knod];
		BC_normal_vel[i] = new double[Knod];
	}
	
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
			BC_normal_vel[el_b][m] = 0.;
		}
			

	if (problem_type == 1) {
		// Periodic moving wall problem; bottom velocity is 1; top is 0
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) BC_switch[el_b] = DirichletBC;
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) {
			if (fabs(mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y) < 1.e-6 && fabs(mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y) < 1.e-6) //boundary on the bottom
				for (int m = 0; m < Knod; ++m) BC_parl_vel[el_b][m] = 1.0;    //Wall parallel velocity
		}
	}

	else if (problem_type == 3) {
		//Square Cavity Problem; top velocity is 1
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) BC_switch[el_b] = DirichletBC;
		for (int el_b=0; el_b < mesh.N_edges_boundary; el_b++) {
			if (fabs(mesh.nodes[mesh.edges[el_b].nodes[0]].coor.y - 1.) < 1.e-6 && fabs(mesh.nodes[mesh.edges[el_b].nodes[1]].coor.y - 1.) < 1.e-6) //top boundary
				for (int m = 0; m < Knod; ++m) BC_parl_vel[el_b][m] = 1.0; //top wall parallel velocity
		}
	}



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
			!Haji : The Face_Acoef and Face_Bcoef, forms the covariant metric tensor g_ij = [g11 g12; g12, g22]] (symmetric matrix))
	*/
	unsigned int Lnod = mesh.Lnod;
	double dx_dxsi, dy_dxsi, dx_deta, dy_deta;
	Cmpnts2** local_coor = new Cmpnts2 * [Lnod]; //local array tp store the coor of the sps in an element, dont really need it just for convenience
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
			}


		for (int k = 0; k < Knod; ++k) { //loop on all sps on the faces to calculate metrics
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
			face_Anorm[el][0][k] = sqrt(face_Acoef[el][0][k]); //||g_2||
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
			face_Anorm[el][1][k] = sqrt(face_Acoef[el][1][k]); //||g_2||
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
			face_Anorm[el][2][k] = sqrt(face_Acoef[el][2][k]); //||g_1||
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
			face_Anorm[el][3][k] = sqrt(face_Acoef[el][3][k]); //||g_1||
		}  //for k=0; k<Knod
	} //for el
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

	if (Lnod_in == 1) t2f = 2 * (1 - i) * (j > 0) + j > 0 + i > 0;   //linear element
		
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
		



	for (unsigned int ti = 0; ti < num_time_steps; ++ti) {
		std::cout << "timestep  " << ti << std::endl;
		solve_advection_diffusion();
		solve_Poisson();
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
	//solves the Poisson's equation for the stream function field

	return 0;
}

char HO_2D::calc_RHS_advection() {
	// ********* This subroutine calculates the - div(VW), V is Cartesian velocity vector, W is vorticity *********
	int ijp, ijpm, eln, neighbor_side;
	double bndr_vel;
	double** local_vort = new double* [Knod]; //local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	for (int i = 0; i < Knod; ++i) local_vort[i] = new double[2];
	double** local_vel  = new double* [Knod]; //local array to hold the velocity along a row of csi [0], and row of eta [1] direction
	for (int i = 0; i < Knod; ++i) local_vel[i] = new double[2];

	double ***bndr_vort, ***bndr_flx, ***upwnd_flx; //quantities per element, per face (west, east, south, north) per sps index
	bndr_vort = new double** [mesh.N_el];
	bndr_flx = new double** [mesh.N_el];
	upwnd_flx = new double** [mesh.N_el];

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
		disc_flx[i] = new double** [2]; //0 for xsi and 1: eta directions
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
			for (int side = 0; side < 4; ++side) bndr_vort[el][side][row_col] = bndr_flx[el][side][row_col] = upwnd_flx[el][side][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_vort[i][0] = vorticity[el][row_col][i];  // xsi direction
				local_vort[i][1] = vorticity[el][i][row_col];  //eta direction
				local_vel[i][0] = velocity_cart[el][row_col][i].x * vol_Dy_Dxsi[el][row_col][i].y - velocity_cart[el][row_col][i].y * vol_Dx_Dxsi[el][row_col][i].y;  //contravariant xsi FLUX component along xsi direction, based on eq. 2.11 in fotis class notes
				local_vel[i][1] = -velocity_cart[el][i][row_col].x * vol_Dy_Dxsi[el][i][row_col].x + velocity_cart[el][i][row_col].y * vol_Dx_Dxsi[el][i][row_col].x;  //contravariant eta FLUX component along eta direction,
			}

			ijp = 0; //ijp=0 (ibnd=idir=0)
			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					for (int m = 0; m < Knod; m++) bndr_vort[el][ijp][row_col] += local_vort[m][idir] * sps_boundary_basis[m][ibnd];
					int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
					if (mesh.elem_neighbor[el].is_on_boundary[elem_side])
						bndr_flx[el][ijp][row_col] = bndr_vort[el][ijp][row_col] * BC_normal_vel[mesh.elem_neighbor[el].boundary_index[elem_side]][row_col];
					else {
						eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
						//value of u* w; ijp = 0: on left boundary at the row_col^th row of sps, ijp = 1 : on right boundary at the row_col^th row of sps, ijp = 2 : on bottom boundary at the row_col^th column of sps, ijp = 3, top at row_col^th column
						for (int m = 0; m < Knod; ++m) bndr_flx[el][ijp][row_col] += local_vel[m][idir] * sps_boundary_basis[m][ibnd]; //flux at the ibnd boundary in idir direction
						bndr_flx[el][ijp][row_col] *= bndr_vort[el][ijp][row_col];
					}
					ijp++;
				}
				//discontinuous flux values on (row_col,1:Knod) solution points for xsi direction and on(row_col, 1:Knod) solution points for eta direction
				for (int i = 0; i < Knod; ++i) disc_flx[el][idir][row_col][i] = local_vel[i][idir] * local_vort[i][idir];
			}
		}
	}

	for (int el = 0; el < mesh.N_el; el++) {
		for (ijp = 0; ijp < 4; ++ijp) {
			ijpm = nbr[ijp]; //neighbor of ijp face: face to the left/south of left/south face or right/north of right/north face
			int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
			if (mesh.elem_neighbor[el].is_on_boundary[elem_side])
				for (int j = 0; j < Knod; ++j)
					upwnd_flx[el][ijp][j] = bndr_flx[el][ijp][j];
			else {
				eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
				for (int j = 0; j < Knod; ++j) {
					if (fabs(bndr_vort[el][ijp][j] - bndr_vort[eln][ijpm][j]) > 1.e-6) {
						bndr_vel = (bndr_flx[el][ijp][j] - bndr_flx[eln][ijpm][j]) / (bndr_vort[el][ijp][j] - bndr_vort[eln][ijpm][j]); //a_tilda
						upwnd_flx[el][ijp][j] = 0.5 * (bndr_flx[el][ijp][j] + bndr_flx[eln][ijpm][j]
							+ sgn[ijp] * fabs(bndr_vel) * (bndr_vort[el][ijp][j] - bndr_vort[eln][ijpm][j]));
					}
					else upwnd_flx[el][ijp][j] = 0.5 * (bndr_flx[el][ijp][j] + bndr_flx[eln][ijpm][j]);
				}
			}
		}
	}

	for (int el = 0; el < mesh.N_el; el++) {

		ijp = 0;
		for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
			for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
				for (int j = 0; j < Knod; ++j) {
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
				for (int m = 0; m < Knod; ++m) BC_values[edge_index][m] = 0.;
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
						ijpm = nbr[ijp];  //local ijp of the neighboring element (eln)
						calc_internal_comm_vals_meth2(el, ijp, ibnd, bndr_vort[el][ijp], bndr_vort[eln][ijpm], bndr_grad_vort[el][ijp], comm_vort[el][ijp]);
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
					int ijpm = nbr[ijp];
					//comm_grad_vort[el][ijp][Knod] is Average of bndr_grad_vort[el][ijp][Knod] + [eln][ijpm][Knod]; same works for Left, South, and North
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
	return 0;
}

void HO_2D::calc_internal_comm_vals_meth2(int el, int ijp, int ibnd, double* vort, double* vort_nbr, double* Dvort, double* com_vort) {
	/*
	vort[0:Knod-1]: in: Knod vorticity values on the ijp internal boundary of a cell which stores extrapolated vorticity from the internal sps nodes on xsi or eta constant sps
	vort_nbr[0:Knod-1]: in: Knod vorticity values on the same local boundary as in vort, but extrapolated from the vorticity of the neighboring sps. so vort and vort_nbr correspond to the same flux points on the boundary, but extrapolated from two separate sides
	Dvort[0:Knod-1]: in: gradient of Knod vorticity at the flux points of ijp side of element el. These gradients are calculated based extrapolation of sps vorticity*derivative of shape functions; out: normal flux to the global boundary on the ijp side
	com_vort[0:Knod-1]: out: The common vorticity (average of left and right (top+bottom) at all flux points on the ijp side of element el
	*/

	double* cross_Dvort = new double[Knod];   //cross derivatives using common values com_vort

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

		//calc the corrected values of grad(Phi) along the mesh interface;
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


}