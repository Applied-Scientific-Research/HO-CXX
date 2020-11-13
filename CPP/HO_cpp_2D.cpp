#include "calculation.h"

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
}

char HO_2D::allocate_arrays() {
	//allocate the memory for the arrays and resize them.
	sps_local_coor = new double[Knod];
	sps_weight = new double[Knod];
	gps_local_coor = new double[mesh.Lnod];
	sps_boundary_basis = new LR_boundary [Knod];
	sps_boundary_grad_basis = new LR_boundary[Knod];

}

void HO_2D::form_bases() {
	//This method froms the Lagrange shape functions: form the sps lagrangian shape function (&derivatives), gps shape functions (&derivatives), form Radau bases
	//--------------------------------------------------------------------------------------------------------------
	// *************************************** in This sub section ***************************************
	// calculate the shape functions of Knod solution points on the local left (-1) and right (+1) boundaries(sps_boundary_basis).
	// calculates the derivative of the Lagrange shape functions of Knod solution points on the local left (-1) and right (+1) boundaries(sps_boundary_grad_basis).
	// calculates the derivative of Knod solution points on the Knod solution points(SolnNodesGradLgrangeBasis)

	double denominator, Grad_numerator[2];
	double* nodes_grad_numerator = new double [Knod];

		SolnBndryLgrangeBasis = 0.d0 !Haji: the value of the shape function of the Knod sps on the left and right boundary(-1, 1))
		SolnBndryGradLgrangeBasis = 0.d0 !Haji: the derivative of shape functions of Knode sps on the left and right boundaries
		SolnNodesGradLgrangeBasis = 0.d0 !Haji : the derivative of shape functions of Knod sps on all other sps.the Knod sps is the first index(say i) and other nodes are second index(say j) : SolnNodesGradLgrangeBasis(i, j)
		DO k = 1, Knod
		SolnBndryLgrangeBasis(k, 0:1) = 1.d0
		Denom = 1.d0
		DO j = 1, Knod
		IF(j.eq.k) CYCLE
		Denom = Denom * (sps(k) - sps(j))       !Basis denominator is common to all evaluations
		!Get the numerators for the extrapolations to L + R
		SolnBndryLgrangeBasis(k, 1) = SolnBndryLgrangeBasis(k, 1) * (1.d0 - sps(j)) !Haji: local right boundary, i.e.xsi = 1
		SolnBndryLgrangeBasis(k, 0) = SolnBndryLgrangeBasis(k, 0) * (-1.d0 - sps(j)) !Haji : local left  boundary, i.e.xsi = -1
		GradNumer = 1.d0  !Haji : 1 sentence in the numerator.The derivative has Knod - 1 sentences.it has 2 indices 0 : left boundary and 1 : right boundary
		NodesGradNumer = 1.d0
		DO i = 1, Knod
		IF(i.eq.k . or .i.eq.j) CYCLE
		!Get the numerators for derivatives of extrpolations to L + R
		GradNumer(1) = GradNumer(1) * (1.d0 - sps(i))
		GradNumer(0) = GradNumer(0) * (-1.d0 - sps(i))
		!Get the numerators for derivatives of interpolation to interior nodes
		NodesGradNumer(1:Knod) = NodesGradNumer(1:Knod) * (sps(1:Knod) - sps(i))
		ENDDO
		SolnBndryGradLgrangeBasis(k, 0:1) = SolnBndryGradLgrangeBasis(k, 0:1) + GradNumer(0:1)
		SolnNodesGradLgrangeBasis(k, 1:Knod) = SolnNodesGradLgrangeBasis(k, 1:Knod) + NodesGradNumer(1:Knod)
		ENDDO
		SolnBndryLgrangeBasis(k, 0:1) = SolnBndryLgrangeBasis(k, 0:1) / Denom
		SolnBndryGradLgrangeBasis(k, 0:1) = SolnBndryGradLgrangeBasis(k, 0:1) / Denom
		!Get grads DOT_PROD[SolnNodesGradLgrangeBasis(k, j), u(k)] at a node / point j using data values u(k) at nodes k
		SolnNodesGradLgrangeBasis(k, 1:Knod) = SolnNodesGradLgrangeBasis(k, 1:Knod) / Denom
		ENDDO





}