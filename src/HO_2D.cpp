/*
 * HO_2D.cpp - Class methods for HO_2D
 *
 * (c)2020-1 Applied Scientific Research, Inc.
 *           Mohammad Hajit <mhajit@gmail.com>
 *           Mark J Stock <markjstock@gmail.com>
 */

#include "HO_2D.hpp"
#include "MathHelpers.hpp"
#include "MemoryHelper.hpp"

#include <cmath>

void HO_2D::release_memory() { //release the memory as destructor

	free_array<double>(vorticity);
	free_array<double>(stream_function);
	free_array<double>(initial_vorticity);
	free_array<Cmpnts2>(velocity_cart);

	delete[] sps_local_coor;
	delete[] sps_weight;
	delete[] gps_local_coor;

	free_array<double>(sps_boundary_basis);
	free_array<double>(sps_boundary_grad_basis);
	free_array<double>(sps_sps_grad_basis);
	free_array<double>(sps_gps_basis);
	free_array<double>(sps_radau);
	free_array<double>(sps_grad_radau);

	free_array<double>(gps_boundary_basis);
	free_array<double>(gps_boundary_grad_basis);
	free_array<double>(gps_sps_basis);
	free_array<double>(gps_sps_grad_basis);

	free_array<Cmpnts2>(vol_Dx_Dxsi);
	free_array<Cmpnts2>(vol_Dy_Dxsi);
	free_array<double>(vol_jac);

    for (int i = 0; i < mesh.N_el; ++i) {
        for (int j = 0; j < Knod; ++j) {
            for (int r = 0; r < Knod; ++r) {
                for (int s = 0; s < 2; ++s) delete[] G[i][j][r][s];
                delete[] G[i][j][r];
            }
            delete[] G[i][j];
        }
        delete[] G[i];
    }
    delete[] G;

    for (int i = 0; i < mesh.N_el; ++i) {
        for (int r = 0; r < 2; ++r) {
            for (int s = 0; s < 2; ++s) {
                for (int j = 0; j < Knod; ++j) delete[] GB[i][r][s][j];
                delete[] GB[i][r][s];
            }
            delete[] GB[i][r];
        }
        delete[] GB[i];
    }
    delete[] GB;

	free_array<double>(face_Acoef);
	free_array<double>(face_Bcoef);
	free_array<double>(face_jac);
	free_array<double>(face_Anorm);

	free_array<Cmpnts2>(face_Dx_Dxsi);
	free_array<Cmpnts2>(face_Dy_Dxsi);

	free_array<double>(RHS_advective);
	free_array<double>(RHS_diffusive);

	free_array<double>(BC_Poisson);
	free_array<double>(BC_advection);
	free_array<double>(BC_diffusion);
	free_array<double>(BC_vorticity);
	free_array<double>(BC_parl_vel);
	free_array<double>(BC_normal_vel);
	free_array<double>(velocity_jump);

	free_array<double>(boundary_source);

	delete[] BC_no_slip;
	delete[] BC_switch_Poisson;
	//delete[] BC_switch_advection;
	delete[] BC_switch_diffusion;

	free_array<Cmpnts2>(BC_cart_vel);

    free_array<double>(k1); //these are used in the RK time integration
    free_array<double>(k2);
    free_array<double>(k3);
    free_array<double>(k4);
    free_array<double>(vort);
}

int HO_2D::read_input_file(const std::string filename) {
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
	  // retrieve number of quadrilateral elements in case of structured mesh (prespecified problem)
	  std::getline(file_handle, temp_string);
	  file_handle >> mesh.N_el_i >>temp_char>> mesh.N_el_j; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');   // ignore until next line

	  // retrieve Knod: number of solution points in each element in each direction
	  std::getline(file_handle, temp_string); //read stoff explaining Knod
	  file_handle >> Knod; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read on Lnod_in: degree of the polynomial of the geometry edges: number of nodes on the elements edges -1
	  getline(file_handle, temp_string);
	  file_handle >> mesh.Lnod_in; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	  mesh.Lnod = mesh.Lnod_in + 1; //number of geometry nodes on each edge

	  // read the name of the mesh file and store it in the mesh.input_msh_file
	  std::getline(file_handle,temp_string);
	  file_handle >> mesh.input_msh_file; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      // read if the mesh file is written in OCC format or buildin GMSH format . 1; means OCC format, zero means original
	  std::getline(file_handle,temp_string);
	  file_handle >> mesh.OCC_format; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  //mesh.input_msh_file = temp_string.c_str();
	  //ignore the next 2 lines as they correspond to the exact geometry for concentric cylinders
	  std::getline(file_handle, temp_string); getline(file_handle, temp_string);

	  // read in the Reynolds number and store the inverse of Reynolds number in Reyn_inv
	  std::getline(file_handle, temp_string);
	  file_handle >> Reyn_inv; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	  Reyn_inv = Reyn_inv > 1.e-3 ? 1. / Reyn_inv : -1.; //if negative then ignore the vicous term

	  // read inmultiplier to expand the grid geometrically by factor dx_ratio
	  std::getline(file_handle, temp_string);
	  file_handle >> mesh.dx_ratio; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in multiplier to make grids randomly non-uniform; uniform: fac=0.
	  std::getline(file_handle, temp_string);
	  file_handle >> mesh.fac; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in HuynhSolver_type
	  std::getline(file_handle, temp_string);
	  file_handle >> HuynhSolver_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in time integration method, 1 = Euler; 2 = Mod.Euler; 4 = RK4
	  std::getline(file_handle, temp_string);
	  file_handle >> time_integration_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in the advective term discretization method, 1 = original discontinuous mass flux field; 2 = Modified continuous mass flux across cell interfaces
	  std::getline(file_handle, temp_string);
	  file_handle >> advection_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in problem type
	  std::getline(file_handle, temp_string);
	  file_handle >> problem_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in total number of time steps
	  std::getline(file_handle, temp_string);
	  file_handle >> num_time_steps; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in time step size
	  std::getline(file_handle, temp_string);
	  file_handle >> dt; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in file saving frequency
	  std::getline(file_handle, temp_string);
	  file_handle >> dump_frequency; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in if it uses large stencil or compact stencil. fast-true means compact
	  std::getline(file_handle, temp_string);
	  file_handle >> fast; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in the type of matrix formation for poisson equation (1=Eigen, 2=Hypre, 3=amgcl)
	  std::getline(file_handle, temp_string);
	  file_handle >> LHS_type; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	  // read in the name of the sample points file
	  std::getline(file_handle, temp_string);
	  file_handle >> sample_points_file; file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

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

void HO_2D::read_process_sample_points() {
	//This method reads the sample points file and process the data inside the file, then finds the element number that contains each sample point
	int N_el = mesh.N_el;
	int Lnod = mesh.Lnod;
	// *********** First read in the sample points coordinates, if it exists ***********
	std::ifstream file_handle(sample_points_file);
	if (file_handle.fail())
		return;

	Cmpnts2 tmp_coor;
	do {
		file_handle >> tmp_coor.x >> tmp_coor.y; //file_handle.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		sample_points_coor.push_back(tmp_coor);
		file_handle >> std::skipws;
	} while (!file_handle.eof());
	tmp_coor.subtract(sample_points_coor.back(), sample_points_coor.rbegin()[1]); // the subtract of last and second to the last element. This is to check if the last element is not repeated
	if (tmp_coor.norm2() < 1e-6) sample_points_coor.pop_back();
	N_sample_points = sample_points_coor.size();
	sample_points_containing_elements.resize(N_sample_points); // the vector to store the element index that contains each sample point
	gps_shapefunc_on_sample_points = new double* [Lnod*Lnod]; //the shape function value of the gps nodes on the sample points
	for (int i = 0; i < Lnod * Lnod; ++i) gps_shapefunc_on_sample_points[i] = new double[N_sample_points];

	//  ***** Now find and save the elements min-max of coordinates and store in min_coor and max_coor *****
	std::vector<Cmpnts2> min_coor(N_el), max_coor(N_el);
	for (int el = 0; el < N_el; el++) {
		max_coor[el] = min_coor[el] = mesh.nodes[mesh.elements[el].nodes[south]].coor; // assume the SW corner node (node index 0) is the min and max (just to fill in values)
		for (int edge_index : mesh.elements[el].edges) { //loop over the 4 edges of each element. doing this because the element may be higher order.
			for (int node_index : mesh.edges[edge_index].nodes) { //loop over the Lnod nodes of the edge edge_index
				if (mesh.nodes[node_index].coor.x < min_coor[el].x) min_coor[el].x = mesh.nodes[node_index].coor.x;
				if (mesh.nodes[node_index].coor.x > max_coor[el].x) max_coor[el].x = mesh.nodes[node_index].coor.x;
				if (mesh.nodes[node_index].coor.y < min_coor[el].y) min_coor[el].y = mesh.nodes[node_index].coor.y;
				if (mesh.nodes[node_index].coor.y > max_coor[el].y) max_coor[el].y = mesh.nodes[node_index].coor.y;
			}
		}
	}

	// **************** Now find the element index that contains each sample point ******************
	std::srand(time(0)); //initiate the randon number generation
	for (int sp = 0; sp < N_sample_points; ++sp) { //loop over all sample points
		bool found_element = false;
		for (int el = 0; el < N_el; ++el) {
			if (sample_points_coor[sp].x < min_coor[el].x || sample_points_coor[sp].x > max_coor[el].x ||
				sample_points_coor[sp].y < min_coor[el].y || sample_points_coor[sp].y > max_coor[el].y) continue; //the point is not inside this element for sure

			// *********** test the casted ray against all sides of the element ***********
			double slope = 2.* static_cast <double> (std::rand()) / static_cast <double> (RAND_MAX) -1.; //generate a random slope between -1 and 1
			int intersections = 0;
			for (int edge_index : mesh.elements[el].edges) { // check if current side intersects with ray. if yes, intersections++;
				for (int segment = 0; segment < mesh.Lnod_in; segment++) { //number of segments is one less than the number of nodes on the edge. This is basically loop over sub-edges on each edge: straight lines
					Cmpnts2 node1 = mesh.nodes[mesh.edges[edge_index].nodes[mesh.tensor2FEM(segment)]].coor;  //node index of the start of the sub-edge (straight line in an edge between two consecutive nodes)
					Cmpnts2 node2 = mesh.nodes[mesh.edges[edge_index].nodes[mesh.tensor2FEM(segment+1)]].coor;//node index of the end of the sub-edge (straight line in an edge between two consecutive nodes)
					intersections += are_intersecting(sample_points_coor[sp], slope, node1, node2);
				}
			}
			if ((intersections & 1) == 1) { //odd number of intersections
				sample_points_containing_elements[sp] = el;
				found_element = true;
				break;
			}
		}

		if (!found_element) { //no element is found to contain the sampling point sp. so loop over all nodes of all elements and take element with the closest node distance to sampling point sp
			tmp_coor.subtract(sample_points_coor[sp], mesh.nodes[mesh.elements[0].nodes[0]].coor);
			double distance = tmp_coor.norm2()+1.;  //initiated distance as distance between sp and node 0 of the element 0 +1. +1 is to make sure it will change in the search process below
			int containing_element = -1;
			for (int el = 0; el < N_el; ++el) {
				for (int i = 0; i < mesh.Lnod * mesh.Lnod; ++i) {
					tmp_coor.subtract(sample_points_coor[sp], mesh.nodes[mesh.elements[el].nodes[i]].coor);
					if (distance > tmp_coor.norm2()) {
						distance = tmp_coor.norm2();
						containing_element = el;
					}
				}
			}
			sample_points_containing_elements[sp] = containing_element;
		}
	}

	// ******* Now find the local coordinate of each sample point in xsi and eta direction, using Newton method *********
	// ******************************************************************************************************************
	std::vector<int> loc_ID;//sweeping the cell indices row by row in the eta direction
	if (mesh.Lnod == 4) loc_ID = { 0, 4, 5, 1, 11, 12, 13, 6, 10, 15, 14, 7, 3, 9, 8, 2 };
	else if (mesh.Lnod == 3) loc_ID = { 0, 4, 1, 7, 8, 5, 3, 6, 2 };
	else if (mesh.Lnod == 2) loc_ID = { 0, 1, 3, 2 };
	std::vector<double> gps_sp_grad_basis_xsi(Lnod), gps_sp_basis_xsi(Lnod), gps_sp_grad_basis_eta(Lnod), gps_sp_basis_eta(Lnod);
	std::vector<Cmpnts2> sample_points_local_coor(N_sample_points); //temporary vector to store the local coordinate (xsi, eta) of the sampling points
	double del_xsi, del_eta; //the increment in xsi and eta in the Newton method
	for (int sp = 0; sp < N_sample_points; ++sp) {
		int el = sample_points_containing_elements[sp];
		double xsi = 0.08;
		double eta = 0.08; //initial guess for the xsi and eta (located at the cell center)

		do {
			for (int k = 0; k < Lnod; ++k) { //find the shape function of all the gps in an element on the sampling points
				gps_sp_grad_basis_xsi[k] = gps_sp_grad_basis_eta[k] = 0.0;
				gps_sp_basis_xsi[k] = gps_sp_basis_eta[k] = 1.;

				double denominator = 1.;
				for (int j = 0; j < Lnod; ++j) {
					if (j == k) continue;
					denominator *= gps_local_coor[k] - gps_local_coor[j];
					gps_sp_basis_xsi[k] *= xsi - gps_local_coor[j];
					gps_sp_basis_eta[k] *= eta - gps_local_coor[j];

					double sp_grad_numerator_xsi = 1.;
					double sp_grad_numerator_eta = 1.;

					for (int i = 0; i < Lnod; ++i) {
						if (i == k || i == j) continue;
						sp_grad_numerator_xsi *= xsi - gps_local_coor[i];
						sp_grad_numerator_eta *= eta - gps_local_coor[i];
					}
					gps_sp_grad_basis_xsi[k] += sp_grad_numerator_xsi;
					gps_sp_grad_basis_eta[k] += sp_grad_numerator_eta;
				}

				gps_sp_grad_basis_xsi[k] /= denominator;
				gps_sp_basis_xsi[k] /= denominator;
				gps_sp_grad_basis_eta[k] /= denominator;
				gps_sp_basis_eta[k] /= denominator;
			}

			double f = sample_points_coor[sp].x;
			double g = sample_points_coor[sp].y;
			double df_dxsi=0.,df_deta=0., dg_dxsi=0., dg_deta = 0.;
			for (int j = 0; j < Lnod; j++) {
				for (int i = 0; i < Lnod; ++i) {
					int gps_index = loc_ID[j*Lnod + i];
					int node_index = mesh.elements[el].nodes[gps_index];
					f -= gps_sp_basis_xsi[i] * gps_sp_basis_eta[j] * mesh.nodes[node_index].coor.x;
					g -= gps_sp_basis_xsi[i] * gps_sp_basis_eta[j] * mesh.nodes[mesh.elements[el].nodes[gps_index]].coor.y;
					df_dxsi -= gps_sp_basis_eta[j] * gps_sp_grad_basis_xsi[i] * mesh.nodes[node_index].coor.x;
					dg_dxsi -= gps_sp_basis_eta[j] * gps_sp_grad_basis_xsi[i] * mesh.nodes[node_index].coor.y;
					df_deta -= gps_sp_grad_basis_eta[j] * gps_sp_basis_xsi[i] * mesh.nodes[node_index].coor.x;
					dg_deta -= gps_sp_grad_basis_eta[j] * gps_sp_basis_xsi[i] * mesh.nodes[node_index].coor.y;
				}
			}

			del_xsi = -1. / (df_dxsi * dg_deta - df_deta * dg_dxsi) * (f * dg_deta - g * df_deta);
			del_eta = -1. / (df_dxsi * dg_deta - df_deta * dg_dxsi) * (g * df_dxsi - f * dg_dxsi);

			xsi += del_xsi;
			eta += del_eta;
		} while (std::fabs(del_xsi) > 0.001 || std::fabs(del_eta) > 0.001);

		sample_points_local_coor[sp].set_coor(xsi, eta);

		// **************** Now that the local coor is found, then find the shape function values of all the gps nodes on this sp point and save it in gps_shapefunc_on_sample_points ***************
		for (int j = 0; j < Lnod; j++) {
			for (int i = 0; i < Lnod; ++i) {
				int gps_index = loc_ID[j * Lnod + i]; //the gps indices are the same as Gmsh indices
				gps_shapefunc_on_sample_points[gps_index][sp] = gps_sp_basis_xsi[i] * gps_sp_basis_eta[j];
			}
		}
	}  //for sp
}

char HO_2D::setup_mesh() {
	//based on the problem type (problem_type variable) creates the proper mesh. if problem_type is 10 then read from msh file
	char success;
	success = mesh.read_msh_file();
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
		g_prime[0] = -0.5 * Knod * (Knod+1.);
		g_prime[1] = 0.5 * Knod * (Knod+1.); //0.5 * minus_one_to_power(Knod) * Knod;
	}

}

char HO_2D::allocate_arrays() {
	//allocate the memory for the arrays and resize them.
	int N_el = mesh.N_el;
	int Lnod = mesh.Lnod;

	initial_vorticity = allocate_array<double>(N_el, Knod, Knod);
	vorticity = allocate_array<double>(N_el, Knod, Knod);
	stream_function = allocate_array<double>(N_el, Knod, Knod);

	velocity_cart = allocate_array<Cmpnts2>(N_el, Knod, Knod);

	sps_local_coor = new double[Knod];
	sps_weight = new double[Knod];
	gps_local_coor = new double[Lnod];

	sps_boundary_basis = allocate_array<double>(Knod, 2);
	sps_boundary_grad_basis = allocate_array<double>(Knod, 2);
	sps_sps_grad_basis = allocate_array<double>(Knod, Knod);

	gps_boundary_basis = allocate_array<double>(Lnod, 2);
	gps_boundary_grad_basis = allocate_array<double>(Lnod, 2);
	gps_sps_basis = allocate_array<double>(Lnod, Knod);
	sps_gps_basis = allocate_array<double>(Knod, Lnod);
	gps_sps_grad_basis = allocate_array<double>(Lnod, Knod);

	sps_radau = allocate_array<double>(Knod, 2);
	sps_grad_radau = allocate_array<double>(Knod, 2);

	vol_Dx_Dxsi = allocate_array<Cmpnts2>(N_el, Knod, Knod);
	vol_Dy_Dxsi = allocate_array<Cmpnts2>(N_el, Knod, Knod);
	vol_jac = allocate_array<double>(N_el, Knod, Knod);

	G = new double**** [N_el];
	GB = new double**** [N_el];
	for (int i = 0; i < N_el; ++i) {
		G[i] = new double*** [Knod];
		GB[i] = new double*** [2]; //for the direction d
		for (int j = 0; j < Knod; ++j) {
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

	face_Acoef = allocate_array<double>(N_el, 4, Knod);
	face_Bcoef = allocate_array<double>(N_el, 4, Knod);
	face_jac   = allocate_array<double>(N_el, 4, Knod);
	face_Anorm = allocate_array<double>(N_el, 4, Knod);
	face_Dx_Dxsi = allocate_array<Cmpnts2>(N_el, 4, Knod);
	face_Dy_Dxsi = allocate_array<Cmpnts2>(N_el, 4, Knod);

	RHS_advective = allocate_array<double>(N_el, Knod, Knod);
	RHS_diffusive = allocate_array<double>(N_el, Knod, Knod);

	BC_no_slip = new bool[mesh.N_Gboundary]; //Set NoSlip to all walls (as default)

	BC_switch_Poisson = new unsigned char[mesh.N_edges_boundary];
	//BC_switch_psi = new unsigned char[mesh.N_edges_boundary];
	//BC_switch_advection = new unsigned char[mesh.N_edges_boundary];
	BC_switch_diffusion = new unsigned char[mesh.N_edges_boundary];

	boundary_source = allocate_array<double>(mesh.N_edges_boundary, Knod*Knod, Knod);

	BC_Poisson = allocate_array<double>(mesh.N_edges_boundary, Knod);
	//BC_psi = allocate_array<double>(mesh.N_edges_boundary, Knod);
	BC_advection = allocate_array<double>(mesh.N_edges_boundary, Knod);
	BC_diffusion = allocate_array<double>(mesh.N_edges_boundary, Knod);

	BC_vorticity = allocate_array<double>(mesh.N_edges_boundary, Knod);
	velocity_jump = allocate_array<double>(mesh.N_edges_boundary, Knod);
	BC_parl_vel = allocate_array<double>(mesh.N_edges_boundary, Knod);
	BC_normal_vel = allocate_array<double>(mesh.N_edges_boundary, Knod);

	BC_cart_vel = allocate_array<Cmpnts2>(mesh.N_edges_boundary, Knod);

    k1 = allocate_array<double>(N_el, Knod, Knod); //these are used in the RK time integration
    k2 = allocate_array<double>(N_el, Knod, Knod);
    k3 = allocate_array<double>(N_el, Knod, Knod);
    k4 = allocate_array<double>(N_el, Knod, Knod);
    vort = allocate_array<double>(N_el, Knod, Knod);


	return 0;
}

void HO_2D::form_bases() {
	//This method froms the Lagrange shape functions: form the sps lagrangian shape function (&derivatives), gps shape functions (&derivatives), form Radau bases
	//--------------------------------------------------------------------------------------------------------------
	// *************************************** in This sub section ***************************************
	// calculate the value of shape functions of Knod solution points on the local left (-1) and right (+1) boundaries(sps_boundary_basis).
	// calculates the derivative of the Lagrange shape functions of Knod solution points on the local left (-1) and right (+1) boundaries(sps_boundary_grad_basis).
	// calculates the derivative of Knod solution points on the Knod solution points(sps_sps_grad_basis). the Knod sps is the first index(say i) and other nodes are second index(say j) : sps_sps_grad_basis(i, j)

	int Lnod = mesh.Lnod;
	double denominator;
	double grad_numerator[2];
	double* sps_grad_numerator = new double [Knod];

	for (int k = 0; k < Knod; ++k) {
		sps_boundary_basis[k][0] = sps_boundary_basis[k][1] = 1.;
		sps_boundary_grad_basis[k][0] = sps_boundary_grad_basis[k][1] = 0.0;
		for (int m = 0; m < Knod; ++m) sps_sps_grad_basis[k][m] = 0.;
		for (int m = 0; m < Lnod; ++m) sps_gps_basis[k][m] = 1.0;
		denominator = 1.;

		for (int j = 0; j < Knod; ++j) {
			if (j == k) continue;
			denominator *= sps_local_coor[k] - sps_local_coor[j];
			sps_boundary_basis[k][0] *= -1.0 - sps_local_coor[j]; //local left  boundary, i.e.xsi = -1
			sps_boundary_basis[k][1] *= 1.0 - sps_local_coor[j]; //local right  boundary, i.e.xsi = +1
			grad_numerator[0] = grad_numerator[1] = 1.;  //1 sentence in the numerator. The derivative has Knod-1 sentences. it has 2 indices 0:left boundary and 1 :right boundary
			for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] = 1.;
			for (int m=0; m<Lnod; ++m) sps_gps_basis[k][m] *= gps_local_coor[m] - sps_local_coor[j];

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
		for (int m = 0; m < Lnod; ++m) sps_gps_basis[k][m] /= denominator;
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

	/* in this version, the Cartesian velocity BC (BC_cart_vel) are replaced with the slip and normal components.
	Slip velocity: is always in the CCW to the boundary edges. For instance, if the top lid of cavity goes to the right, then the slip velocity is -1
	Normal velocity: always measured in the outward to the boundary. So, incoming flux is negative and outgoing flux is positive

	Three problems are specified. a general problem specification is also possible with setting the right BCs.
	Problem_type=1 : Cavity flow
	problem_type=2: Backward facing step (BFS)
	problem_type=3: Flow over a cylinder
	problem_type=11: Hybrid run driven by Omega2D
	other problems: change the proper boundary conditions as required.
	*/

	for (int el = 0; el < mesh.N_el; ++el) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				initial_vorticity[el][j][i] = 0.;
				vorticity[el][j][i] = initial_vorticity[el][j][i];
				stream_function[el][j][i] = 0.;
			}
		}
	}

	if (problem_type==1) { //cavity flow: the mass flux is zero and the slip velocity is -1 on the top plane (lid moves to the right, so CCW)
        	for (int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b) {
            		for (int j = 0; j < Knod; ++j) {
                		BC_Poisson[el_b][j] = 0.; //psi=0
                		BC_normal_vel[el_b][j] = 0.;
                		BC_parl_vel[el_b][j] = 0.;
            		}
        	}

		// std::fill works if all the data is contiguous in memory, which is currently is not
		std::fill(BC_no_slip, BC_no_slip + mesh.N_Gboundary, true); // if a boundary is no-slip wall or not.
		std::fill(BC_switch_Poisson, BC_switch_Poisson + mesh.N_edges_boundary, DirichletBC);
		//std::fill(&BC_Poisson[0][0], &BC_Poisson[0][0] + mesh.N_edges_boundary * Knod, 0.); //psi=0
		std::fill(BC_switch_diffusion, BC_switch_diffusion + mesh.N_edges_boundary, NeumannBC);
		//std::fill(&BC_normal_vel[0][0], &BC_normal_vel[0][0] + mesh.N_edges_boundary * Knod, 0.);
		//std::fill(&BC_parl_vel[0][0], &BC_parl_vel[0][0] + mesh.N_edges_boundary * Knod, 0.);

		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
			if (mesh.boundaries[G_boundary].name == "top") { // the proper boundary is selected
				for (int edge_index : mesh.boundaries[G_boundary].edges) { //loop over all the edges forming the top boundary
					for (int i = 0; i < Knod; ++i) { // loop over all sps on each edge on the proper boundary
						BC_parl_vel[edge_index][i] = -1.; //the lid moves CW on the global boundary
					}
				}
			}
		}
	}

    else if (problem_type==2) {  //Backward Facing step
        //BFS: inlet: specific normal vel ---> calculate psi based on the mass flux
        double y_min = 0, y_max = 0;
        for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) { //find the min and max y at the inlet
            if (mesh.boundaries[G_boundary].name != "inlet") continue;
            y_min = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
            y_max = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
            for (int edge_index: mesh.boundaries[G_boundary].edges) { //loop over all the edges to find min y and max y of the inlet boundary
                y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_min); // compare with the first node of the edge
                y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_min); // compare with the last  node of the edge
                y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_max); // compare with the first node of the edge
                y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_max); // compare with the last  node of the edge
            }
        }

        for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
            if (mesh.boundaries[G_boundary].name == "inlet") { //vorticity, psi and normal velocity are defined
                BC_no_slip[G_boundary] = false;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_diffusion[edge_index] = DirichletBC;
                    BC_switch_Poisson[edge_index] = DirichletBC;

                    for (int k=0; k<Knod; ++k) {
                        double y=0.;
                        for (int j=0; j<mesh.Lnod; ++j) y += gps_sps_basis[j][k] * mesh.nodes[ mesh.edges[edge_index].nodes[mesh.tensor2FEM(j)] ].coor.y; //y at k sps on the edge_index boundary edge

                        if (false) { // uniform velocity profile = 1 into the domain
                            BC_Poisson[edge_index][k] = y - y_min;
                            BC_diffusion[edge_index][k] =0.; //zero vorticity
                            BC_normal_vel[edge_index][k] = -1.; //In the Fortran code, it is considered as the volume flux but here it is the normal out of domain VELOCITY, to be consistent with the BC_parl_vel
                        }
                        else { //parabolic inlet profile
                            BC_Poisson[edge_index][k] = pow(y - y_min,2.) * (3.*y_max-y_min-2.*y)/ pow(y_max-y_min,2.); //psi BC
                            BC_diffusion[edge_index][k] = 6.*(2.*y-y_min-y_max)/ pow(y_max-y_min,2.);  //vorticity BC
                            BC_normal_vel[edge_index][k] = 6.*(y-y_min)*(y-y_max)/pow(y_max-y_min,2.);
                        }
                    }
                }
            }  //case inlet

            else if (mesh.boundaries[G_boundary].name == "outlet") {  //Neumann for vorticity, Neumann for psi, and slip velocity are known
                // set slip BC to not allow diffusion
                BC_no_slip[G_boundary] = false;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_diffusion[edge_index] = NeumannBC;
                    BC_switch_Poisson[edge_index] = NeumannBC;
                    for (int k=0; k<Knod; ++k) {
                        BC_diffusion[edge_index][k] = 0.;  //zero domega/dn
                        BC_Poisson[edge_index][k] = 0.;
                    }
                }
            }

            else if (mesh.boundaries[G_boundary].name == "bottom") {  //psi is given on the bottom wall. domega/dn is derived in the code while running, and both the tangential and normal velocities are known
                BC_no_slip[G_boundary] = true;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_Poisson[edge_index] = DirichletBC;
                    BC_switch_diffusion[edge_index] = NeumannBC; //BC_diffusion is calculated in the code during runtime when BC_no_slip is true
                    for (int k=0; k<Knod; ++k) {
                        BC_Poisson[edge_index][k] = 0.;
                        BC_normal_vel[edge_index][k] = 0.;
                        BC_parl_vel[edge_index][k] = 0.; // it is used to calculate the normal derivative of vorticity
                    }
                }
            }

            else if (mesh.boundaries[G_boundary].name == "top") {
                BC_no_slip[G_boundary] = true;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_Poisson[edge_index] = DirichletBC;
                    BC_switch_diffusion[edge_index] = NeumannBC; //BC_diffusion is calculated in the code during runtime
                    for (int k=0; k<Knod; ++k) {
                        BC_Poisson[edge_index][k] = y_max - y_min; //the average velocity is 1, so the volume flux in equal to the area of the inlet = psi(top) - psi(bottom)
                        BC_normal_vel[edge_index][k] = 0.;
                        BC_parl_vel[edge_index][k] = 0.; // it is used to calculate the normal derivative of vorticity
                    }
                }
            }
        }  //for int Gboundary
    }  //else if problem_type ==2

    else if (problem_type==3) { //flow over an object
        double y_min = 0, y_max = 0;
        for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) { //find the min and max y at the inlet
            if (mesh.boundaries[G_boundary].name != "inlet") continue;
            y_min = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
            y_max = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
            for (int edge_index: mesh.boundaries[G_boundary].edges) { //loop over all the edges to find min y and max y of the inlet boundary
                y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_min); // compare with the first node of the edge
                y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_min); // compare with the last  node of the edge
                y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_max); // compare with the first node of the edge
                y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_max); // compare with the last  node of the edge
            }
        }

	    for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
	        if (mesh.boundaries[G_boundary].name == "inlet") { //vorticity, psi and normal velocity are defined
	            BC_no_slip[G_boundary] = false;
	            for (int edge_index : mesh.boundaries[G_boundary].edges) {
	                BC_switch_diffusion[edge_index] = DirichletBC;
	                BC_switch_Poisson[edge_index] = DirichletBC;
	                for (int k=0; k<Knod; ++k) {
	                    double y=0.;
	                    for (int j=0; j<mesh.Lnod; ++j) y += gps_sps_basis[j][k] * mesh.nodes[ mesh.edges[edge_index].nodes[mesh.tensor2FEM(j)] ].coor.y; //y at k sps on the edge_index boundary edge
	                    BC_Poisson[edge_index][k] = y - y_min; // uniform velocity profile = 1 into the domain
	                    BC_diffusion[edge_index][k] =0.; //zero vorticity
	                    BC_normal_vel[edge_index][k] = -1.; //In the Fortran code, it is considered as the volume flux but here it is the normal out of domain VELOCITY, to be consistent with the BC_parl_vel
	                }
	            }
	        }  //case inlet

            else if (mesh.boundaries[G_boundary].name == "outlet") {  //Neumann for vorticity, Neumann for psi, and slip velocity are known
                // set slip BC to not allow diffusion
                BC_no_slip[G_boundary] = false;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_diffusion[edge_index] = NeumannBC;
                    BC_switch_Poisson[edge_index] = NeumannBC;
                    for (int k=0; k<Knod; ++k) {
                        BC_diffusion[edge_index][k] = 0.;  //zero domega/dn
                        BC_Poisson[edge_index][k] = 0.;
                    }
                }
            }

            else if (mesh.boundaries[G_boundary].name == "top") {
                BC_no_slip[G_boundary] = false;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_Poisson[edge_index] = DirichletBC;
                    BC_switch_diffusion[edge_index] = NeumannBC; //BC_diffusion is calculated in the code during runtime
                    for (int k=0; k<Knod; ++k) {
                        BC_Poisson[edge_index][k] = y_max - y_min; //the average velocity is 1, so the volume flux in equal to the area of the inlet = psi(top) - psi(bottom)
                        BC_normal_vel[edge_index][k] = 0.;
                        BC_parl_vel[edge_index][k] = -1.; // goes to right, so CW. it is used to calculate the normal derivative of vorticity
                    }
                }
            }

            else if (mesh.boundaries[G_boundary].name == "bottom") {
                BC_no_slip[G_boundary] = false;
                for (int edge_index : mesh.boundaries[G_boundary].edges) {
                    BC_switch_Poisson[edge_index] = DirichletBC;
                    BC_switch_diffusion[edge_index] = NeumannBC; //BC_diffusion is calculated in the code during runtime
                    for (int k=0; k<Knod; ++k) {
                        BC_Poisson[edge_index][k] = 0.; //the average velocity is 1, so the volume flux in equal to the area of the inlet = psi(top) - psi(bottom)
                        BC_normal_vel[edge_index][k] = 0.;
                        BC_parl_vel[edge_index][k] = 1.; // Goes to the right, so CCW. it is used to calculate the normal derivative of vorticity
                    }
                }
            }

			else if (mesh.boundaries[G_boundary].name == "object") {
				BC_no_slip[G_boundary] = true; // if a boundary is no-slip wall or not.
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					BC_switch_Poisson[edge_index] = DirichletBC;
					BC_switch_diffusion[edge_index] = NeumannBC;  //calculate domega/dn during runtime
					for (int k=0; k<Knod; ++k) {
						BC_Poisson[edge_index][k] = 0.5*(y_max-y_min); //psi=0
						BC_normal_vel[edge_index][k] = 0.;
						BC_parl_vel[edge_index][k] = 0.;
					}
                }
            }
        } //for Gboundary
    } //if problem_type==3

	else if (problem_type==11) { //hybrid run

		// first, set switches for entire boundaries
		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {

			if (mesh.boundaries[G_boundary].name == "open") {
				BC_no_slip[G_boundary] = false;
				for (auto edge_index : mesh.boundaries[G_boundary].edges) {
					//std::cout << "open edge index " << edge_index << std::endl;
					BC_switch_Poisson[edge_index] = NeumannBC;
					BC_switch_diffusion[edge_index] = DirichletBC;
					//BC_switch_advection[edge_index] = DirichletBC;
				}

			} else if (mesh.boundaries[G_boundary].name == "wall") {
				BC_no_slip[G_boundary] = true;
				for (auto edge_index : mesh.boundaries[G_boundary].edges) {
					//std::cout << "wall edge index " << edge_index << std::endl;
					BC_switch_Poisson[edge_index] = DirichletBC;
					BC_switch_diffusion[edge_index] = NeumannBC;
					//BC_switch_advection[edge_index] = DirichletBC;
				}

			} else if (mesh.boundaries[G_boundary].name == "slipwall") {
				assert(false && "HO_2D::setup_IC_BC_SRC - slipwall BC unsupported");
				BC_no_slip[G_boundary] = false;

			} else if (mesh.boundaries[G_boundary].name == "inlet") {
				assert(false && "HO_2D::setup_IC_BC_SRC - inlet BC unsupported");
				BC_no_slip[G_boundary] = false;

			} else if (mesh.boundaries[G_boundary].name == "outlet") {
				assert(false && "HO_2D::setup_IC_BC_SRC - outlet BC unsupported");
				BC_no_slip[G_boundary] = false;

			} else {
				assert(false && "HO_2D::setup_IC_BC_SRC - unknown BC in mesh");
			}
		}
    } //if problem_type==11

    else {
		//for any arbitrary problem, specify the BC's of the problem here:
		//  modify the below terms as the problem is designed

        std::fill(BC_no_slip, BC_no_slip + mesh.N_Gboundary, true); // if a boundary is no-slip wall or not.
        std::fill(BC_switch_Poisson, BC_switch_Poisson + mesh.N_edges_boundary, DirichletBC);
        std::fill(&BC_Poisson[0][0], &BC_Poisson[0][0] + mesh.N_edges_boundary * Knod, 0.); //psi=0
        std::fill(BC_switch_diffusion, BC_switch_diffusion + mesh.N_edges_boundary, NeumannBC);
        std::fill(&BC_normal_vel[0][0], &BC_normal_vel[0][0] + mesh.N_edges_boundary * Knod, 0.);
        std::fill(&BC_parl_vel[0][0], &BC_parl_vel[0][0] + mesh.N_edges_boundary * Knod, 0.);
    }

    // Now convert the outward normal directions into the local h direction
    for (int el_b=0; el_b<mesh.N_edges_boundary; ++el_b) {
        //int element_index = mesh.boundary_elem_ID[el_b].element_index; //index of the element that has a face on the global boundary
		const int elem_side = mesh.boundary_elem_ID[el_b].side; //side of the element that has the edge on the global boundary
		const int ijp = f2i[elem_side];
        const double sign = 2.*(ijp%2) - 1.;
        for (int k=0; k<Knod; ++k) {
            BC_normal_vel[el_b][k] *= sign;
            if (BC_switch_diffusion[el_b]==NeumannBC) BC_diffusion[el_b][k] *= sign;
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
			!Haji : The Face_Acoef and Face_Bcoef, forms the covariant metric tensor g_ij = [g11 g12; g12, g22]] (symmetric matrix))\

	G[el][j][i]][d1][d2] is (g_d1 dot g_d2)/J at j,i of element el, so the size is[N_el][Knod][[Knod][2][2]. here g_0=(-partial(x)/partial(eta),partial(y)/partial(eta)), g_1=(partial(x)/partial(csi),-partial(y)/partial(csi))
	GB[el][d][t]][j][d2] is (g_d dot g_d2)/J at the jth flux point on the d direction , t side of the element el, so the size is[N_el][2][[2][Knod][2]. here g_0=(-partial(x)/partial(eta),partial(y)/partial(eta)), g_1=(partial(x)/partial(csi),-partial(y)/partial(csi)) on the face

	*/


	int Lnod = mesh.Lnod;
	double dx_dxsi, dy_dxsi, dx_deta, dy_deta;
	Cmpnts2 g[2]; //to store g_0=(-partial(x)/partial(eta),partial(y)/partial(eta)) and g_1=(partial(x)/partial(csi),-partial(y)/partial(csi))
	Cmpnts2** local_coor = new Cmpnts2 * [Lnod]; //local array tp store the coor of the gps in an element, dont really need it just for convenience
	for (int i = 0; i < Lnod; ++i) local_coor[i] = new Cmpnts2 [Lnod];

	for (int el = 0; el < mesh.N_el; ++el) {
		for (int j = 0; j < Lnod; ++j)
			for (int i = 0; i < Lnod; ++i)
				local_coor[j][i] = mesh.nodes[mesh.elements[el].nodes[mesh.tensor2FEM(i, j)]].coor;

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
			face_Dx_Dxsi[el][0][k].x = dx_dxsi; face_Dx_Dxsi[el][0][k].y = dx_deta;
			face_Dy_Dxsi[el][0][k].x = dy_dxsi; face_Dy_Dxsi[el][0][k].y = dy_deta;
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
			face_Dx_Dxsi[el][1][k].x = dx_dxsi; face_Dx_Dxsi[el][1][k].y = dx_deta;
			face_Dy_Dxsi[el][1][k].x = dy_dxsi; face_Dy_Dxsi[el][1][k].y = dy_deta;
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
			face_Dx_Dxsi[el][2][k].x = dx_dxsi; face_Dx_Dxsi[el][2][k].y = dx_deta;
			face_Dy_Dxsi[el][2][k].x = dy_dxsi; face_Dy_Dxsi[el][2][k].y = dy_deta;
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
			face_Dx_Dxsi[el][3][k].x = dx_dxsi; face_Dx_Dxsi[el][3][k].y = dx_deta;
			face_Dy_Dxsi[el][3][k].x = dy_dxsi; face_Dy_Dxsi[el][3][k].y = dy_deta;
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

char HO_2D::solve_vorticity_streamfunction() {

	const bool debug = false;

	//solves the two set of equations: 1) advection diffusion of vorticity + 2) Poisson eq for streamfunction
	for (int el = 0; el<mesh.N_el; ++el)
		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i) {
				vorticity[el][j][i] = initial_vorticity[el][j][i];
			}


	// for validation of cases, it exports the quantities and calculates the error norm versus analytical results.
	if (debug) save_output(0);

	// perform time integration, note ti is global step count
	for (ti = 0; ti <= num_time_steps; ++ti) {
		std::cout << "timestep  " << ti << std::endl;

		// do all the work
		solve_advection_diffusion();

		// for debugging purposes and testing only
		if (debug) {
			//temporarily set the right hand side (vorticity) of the poisson equation
			for (int el = 0; el < mesh.N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
				Cmpnts2 coor(0., 0.); //coordinate of sps[j][i]
				for (int ny = 0; ny < mesh.Lnod; ++ny) for (int nx = 0; nx < mesh.Lnod; ++nx) {
					const int node_index = mesh.elements[el].nodes[mesh.tensor2FEM(nx, ny)];
					coor.plus(gps_sps_basis[ny][j] * gps_sps_basis[nx][i], mesh.nodes[node_index].coor);
				}
				//double ym = M_PI * (1. - coor.y);
				//vorticity[el][j][i] = M_PI * std::sin(M_PI * coor.x) * (ym * (1. + 2. * std::cos(ym)) + 2. * std::sin(ym));
				//vorticity[el][j][i] = std::sin(M_PI * coor.y);
				//vorticity[el][j][i] = M_PI*std::sin(M_PI*coor.x)*(2.*std::sin(M_PI*coor.y) + 2.*M_PI*coor.y*std::cos(M_PI*coor.y) + M_PI*coor.y);
				vorticity[el][j][i] = -4. * (coor.x * coor.x + coor.y * coor.y + 1.) * std::exp(coor.x * coor.x + coor.y * coor.y);
				//vorticity[el][j][i] = -4.;
			}
		}

		//if you would like to calculate the streamfunction at n+1 time step
		//solve_Poisson(vorticity);

		if (!((ti+1) % dump_frequency)) {
			//save_output(ti);
			save_vorticity_vtk(ti+1);
			save_smooth_vtk(ti+1);
		}
 	}

	return 0;
}

char HO_2D::solve_advection_diffusion() {
	//solves the advection diffusion eq. for the vorticity field: VTE

	const int N_el = mesh.N_el; 

	//static bool are_allocated = false;
	//static int orig_nel = N_el, orig_knod = Knod;
	//static double*** k1, ***k2, ***k3, ***k4, ***vort;

/*
	if (not are_allocated) {
		// if anything changed since last allocation, free first
		if (N_el != orig_nel or Knod != orig_knod) {
			free_array(k1);
			free_array(k2);
			free_array(k3);
			free_array(k4);
			free_array(vort);
		}

		// allocate here
		k1 = allocate_array<double>(N_el, Knod, Knod);
		k2 = allocate_array<double>(N_el, Knod, Knod);
		k3 = allocate_array<double>(N_el, Knod, Knod);
		k4 = allocate_array<double>(N_el, Knod, Knod);
		vort = allocate_array<double>(N_el, Knod, Knod);

		// and remember what we did here
		orig_nel = N_el;
		orig_knod = Knod;
		are_allocated = true;
	}
*/

	// first order Euler time integration method
	if (time_integration_type == 1) {

		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vorticity, k1, 0.);

		//vorticity at time step n+1
		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vorticity[el][j][i] += k1[el][j][i];
		}
	}

	// explicit RK2 time integration
	else if (time_integration_type == 2) {

		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vorticity, k1, RK2_coeff[0]);

		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vort[el][j][i] = vorticity[el][j][i] + RK2_coeff[1] * k1[el][j][i];
		}
		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vort, k2, RK2_coeff[1]);

		//calc vorticity at time step n+1
		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vorticity[el][j][i] += 0.5 * (k1[el][j][i] + k2[el][j][i]);
		}
	}

	//explicit RK4 time integration https://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
	else if (time_integration_type == 4) {

		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vorticity, k1, RK4_coeff[0]);

		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vort[el][j][i] = vorticity[el][j][i] + RK4_coeff[1] * k1[el][j][i];
		}
		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vort, k2, RK4_coeff[1]);

		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vort[el][j][i] = vorticity[el][j][i] + RK4_coeff[2] * k2[el][j][i];
		}
		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vort, k3, RK4_coeff[2]);

		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vort[el][j][i] = vorticity[el][j][i] + RK4_coeff[3] * k3[el][j][i];
		}
		//calculate the RHS terms for the advective term and diffusive terms, i.e. -div(vw) and Laplacian(w)
		Euler_time_integrate(vort, k4, RK4_coeff[3]);

		//calc vorticity at time step n+1
		for (int el = 0; el < N_el; ++el) for (int j = 0; j < Knod; ++j) for (int i = 0; i < Knod; ++i) {
			vorticity[el][j][i] += (k1[el][j][i] + 2.*k2[el][j][i] + 2.*k3[el][j][i] + k4[el][j][i]) / 6.;
		}
	}

    // do not delete memory, but if we wanted to, we would do this
    //free_array<double>(k1);
    //free_array<double>(k2);
    //free_array<double>(k3);
    //free_array<double>(k4);
    //free_array<double>(vort);

	return 0;
}

void HO_2D:: form_Laplace_operator_matrix() {
	// This subroutine forms the left hand side matrix derived form the Laplace discretization. The matrix is sparse and in Eigen format
	// The type of BC for the psi are defined in BC_switch_psi which is specified in setup_IC_BC_SRC
	int N_el = mesh.N_el;
	int eln, ijpm; //neighbor element; element side, neighbor element side
	int dn, tn; // the direction and boundary side of the neighboring element eln, adjacent to d direction and t side of the element el
	double coeff, tmp2, tmp3;
	const int Ksq = Knod * Knod;
	//unsigned int Km1 = Knod - 1;
	const int K4 = Ksq * Ksq;
	//The coefficients in the left hand side matrix that has contribution from the element itself. it has the coefficients for element el, sps ij=j*Knod+i and rs=r*Knod+s, so laplacian_center[el][ij][rs]
	double*** laplacian_center = allocate_array<double>(N_el, Ksq, Ksq);
	//The coefficients in the left hand side matrix that has contribution from the 4 neighboring element (if are not located on the global boundary). it has the coefficients for element el per cell side, sps ij=j*Knod+i and rs=r*Knod+s of the neighbor element, so laplacian_center[el][4][ij][rs]
	double**** laplacian_neighbor = allocate_array<double>(N_el, 4, Ksq, Ksq);
	//double** NeuMatrix_Orig = new double* [Knod], **NeuMatrix = new double* [Knod]; //small dense matrix to obtain comx(per element) for a given Neumann BC
	Eigen::MatrixXd NeuMatrix_Orig(Knod, Knod), NeuMatrix(Knod, Knod); //I am using Eigen here to save energy and time

	//std::fill(BC_switch, BC_switch + mesh.N_edges_boundary, /*NeumannBC*/ DirichletBC); //to test the poisson solver

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
				const int ij = j * Knod + i; //local cumulative sps index in element k

				for (int m = 0; m < Knod; ++m) {
					for (int r = 0; r < 2; r++) {
						const int a = m * (1 - r) + i * r;
						const int b = m * r + j * (1 - r);
						tmp3 = 0.;
						for (int s = 0; s < 2; s++) tmp3 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s]; //from the M,N,P,Q of ADrins note
						double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
						for (int d = 0; d < 2; ++d) {
							const double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
							const int tmp1 = a * (1 - d) + b * d; //the second index for SNGLB
							for (int p = 0; p < Knod; ++p) {
								const int j1 = p * d + b * (1 - d);
								const int i1 = p * (1 - d) + a * d;
								const int rs = j1 * Knod + i1;
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
				const int ijp = 2 * d + t;
				const cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir

				if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) { // if this face is on the global boundary:
					const int edge_index = mesh.elem_neighbor[el].boundary_index[elem_side];  // the index of the edge that is on global boundary, which is on the d direction and t side of element el
					assert(edge_index < mesh.N_edges_boundary && "boundary_index exceeds array bounds");

					if (BC_switch_Poisson[edge_index] == DirichletBC) { // Dirichlet boundary condition
						//******************************* The dirichlet BC case: *********************************
						for (int j = 0; j < Knod; ++j) {
							for (int i = 0; i < Knod; ++i) {
								const int ij = j * Knod + i; //local cumulative sps index in element k
								for (int m = 0; m < Knod; ++m) {
									for (int r = 0; r < 2; ++r) {
										const int a = m * (1 - r) + i * r;
										const int b = m * r + j * (1 - r);
										const int q = a * d + b * (1 - d);
										double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
										double NGR1 = sps_grad_radau[a * (1 - d) + b * d][t];
										double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
										tmp2 = 0.;
										for (int s = 0; s < 2; s++) tmp2 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s];
										for (int p = 0; p < Knod; ++p) {
											const int j1 = p * d + b * (1 - d);
											const int i1 = p * (1 - d) + a * d;
											const int rs = j1 * Knod + i1; //cumulative index
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
								const int ij = j * Knod + i; //local cumulative sps index in element k
								const int alpha = i * d + j * (1 - d);
								const double A = GB[el][d][t][alpha][d];  //[(g_d . g_d)/J]|k,d,t,alpha
								const double NGR1 = sps_grad_radau[j * d + i * (1 - d)][t];
								for (int m = 0; m < Knod; ++m) {
									const int j1 = m * d + j * (1 - d);
									const int i1 = m * (1 - d) + i * d;
									const int rs = j1 * Knod + i1;
									laplacian_center[el][ij][rs] += A * (sps_boundary_grad_basis[m][t] - g_prime[t] * sps_boundary_basis[m][t]) * NGR1; //page 3 of my notes 3 Feb 2021
								}
								boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * alpha - 1) + alpha] += A * g_prime[t] * NGR1; //page 3 of my notes 3 Feb 2021
								const double B = GB[el][d][t][alpha][1 - d];
								for (int m = 0; m < Knod; ++m) {
									boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * m - 1) + m] += B * sps_sps_grad_basis[m][alpha] * NGR1; //page 3 of my notes 3 Feb 2021
								}
							}
						}
					} // if the edge is located on the Dirichlet-type global boundary

					else if (BC_switch_Poisson[edge_index] == NeumannBC) { //if the edge is located on Neumann type boundary
						//******************************* The Neumann BC case: *********************************
						// find NeuMatrix = inverse(NeuMatrix_Orig)
						for (int i = 0; i < Knod; ++i) { // row index
							for (int j = 0; j < Knod; ++j) //column index of the matrix
								NeuMatrix_Orig(i,j) = -GB[el][d][t][i][1 - d] * sps_sps_grad_basis[j][i];
							NeuMatrix_Orig(i,i) += -GB[el][d][t][i][d] * g_prime[t];
						}

						NeuMatrix = NeuMatrix_Orig.inverse();

						const double sign = t ? 1. : -1.;

						for (int j = 0; j < Knod; ++j) {
							for (int i = 0; i < Knod; ++i) {
								const int ij = j * Knod + i; //local cumulative sps index in element k

								const int alpha = i * d + j * (1 - d);
								boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * alpha - 1) + alpha] += sign * sps_grad_radau[j * d + i * (1 - d)][t]
									 * std::sqrt(GB[el][d][t][alpha][d]*face_jac[el][ijp][alpha]); //page 9 on my notes 17 Feb 2021, effect of NBC, according to the edge[edge_index] direction

								for (int m = 0; m < Knod; ++m) {
									for (int r = 0; r < 2; r++) {
										const int a = m * (1 - r) + i * r;
										const int b = m * r + j * (1 - r);
										const int q = a * d + b * (1 - d);
										double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
										double NGR1 = sps_grad_radau[a * (1 - d) + b * d][t];
										double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
										tmp2 = 0.;
										for (int s = 0; s < 2; s++) tmp2 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s];
										for (int p = 0; p < Knod; ++p) {
											const int j1 = p * d + b * (1 - d);
											const int i1 = p * (1 - d) + a * d;
											const int rs = j1 * Knod + i1; //cumulative index
											laplacian_center[el][ij][rs] -= 0.5 * SNGLB1 * A * sps_boundary_basis[p][t] * NGR1;  //last term of page 8 of my notes 17 Feb 2021
											laplacian_center[el][ij][rs] += 0.5 * tmp2 * A * sps_boundary_basis[p][t] * NGR1;  //last term of page 9 of my notes 17 Feb 2021

											boundary_source[edge_index][ij][(elem_side / 2) * (Knod - 2 * p - 1) + p] += A * NGR1 * sign * NeuMatrix(q, p) *
												std::sqrt(GB[el][d][t][p][d] * face_jac[el][ijp][p]) * (tmp2 - SNGLB1); //NBC terms on page 8,9 of mu notes 17 Feb 2021

											double tmp4 = sps_boundary_grad_basis[p][t] - g_prime[t] * sps_boundary_basis[p][t]; //SNGLB(p,t)-gLp(t)*SBLB(p,t)
											for (int c = 0; c < Knod; ++c) {
												const int j1 = p * d + c * (1 - d);
												const int i1 = p * (1 - d) + c * d;
												const int rs = j1 * Knod + i1; //cumulative index
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
					ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
					dn = ijpm / 2;
					tn = ijpm % 2;

					for (int j = 0; j < Knod; ++j) {
						for (int i = 0; i < Knod; ++i) {
							int ij = j * Knod + i; //local cumulative sps index in element k
							for (int m = 0; m < Knod; ++m) {
								for (int r = 0; r < 2; r++) {
									const int a = m * (1 - r) + i * r;
									const int b = m * r + j * (1 - r);
									const int q = a * d + b * (1 - d);
									const int Mq = ((d + t + dn + tn + 1) % 2) * (Knod - 2 * q - 1) + q; //index of the flux point in the neighbor element
									double SNGLB1 = sps_sps_grad_basis[m][i * (1 - r) + j * r];
									double NGR1 = sps_grad_radau[a * (1 - d) + b * d][t];
									double A = G[el][b][a][r][d]; //(g_r . g_d)/J at k,a,b
									tmp2 = 0.;
									for (int s = 0; s < 2; s++) tmp2 += sps_boundary_basis[m][s] * sps_grad_radau[j * r + i * (1 - r)][s];
									for (int p = 0; p < Knod; ++p) {
										const int j1 = p * dn + Mq * (1 - dn); //j of the neighbor element
										const int i1 = p * (1 - dn) + Mq * dn; //i of the neighbor element
										const int rs = j1 * Knod + i1; //cumulative index in the neighbor element
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
							const int ij = j * Knod + i; //local cumulative sps index in element k
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
								int j1 = m * d + j * (1 - d);
								int i1 = m * (1 - d) + i * d;
								int rs = j1 * Knod + i1;
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

	if (LHS_type==1) Poisson_solver_Eigen_setup(laplacian_center, laplacian_neighbor); //Use Eigen to solve for the linear system
	else if (LHS_type==2) Poisson_solver_Hypre_setup(laplacian_center, laplacian_neighbor);  //use Hypre to solve for the linear system
	else if (LHS_type==3) Poisson_solver_AMGCL_setup(laplacian_center, laplacian_neighbor);//use AMGCL to solve for the linear system

	//  **************** free unnecessary memory ****************
	free_array<double>(laplacian_center);
	free_array<double>(laplacian_neighbor);
}

char HO_2D::Euler_time_integrate(double*** vort_in, double*** vort_out, double coeff) {
	/*uses the Euler explicit time integration to solve the advection-diffusion equation for vorticity
	first it calculates the velocity field (which is needed for advective term). for this purpose, it first solves the psi Poisson quation:
	Laplacian(psi)=vort_in (at time step n+coeff). Then calculates the Cartesian velocity field from psi.
	The calculated velocity field corresponds to the vorticity field of vort_in.
	Then using this velocity field, the advection diffusion equation for the vorticity is solved to calculate vort_out (k2,k3,k4)
	*/

	// for now the BC are setup once in the setup_IC_BC_SRC. This is for a time-dependent BCs
	const double current_time = (ti + coeff) * dt;
	update_BCs(current_time);

	solve_Poisson(vort_in); //calculate the streamfunction corresponding to the vort_in, i.e. Laplacian(psi)=-vort_in. solves for stream_function
	calc_velocity_vector_from_streamfunction(); //calculate the Cartesian velocity field velocity_cart from psi
	update_diffusion_BC();  //sets the BC for the vorticity
	calc_RHS_diffusion(vort_in);  //stores Laplace(w) in RHS_diffusive[el][Knod][Knod]
	update_advection_BC(); //sets the BC for the advection term
	if (advection_type==1) calc_RHS_advection(vort_in); //stores -div(vw) in RHS_advective[el][Knod][Knod]
	else if (advection_type == 2) calc_RHS_advection_continuous(vort_in); //stores -div(vw) in RHS_advective[el][Knod][Knod]

	//update the vorticity field now
	for (int el = 0; el < mesh.N_el; ++el)
		for (int ky = 0; ky < Knod; ++ky)
			for (int kx = 0; kx < Knod; ++kx)
				vort_out[el][ky][kx] = dt * (RHS_advective[el][ky][kx] + RHS_diffusive[el][ky][kx]);

	return 0;
}

char HO_2D::solve_Poisson(double const* const* const* vort_in) {
	//This subroutine solves the Poisson's equation for the stream function field
	// in: vort_in is the right hand side vorticity field
	// out: It solves for the streamfunction and writes the result into the stream_function array
	// The boundary condition for psi is BC_Poisson which is set in update_BCs(time)

	int Ksq = Knod * Knod;

	int N_el = mesh.N_el;
	double* RHS = new double[N_el * Ksq];

	/* these lines were written to just test the poisson solver
	for (int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b)
		for (int m = 0; m < Knod; ++m) {
			//set the proper normal_vel[el_b][m] and BC_parl_vel[el_b][m] to construct proper Neuman BC for sai
			Cmpnts2 coor(0., 0.); //coordinate of sps[j][i]
			for (int n = 0; n < mesh.Lnod; ++n) {
				int node_index = mesh.edges[el_b].nodes[mesh.tensor2FEM(n)];
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
	*/

	for (int el = 0; el < N_el; ++el) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				int ij = j * Knod + i;
				RHS[el * Ksq + ij] = -vort_in[el][j][i] * vol_jac[el][j][i];
				//RHS[el * Ksq + ij] = 0.; //temporary to just test the code for poisson solver performance
			}
		}
	}

	for (int el_b = 0; el_b < mesh.N_edges_boundary; ++el_b) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				int ij = j * Knod + i;
				int el = mesh.boundary_elem_ID[el_b].element_index;
				for (int alpha = 0; alpha < Knod; ++alpha)
					RHS[el * Ksq + ij] -= boundary_source[el_b][ij][alpha] * BC_Poisson[el_b][alpha];
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
		//std::cout << "  iterations:     " << bicg_Eigen.iterations() << std::endl;
		//std::cout << "  estimated error: " << bicg_Eigen.error() << std::endl;
		std::cout << "  iters:  " << bicg_Eigen.iterations() << "  and error:  " << bicg_Eigen.error() << std::endl;


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

	else if (LHS_type == 2) { //Hypre
		// die - there is no Hypre
		assert("Hypre solver is unsupported in HO-CXX, quitting.");

	}  //if LHS_type==2 (Hypre)

	else if (LHS_type == 3) { //AMGCL
		std::vector<double> AMGCL_sol(N_el * Ksq);
		for (int el = 0; el < N_el; ++el)
			for (int j = 0; j < Knod; ++j)
				for (int i = 0; i < Knod; ++i)
					AMGCL_sol[el * Ksq + j * Knod + i] = stream_function[el][j][i];  //initial guess is the solution from previous time step

		//std::fill(RHS_AMGCL.begin(), RHS_AMGCL.end(), 1.);
		std::copy(&RHS[0], &RHS[N_el * Ksq], RHS_AMGCL.begin());
		int iters;
		double error;
		std::tie(iters, error) = (*AMGCL_solver)(RHS_AMGCL, AMGCL_sol);
		printf("#Iters=%d Error=%lf\n", iters, error);


		for (int el = 0; el < N_el; ++el)
			for (int j = 0; j < Knod; ++j)
				for (int i = 0; i < Knod; ++i) {
					int index = el * Ksq + j * Knod + i;
					stream_function[el][j][i] = AMGCL_sol[index];
				}
	}

	delete[] RHS;

	return 0;
}

void HO_2D::update_BCs(const double current_time) {
	/* updates the BC for the slip and normal velocity for the time "time".
	typically velocity BC is fixed over time,
	so this gives more flexibility for velocity BC changes in time */

	if (problem_type==1) {  // for testing: Assume the lid moves to the right from stall and reaches velocity 1 and then again goes back to zero, with a period of 20
		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
			if (mesh.boundaries[G_boundary].name == "top") { // the proper boundary is selected
				for (int edge_index : mesh.boundaries[G_boundary].edges) { //loop over all the edges forming the top boundary
					for (int i = 0; i < Knod; ++i) // loop over all sps on each edge on the proper boundary
						BC_parl_vel[edge_index][i] = -std::fabs(std::sin(M_PI/20.*current_time)); //the lid moves CW on the global boundary
				}
			}
		}
	}

	else if (problem_type==2) {  //Backward Facing step: the inlet velocity goes from 0 to 2 and reverts back
		//BFS: inlet: specific normal vel ---> calculate psi based on the mass flux
		double y_min = 0, y_max = 0;
		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) { //find the min and max y at the inlet
			if (mesh.boundaries[G_boundary].name != "inlet") continue;
			y_min = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
			y_max = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
			for (int edge_index: mesh.boundaries[G_boundary].edges) { //loop over all the edges to find min y and max y of the inlet boundary
				y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_min); // compare with the first node of the edge
				y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_min); // compare with the last  node of the edge
				y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_max); // compare with the first node of the edge
				y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_max); // compare with the last  node of the edge
			}
		}

		const double time_coeff = 2.*std::fabs(std::sin(M_PI/40.*current_time));

		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
			if (mesh.boundaries[G_boundary].name == "inlet") { //vorticity, psi and normal velocity are defined
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					for (int k=0; k<Knod; ++k) {
						double y=0.;
						for (int j=0; j<mesh.Lnod; ++j) y += gps_sps_basis[j][k] * mesh.nodes[ mesh.edges[edge_index].nodes[mesh.tensor2FEM(j)] ].coor.y; //y at k sps on the edge_index boundary edge

						if (false) { // uniform velocity profile into the domain
							BC_Poisson[edge_index][k] = time_coeff*(y - y_min);
							BC_diffusion[edge_index][k] =0.; //zero vorticity
							BC_normal_vel[edge_index][k] = -time_coeff;
						}
						else { //parabolic inlet profile
							BC_Poisson[edge_index][k] = time_coeff * pow(y - y_min,2.) * (3.*y_max-y_min-2.*y)/ pow(y_max-y_min,2.); //psi BC
							BC_diffusion[edge_index][k] = time_coeff * 6.*(2.*y-y_min-y_max)/ pow(y_max-y_min,2.);  //vorticity BC
							BC_normal_vel[edge_index][k] = time_coeff * 6.*(y-y_min)*(y-y_max)/pow(y_max-y_min,2.);
						}
					}
				}
			}  //case inlet


			else if (mesh.boundaries[G_boundary].name == "top") {
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					for (int k=0; k<Knod; ++k) {
						BC_Poisson[edge_index][k] = time_coeff * (y_max - y_min);
					}
				}
			}
		}  //for int Gboundary
	}  //else if problem_type ==2

	else if (problem_type==3) { //flow over an object
		double y_min = 0, y_max = 0;
		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) { //find the min and max y at the inlet
			if (mesh.boundaries[G_boundary].name != "inlet") continue;
			y_min = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
			y_max = mesh.nodes[ mesh.edges[ mesh.boundaries[G_boundary].edges[0] ].nodes[0] ].coor.y;
			for (int edge_index: mesh.boundaries[G_boundary].edges) { //loop over all the edges to find min y and max y of the inlet boundary
				y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_min); // compare with the first node of the edge
				y_min = std::min(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_min); // compare with the last  node of the edge
				y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[0] ].coor.y, y_max); // compare with the first node of the edge
				y_max = std::max(mesh.nodes[ mesh.edges[edge_index].nodes[1] ].coor.y, y_max); // compare with the last  node of the edge
			}
		}

		double time_coeff = 2.*std::fabs(std::sin(M_PI/40.*current_time));

		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
			if (mesh.boundaries[G_boundary].name == "inlet") { //vorticity, psi and normal velocity are defined
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					for (int k=0; k<Knod; ++k) {
						double y=0.;
						for (int j=0; j<mesh.Lnod; ++j) y += gps_sps_basis[j][k] * mesh.nodes[ mesh.edges[edge_index].nodes[mesh.tensor2FEM(j)] ].coor.y; //y at k sps on the edge_index boundary edge
						BC_Poisson[edge_index][k] = time_coeff * (y - y_min); // uniform velocity profile = 1 into the domain
						BC_normal_vel[edge_index][k] = -time_coeff; //In the Fortran code, it is considered as the volume flux but here it is the normal out of domain VELOCITY, to be consistent with the BC_parl_vel
					}
				}
			}  //case inlet

			else if (mesh.boundaries[G_boundary].name == "top") {
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					for (int k=0; k<Knod; ++k) {
						BC_Poisson[edge_index][k] = time_coeff * (y_max - y_min); //the average velocity is 1, so the volume flux in equal to the area of the inlet = psi(top) - psi(bottom)
						BC_diffusion[edge_index][k] = 0.;
						BC_parl_vel[edge_index][k] = -time_coeff; // goes to right, so CW. it is used to calculate the normal derivative of vorticity
					}
				}
			}

			else if (mesh.boundaries[G_boundary].name == "bottom") {
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					for (int k=0; k<Knod; ++k) {
						BC_parl_vel[edge_index][k] = time_coeff; // Goes to the right, so CCW, so positive sign
					}
				}
			}

			else if (mesh.boundaries[G_boundary].name == "object") {
				for (int edge_index : mesh.boundaries[G_boundary].edges) {
					for (int k=0; k<Knod; ++k) {
						BC_Poisson[edge_index][k] = time_coeff * 0.5 * (y_max-y_min); //psi=0
					}
				}
			}
		} //for Gboundary
	} //if problem_type==3

	else if (problem_type==11) { //hybrid run

		// what fraction through the outer time step are we?
		// for the beginning of the first step, weights are fac1=1.0, fac2=0.0
		// for the end of the last step, weights are fac1=0.0, fac2=1.0
		const double fac1 = (current_time-time_start) / (time_end-time_start);
		const double fac2 = 1.0-fac1;
		//std::cout << "  substep normalized time " << fac1 << std::endl;

		// then set BC values for each element
		// note that the hybrid BC_VelNorm_start, etc. arrays are indexed only for the open bdry
		// but the simulation BCs BC_normal_vel, etc. are indexed over all boundry edges!
		unsigned int eoffset = 0;
		for (int G_boundary=0; G_boundary<mesh.N_Gboundary; G_boundary++) {
			auto& bdry = mesh.boundaries[G_boundary];

			if (bdry.name == "open") {

				for (unsigned int i=0; i<bdry.N_edges; ++i) {
					const int eidx = eoffset + i;	// edge index in the master list
					for (int j=0; j<Knod; ++j) {
						// set normal and parallel vels
						BC_normal_vel[eidx][j] = -(fac1*BC_VelNorm_start[i][j] + fac2*BC_VelNorm_end[i][j]);
						BC_parl_vel[eidx][j]   = -(fac1*BC_VelParl_start[i][j] + fac2*BC_VelParl_end[i][j]);
						//if (i>150) std::cout << "set bc vel on " << eidx << " " << i << " " << j << " to " << BC_normal_vel[eidx][j] << " " << BC_parl_vel[eidx][j] << std::endl;
						BC_Poisson[eidx][j]    = BC_parl_vel[eidx][j];

						// and vorticity for diffusion - is this right? BC_vorticity is set in calc_RHS_diffusion:2407
						//BC_vorticity[eidx][j] = fac1*BC_Vort_start[i][j] + fac2*BC_Vort_end[i][j];
						BC_diffusion[eidx][j] = fac1*BC_Vort_start[i][j] + fac2*BC_Vort_end[i][j];
						//if (i>150) std::cout << "set bc vort on " << eidx << " " << i << " " << j << " to " << BC_vorticity[eidx][j] << std::endl;
					}
					// these are already reprojected into eta, xsi components
				}

			} else if (bdry.name == "wall") {

				// apply parallel and normal velocities on the wall boundary
				for (unsigned int i=0; i<bdry.N_edges; ++i) {
					const int eidx = eoffset + i;	// edge index in the master list
					for (int j=0; j<Knod; ++j) {
						BC_normal_vel[eidx][j] = 0.0;
						BC_parl_vel[eidx][j]   = 0.0;
						BC_Poisson[eidx][j]    = 0.0;
						// let calc_RHS_diffusion calculate BC_vorticity
					}
					// do we need to reproject into eta, xi components?
				}
				// no need to set BC_vorticity - the solver takes care of the walls

			} else if (bdry.name == "inlet") {

			} else if (bdry.name == "outlet") {

			}

			// keep track of where in the full list of boundary edges this boundary's edges begin
			eoffset += bdry.N_edges;
		}

		// finally, apply the vorticity times the weight to all of the elements
		for (int i=0; i<mesh.N_el; ++i) {
			for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
				vorticity[i][j][k] = (1.0-Vort_wgt[i][j][k]) * vorticity[i][j][k]
										+ Vort_wgt[i][j][k]*(fac1*Vort_start[i][j][k] + fac2*Vort_end[i][j][k]);
			}
		}

	}  //else if problem_type ==11

	else {  //for any arbitrary problem, specify the BC's of the problem here: modify the below terms as the problem is designed
		std::fill(BC_no_slip, BC_no_slip + mesh.N_Gboundary, true);
		std::fill(BC_switch_Poisson, BC_switch_Poisson + mesh.N_edges_boundary, DirichletBC);
		std::fill(&BC_Poisson[0][0], &BC_Poisson[0][0] + mesh.N_edges_boundary * Knod, 0.);
		std::fill(BC_switch_diffusion, BC_switch_diffusion + mesh.N_edges_boundary, NeumannBC);
		std::fill(&BC_normal_vel[0][0], &BC_normal_vel[0][0] + mesh.N_edges_boundary * Knod, 0.);
		std::fill(&BC_parl_vel[0][0], &BC_parl_vel[0][0] + mesh.N_edges_boundary * Knod, 0.);
	}

	// Now convert the outward normal directions into the local h direction
	for (int el_b=0; el_b<mesh.N_edges_boundary; ++el_b) {
		//int element_index = mesh.boundary_elem_ID[el_b].element_index; //index of the element that has a face on the global boundary
		const int elem_side = mesh.boundary_elem_ID[el_b].side; //side of the element that has the edge on the global boundary
		const int ijp = f2i[elem_side];
		const double sign = 2.*(ijp%2) - 1.;
		for (int k=0; k<Knod; ++k) {
			// set and check normal velocity
			BC_normal_vel[el_b][k] *= sign;
			if (BC_switch_diffusion[el_b]==NeumannBC) BC_diffusion[el_b][k] *= sign;
			// what about parallel/tangential velocity?
		}
	}

}

void HO_2D::update_advection_BC() {
	// ************************ form the BC_advection *************************
	//std::fill(BC_switch_advection, BC_switch_advection + mesh.N_edges_boundary, DirichletBC); //pick Dirichlet BC for advection term based on the Fortran code
	for (int edge_index = 0; edge_index < mesh.N_edges_boundary; ++edge_index) {
		const int element_index = mesh.boundary_elem_ID[edge_index].element_index; //index of the element that has a face on the global boundary
		const int elem_side = mesh.boundary_elem_ID[edge_index].side; //side of the element that has the edge on the global boundary
		const int ijp = f2i[elem_side];
		//int direction = ijp / 2; //0 means xsi, 1 means eta direction

		for (int i = 0; i < Knod; ++i) {
			const int edge_flux_point_index = (elem_side / 2) * (Knod - 2 * i - 1) + i;
			/*double volumetric_flux_boundary[2] = {
							BC_cart_vel[edge_index][edge_flux_point_index].x * face_Dy_Dxsi[element_index][ijp][i].y - BC_cart_vel[edge_index][edge_flux_point_index].y * face_Dx_Dxsi[element_index][ijp][i].y,
							-BC_cart_vel[edge_index][edge_flux_point_index].x * face_Dy_Dxsi[element_index][ijp][i].x + BC_cart_vel[edge_index][edge_flux_point_index].y * face_Dx_Dxsi[element_index][ijp][i].x };*/

			//BC_advection[edge_index][edge_flux_point_index] = BC_vorticity[edge_index][edge_flux_point_index] * volumetric_flux_boundary[direction];
			//
			// NOTE Adrin thinks this following line should not have the face_Anorm multiplier (05-14),
			//      but Mohammad explained in Discord on 05/16
			//
			BC_advection[edge_index][edge_flux_point_index] = BC_vorticity[edge_index][edge_flux_point_index]  * BC_normal_vel[edge_index][edge_flux_point_index] * face_Anorm[element_index][ijp][i];
		}
	}

}

char HO_2D::calc_RHS_advection(double*** vort_in) {
	// ********* This subroutine calculates the - div(VW), V is Cartesian velocity vector, W is vorticity *********

	// *********************************************************************

	//local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	double** local_vort = allocate_array<double>(Knod, 2);
	//local array to hold the contravariant flux on a row of csi [0] (xsi dir), and row of eta [1] direction (eta dir)
	double** local_vel = allocate_array<double>(Knod, 2);

	//quantities per element, per ijp face per sps index
	//interpolated vorticity from a row/column of K nodes on the ijp faces
	double*** bndr_vort = allocate_array<double>(mesh.N_el, 4, Knod);
	//interpolated flux of uw and vw in ijp faces, i.e. w*f^tilda on ijp=0,1 faces and w*g^tilda on ijp=2,3 faces
	double*** bndr_flx = allocate_array<double>(mesh.N_el, 4, Knod);
	//flux values at el element, in idir direction, at row_col row or column at k sps, hence [el][idir][row_col][k]= f^tilda[k] w[k] if idir=0 AND g^tilda[k] w[k] if idir=1
	double*** upwnd_flx = allocate_array<double>(mesh.N_el, 4, Knod);

	//0 for xsi and 1: eta directions, so 0 means f^tilda and 1 means g^tilda
	double**** disc_flx = allocate_array<double>(mesh.N_el, 2, Knod, Knod);

	for (int i = 0; i < mesh.N_el; ++i)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < Knod; k++)
				for (int m = 0; m < Knod; m++) disc_flx[i][j][k][m] = 0.;

	double** bndr_disc_flx = allocate_array<double>(4,Knod);	//4 side of a cells

	// ****************************************************************

	//Extrapolate the unknown, Phi and the Flux to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) { //form flux (U^1w/U^2w on xsi/eta dirs, respectively) on the 4 boundaries (bndr_flx) and all interior sps of all elements (disc_flx)
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (int ijp = 0; ijp < 4; ++ijp) bndr_vort[el][ijp][row_col] = bndr_flx[el][ijp][row_col] = upwnd_flx[el][ijp][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_vort[i][0] = vort_in[el][row_col][i];  // xsi direction
				local_vort[i][1] = vort_in[el][i][row_col];  //eta direction
				local_vel[i][0] =  velocity_cart[el][row_col][i].x * vol_Dy_Dxsi[el][row_col][i].y - velocity_cart[el][row_col][i].y * vol_Dx_Dxsi[el][row_col][i].y;  //U^1/J = volumetric flux  = contravariant xsi FLUX component along xsi direction, based on eq. 2.11 in fotis class notes
				local_vel[i][1] = -velocity_cart[el][i][row_col].x * vol_Dy_Dxsi[el][i][row_col].x + velocity_cart[el][i][row_col].y * vol_Dx_Dxsi[el][i][row_col].x;  //U^2/J = volumetric flux = contravariant eta FLUX component along eta direction,
			}

			int ijp = 0; //ijp=0 (ibnd=idir=0)
			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				//discontinuous flux values (f^tilda, g^tilda) on (row_col,1:Knod) solution points for xsi direction and on(row_col, 1:Knod) solution points for eta direction
				for (int i = 0; i < Knod; ++i) disc_flx[el][idir][row_col][i] = local_vel[i][idir] * local_vort[i][idir];

				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					for (int m = 0; m < Knod; m++) bndr_vort[el][ijp][row_col] += local_vort[m][idir] * sps_boundary_basis[m][ibnd];
					const int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
					if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) {

						int edge_index = mesh.elem_neighbor[el].boundary_index[elem_side];
						int edge_flux_point_index = (elem_side / 2) * (Knod - 2 * row_col - 1) + row_col; // the local index on the boundary edge
						/*
						double volumetric_flux_boundary[2] = {
							BC_u_vel[edge_index][edge_flux_point_index] * face_Dy_Dxsi[el][ijp][row_col].y - BC_v_vel[edge_index][edge_flux_point_index] * face_Dx_Dxsi[el][ijp][row_col].y,
							-BC_u_vel[edge_index][edge_flux_point_index] * face_Dy_Dxsi[el][ijp][row_col].x + BC_v_vel[edge_index][edge_flux_point_index] * face_Dx_Dxsi[el][ijp][row_col].x};
						bndr_flx[el][ijp][row_col] = bndr_vort[el][ijp][row_col] * volumetric_flux_boundary[idir];
							//instead of the u,v vel BC it was BC_values[mesh.elem_neighbor[el].boundary_index[elem_side]][(elem_side / 2) * (Knod - 2 * row_col - 1) + row_col];
						*/
						bndr_flx[el][ijp][row_col] = BC_advection[edge_index][edge_flux_point_index];
					}
					else {
						//value of u* w; ijp = 0: on left boundary at the row_col^th row of sps, ijp = 1 : on right boundary at the row_col^th row of sps, ijp = 2 : on bottom boundary at the row_col^th column of sps, ijp = 3, top at row_col^th column
						for (int m = 0; m < Knod; ++m) bndr_flx[el][ijp][row_col] += local_vel[m][idir] * sps_boundary_basis[m][ibnd]; //volumetric flux at the ibnd boundary in idir direction
						bndr_flx[el][ijp][row_col] *= bndr_vort[el][ijp][row_col];
					}
					ijp++;
				}
			}
		}  //for row_col
	}

	for (int el = 0; el < mesh.N_el; el++) {
		for (int ijp = 0; ijp < 4; ++ijp) {
			const int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
			if (mesh.elem_neighbor[el].is_on_boundary[elem_side])
				for (int j = 0; j < Knod; ++j)
					upwnd_flx[el][ijp][j] = bndr_flx[el][ijp][j];
			else {
				const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
				const int ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];

				const double sign = 2. * ((ijp + ijpm) % 2) - 1.;  //to see if the normal fluxes are in different orientations
				const int revert_coeff = (ijp + ijpm - ijp / 2 - ijpm / 2 + 1) % 2; //equal to (d+t+dn+tn+1)%2 = (ijp-d+ijpm-dn+1)%2
				//bool revert_tangent = (ijp == ijpm) || (ijp + ijpm == 3);
				//bool revert_normal = (ijp == ijpm) || (std::fabs(ijp - ijpm) == 2);
				//double sign = revert_normal? -1. : 1.;
				for (int j = 0; j < Knod; ++j) {
					//kk = (revert_tangent) ? Knod - 1 - j : j;
					const int jn = revert_coeff * (Knod - 2 * j - 1) + j; // the neighbor side of the common face index (either j or Knod-j-1)
					upwnd_flx[el][ijp][j] = 0.5 * (bndr_flx[el][ijp][j] + sign * bndr_flx[eln][ijpm][jn]);
					double vort_diff = bndr_vort[el][ijp][j] - bndr_vort[eln][ijpm][jn];
					if (std::fabs(vort_diff) > 1.e-4) {
						const double bndr_vel = (bndr_flx[el][ijp][j] - sign * bndr_flx[eln][ijpm][jn]) / vort_diff; //a_tilda
						upwnd_flx[el][ijp][j] += 0.5 * sgn[ijp] * std::fabs(bndr_vel) * vort_diff;
					}
				}
			}
		}
	}

	for (int el = 0; el < mesh.N_el; el++) {
		int ijp = 0;
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
	free_array<double>(bndr_disc_flx);
	free_array<double>(disc_flx);
	free_array<double>(bndr_vort);
	free_array<double>(bndr_flx);
	free_array<double>(upwnd_flx);
	free_array<double>(local_vort);
	free_array<double>(local_vel);

	return 0;
}

char HO_2D::calc_RHS_advection_continuous(double*** vort_in) {
	// ********* This subroutine calculates the - div(VW), V is Cartesian velocity vector, W is vorticity *********
	// *********************************************************************

	int ijp, ijpm, jn;
	double** local_vort = new double* [Knod]; //local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	for (int i = 0; i < Knod; ++i) local_vort[i] = new double[2];
	double** local_vel = new double* [Knod]; //local array to hold the contravariant flux on a row of csi [0] (xsi dir), and row of eta [1] direction (eta dir)
	for (int i = 0; i < Knod; ++i) local_vel[i] = new double[2];

	double*** bndr_vort, *** bndr_flx, *** upwnd_flx; //quantities per element, per ijp face per sps index
	bndr_vort = new double** [mesh.N_el]; //interpolated vorticity from a row/column of K nodes on the ijp faces
	bndr_flx = new double** [mesh.N_el]; //interpolated flux of uw and vw in ijp faces, i.e. w*f^tilda on ijp=0,1 faces and w*g^tilda on ijp=2,3 faces
	upwnd_flx = new double** [mesh.N_el]; //flux values at el element, in idir direction, at row_col row or column at k sps, hence [el][idir][row_col][k]= f^tilda[k] w[k] if idir=0 AND g^tilda[k] w[k] if idir=1

	for (int i = 0; i < mesh.N_el; ++i) {
		bndr_vort[i] = new double* [4];
		bndr_flx[i] = new double* [4];
		upwnd_flx[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			bndr_vort[i][j] = new double[Knod];
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

	double**** disc_vol_flux = new double*** [mesh.N_el]; //the volumetric flux in xsi and eta directions in all elements at i,j sps
	double**** cont_vol_flux = new double*** [mesh.N_el]; //the continuous volumetric flux in xsi and eta directions in all elements at i,j sps
	for (int el = 0; el < mesh.N_el; el++) {
		disc_vol_flux[el] = new double** [Knod];
		cont_vol_flux[el] = new double** [Knod];
		for (int j = 0; j < Knod; ++j) {
			disc_vol_flux[el][j] = new double* [Knod];
			cont_vol_flux[el][j] = new double* [Knod];
			for (int i = 0; i < Knod; ++i) {
				disc_vol_flux[el][j][i] = new double[2]; //in the xsi and eta directions
				cont_vol_flux[el][j][i] = new double[2]; //in the xsi and eta directions
			}
		}
	}

	for (int i = 0; i < mesh.N_el; ++i)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < Knod; k++)
				for (int m = 0; m < Knod; m++) disc_flx[i][j][k][m] = 0.;

	double** bndr_disc_flx = new double* [4];  //4 side of a cells
	for (int i = 0; i < 4; ++i) bndr_disc_flx[i] = new double[Knod];

	// Here I am using the continuous volumetric flux field for the advective term
	double*** bndr_vol_flux; //volumetric flux at sps of the L/R (0:1) and S/N (2:3) boundaries of all elements
	double*** comm_vol_flux; //common volumetric flux per element, per face (west, east, south, north) per sps index
	bndr_vol_flux = new double** [mesh.N_el];
	comm_vol_flux = new double** [mesh.N_el];

	for (int i = 0; i < mesh.N_el; ++i) {
		bndr_vol_flux[i] = new double* [4];
		comm_vol_flux[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			bndr_vol_flux[i][j] = new double[Knod];
			comm_vol_flux[i][j] = new double[Knod];
		}
	}

	// ****************************************************************
	// ******* calculate the original discontinuous volumetric flux at the sps of all elements ********
	for (int el = 0; el < mesh.N_el; el++)
		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i) {
				disc_vol_flux[el][j][i][0] = velocity_cart[el][j][i].x * vol_Dy_Dxsi[el][j][i].y - velocity_cart[el][j][i].y * vol_Dx_Dxsi[el][j][i].y;  //U^1/J = volumetric flux  = contravariant xsi FLUX component, based on eq. 2.11 in fotis class notes
				disc_vol_flux[el][j][i][1] = -velocity_cart[el][j][i].x * vol_Dy_Dxsi[el][j][i].x + velocity_cart[el][j][i].y * vol_Dx_Dxsi[el][j][i].x;  //U^2/J = volumetric flux = contravariant eta FLUX component
			}
	// *************************************************************************
	// ******* calculate the disc_vol_flux on the elements faces by extrapolation, store in bndr_vol_flux *******
	for (int el = 0; el < mesh.N_el; el++) {
		ijp = 0; //ijp=0 (ibnd=idir=0)
		for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
			for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
				for (int alpha = 0; alpha < Knod; ++alpha) {  //alpha is either i or j, depending on the direction idir
					bndr_vol_flux[el][ijp][alpha] = 0.;
					for (int m = 0; m < Knod; ++m)
						bndr_vol_flux[el][ijp][alpha] += disc_vol_flux[el][alpha * (1 - idir) + m * idir][alpha * idir + m * (1 - idir)][idir]
						* sps_boundary_basis[m][ibnd]; //value of discontinuous vol flux on the ijp side of element el at the alpha flux point (alpha measured in the xsi or eta direction)
				}
				ijp++;
			}
		}
	}
	// **********************************************************************************************************
	// ********* calculate the common flux at the mesh faces **********
	for (int el = 0; el < mesh.N_el; ++el) {
		int ijp = 0;
		for (int idir = 0; idir < 2; idir++) {
			for (int ibnd = 0; ibnd < 2; ibnd++) {
				int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
				if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) {
					int edge_index = mesh.elem_neighbor[el].boundary_index[elem_side];
					for (int m = 0; m < Knod; ++m) {
						int edge_flux_point_index = (elem_side / 2) * (Knod - 2 * m - 1) + m;
						/*double volumetric_flux_boundary[2] = {
										BC_cart_vel[edge_index][edge_flux_point_index].x * face_Dy_Dxsi[el][ijp][m].y - BC_cart_vel[edge_index][edge_flux_point_index].y * face_Dx_Dxsi[el][ijp][m].y,
										-BC_cart_vel[edge_index][edge_flux_point_index].x * face_Dy_Dxsi[el][ijp][m].x + BC_cart_vel[edge_index][edge_flux_point_index].y * face_Dx_Dxsi[el][ijp][m].x };*/


						comm_vol_flux[el][ijp][m] = BC_normal_vel[edge_index][edge_flux_point_index] * face_Anorm[el][ijp][m]; //volumetric_flux_boundary[idir];
					}
				}
				else {
					const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
					ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
					double sign = 2. * ((ijp + ijpm) % 2) - 1.;  //to see if the normal fluxes are in different orientations
					int revert_coeff = (ijp + ijpm - ijp / 2 - ijpm / 2 + 1) % 2; //equal to (d+t+dn+tn+1)%2 = (ijp-d+ijpm-dn+1)%2
					for (int m = 0; m < Knod; ++m) {
						int mn = revert_coeff * (Knod - 2 * m - 1) + m; // the neighbor side of the common face index (either j or Knod-j-1)
						comm_vol_flux[el][ijp][m] = 0.5 * (bndr_vol_flux[el][ijp][m] + sign * bndr_vol_flux[eln][ijpm][mn]);
					}
				}
				ijp++;
			}
		}
	}
	// ****************************************************************
	// ******* calculate the continuous volumetric flux field *********
	for (int el = 0; el < mesh.N_el; ++el)
		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i) {
				cont_vol_flux[el][j][i][0] = disc_vol_flux[el][j][i][0] + (comm_vol_flux[el][0][j] - bndr_vol_flux[el][0][j]) * sps_radau[i][0]
					+ (comm_vol_flux[el][1][j] - bndr_vol_flux[el][1][j]) * sps_radau[i][1];
				cont_vol_flux[el][j][i][1] = disc_vol_flux[el][j][i][1] + (comm_vol_flux[el][2][i] - bndr_vol_flux[el][2][i]) * sps_radau[j][0]
					+ (comm_vol_flux[el][3][i] - bndr_vol_flux[el][3][i]) * sps_radau[j][1];
			}
	// ****************************************************************
	//Extrapolate the unknown, Phi and the Flux to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) { //form flux (U^1w/U^2w on xsi/eta dirs, respectively) on the 4 boundaries (bndr_flx) and all interior sps of all elements (disc_flx)
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (ijp = 0; ijp < 4; ++ijp) bndr_vort[el][ijp][row_col] = upwnd_flx[el][ijp][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_vort[i][0] = vort_in[el][row_col][i];  // xsi direction
				local_vort[i][1] = vort_in[el][i][row_col];  //eta direction
				local_vel[i][0] = cont_vol_flux[el][row_col][i][0];
				local_vel[i][1] = cont_vol_flux[el][i][row_col][1];
			}

			ijp = 0; //ijp=0 (ibnd=idir=0)
			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				for (int i = 0; i < Knod; ++i) disc_flx[el][idir][row_col][i] = local_vel[i][idir] * local_vort[i][idir];

				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					for (int m = 0; m < Knod; m++) bndr_vort[el][ijp][row_col] += local_vort[m][idir] * sps_boundary_basis[m][ibnd];
					int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
					if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) {
						int edge_index = mesh.elem_neighbor[el].boundary_index[elem_side];
						int edge_flux_point_index = (elem_side / 2) * (Knod - 2 * row_col - 1) + row_col; // the local index on the boundary edge
						bndr_flx[el][ijp][row_col] = BC_advection[edge_index][edge_flux_point_index];
					}
					else {
						bndr_flx[el][ijp][row_col] = comm_vol_flux[el][ijp][row_col] * bndr_vort[el][ijp][row_col];
					}
					ijp++;
				}
			}
		}  //for row_col
	}

	for (int el = 0; el < mesh.N_el; el++) {
		for (ijp = 0; ijp < 4; ++ijp) {
			const int elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
			if (mesh.elem_neighbor[el].is_on_boundary[elem_side])
				for (int j = 0; j < Knod; ++j)
					upwnd_flx[el][ijp][j] = bndr_flx[el][ijp][j];
			else {
				const double t = (double)(ijp%2); //the ibnd of the current ijp
				const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
				ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
				//double sign = 2. * ((ijp + ijpm) % 2) - 1.;  //to see if the normal fluxes are in different orientations
				const int revert_coeff = (ijp + ijpm - ijp / 2 - ijpm / 2 + 1) % 2; //equal to (d+t+dn+tn+1)%2 = (ijp-d+ijpm-dn+1)%2
				for (int j = 0; j < Knod; ++j) {
					jn = revert_coeff * (Knod - 2 * j - 1) + j; // the neighbor side of the common face index (either j or Knod-j-1)
					if (comm_vol_flux[el][ijp][j] > 0.)
						upwnd_flx[el][ijp][j] = comm_vol_flux[el][ijp][j] * (t * bndr_vort[el][ijp][j] + (1. - t) * bndr_vort[eln][ijpm][jn]);
					else
						upwnd_flx[el][ijp][j] = comm_vol_flux[el][ijp][j] * ((1.-t) * bndr_vort[el][ijp][j] + t * bndr_vort[eln][ijpm][jn]);
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
					for (int m = 0; m < Knod; ++m) bndr_disc_flx[ijp][j] -= cont_vol_flux[el][j * (1 - idir) + m * idir][j * idir + m * (1 - idir)][idir]
						* vort_in[el][j * (1 - idir) + m * idir][j * idir + m * (1 - idir)] * sps_boundary_basis[m][ibnd]; // - f(1 or -1))
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
	for (int el = 0; el < mesh.N_el; ++el) {
		for (int ijp = 0; ijp < 4; ++ijp) {
			delete[] bndr_vol_flux[el][ijp];
			delete[] comm_vol_flux[el][ijp];
		}
		delete[] bndr_vol_flux[el];
		delete[] comm_vol_flux[el];
	}
	delete[] bndr_vol_flux;
	delete[] comm_vol_flux;



	for (int el = 0; el < mesh.N_el; ++el) {
		for (int j = 0; j < Knod; ++j) {
			for (int i = 0; i < Knod; ++i) {
				delete[] disc_vol_flux[el][j][i];
				delete[] cont_vol_flux[el][j][i];
			}
			delete[] disc_vol_flux[el][j];
			delete[] cont_vol_flux[el][j];
		}
		delete[] disc_vol_flux[el];
		delete[] cont_vol_flux[el];
	}
	delete[] disc_vol_flux;
	delete[] cont_vol_flux;



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

void HO_2D::update_diffusion_BC() {
	//updates the BC values for vorticity (BC_diffusion) when BC_no_slip is true

	//std::fill(BC_switch_diffusion, BC_switch_diffusion + mesh.N_edges_boundary, NeumannBC); //pick Neumann BC for diffusion term based on the Fortran code

	for (int Gboundary = 0; Gboundary < mesh.N_Gboundary; ++Gboundary) { //usually 4 in case of square
		if (BC_no_slip[Gboundary]) {  //no slip condition on the global boundary element Gboundary, so if no slip wall then
			for (int edge_index : mesh.boundaries[Gboundary].edges) {  //loop over edges on the Gboundary global boundary
				double sign = (f2i[mesh.boundary_elem_ID[edge_index].side] % 2)*2. - 1.; //it gives -1 if the edge_index belongs to west and south of element, other wise +1
				for (int m = 0; m < Knod; ++m) velocity_jump[edge_index][m] += sign * BC_parl_vel[edge_index][m]; //if S and W of the element is located on the global boundary the sign is negative
				for (int m = 0; m < Knod; ++m) BC_diffusion[edge_index][m] = velocity_jump[edge_index][m] / (dt * Reyn_inv);
			}
		}
		else {
			//for (int edge_index : mesh.boundaries[Gboundary].edges)  //loop over edges on the Gboundary global boundary
				//for (int m = 0; m < Knod; ++m) BC_diffusion[edge_index][m] = 0.; //clear this with Adrin when I get the chance (why the normal derivative is zero? for simple boundary the concavity is not zero)
		}
	}
}

char HO_2D::calc_RHS_diffusion(double*** vort_in) {
	//this method calculates the 1/Re * Laplacian of vorticity at all sps of all elements stored in Lap_vorticity[el][ky][kx]

	// ****************** definition, memory allocation and initialization ****************
	//double du_dxsi, du_deta;
	unsigned char ijp, ijpm; //left boundary(xsi): ijp=0, right_boundary (xsi): ijp=1; south boundary (eta): ijp=2; north boundary (eta): ijp = 3
	//local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	double** local_vort = allocate_array<double>(Knod, 2);

	//vorticity and grad of vorticity at sps of the L/R (0:1) and S/N (2:3) boundaries of all elements
	double*** bndr_vort = allocate_array<double>(mesh.N_el, 4, Knod);
	double*** bndr_grad_vort = allocate_array<double>(mesh.N_el, 4, Knod);

	//common vorticity and common grad vorticity per element, per face (west, east, south, north) per sps index
	double*** comm_vort = allocate_array<double>(mesh.N_el, 4, Knod);
	double*** comm_grad_vort = allocate_array<double>(mesh.N_el, 4, Knod);

	//f_tilda is the derivative of the continuous vorticity function at all sps, i.e. d/dxsi(w^C)
	double** f_tilda = allocate_array<double>(Knod, Knod);
	double** f_tilda_B = allocate_array<double>(2, Knod);
	//g_tilda is the derivative of the continuous vorticity function at all sps, i.e. d/deta(w^C)
	double** g_tilda = allocate_array<double>(Knod, Knod);
	double** g_tilda_B = allocate_array<double>(2, Knod);


	//**********************************************************************

	//Extrapolate the unknown vorticity and its derivative to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) { //form vorticity xsi/eta dirs, respectively) on the 4 boundaries (bndr_flx) and all interior sps of all elements (disc_flx)
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (int side = 0; side < 4; ++side) bndr_vort[el][side][row_col] = bndr_grad_vort[el][side][row_col] = comm_vort[el][side][row_col] = comm_grad_vort[el][side][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_vort[i][0] = vort_in[el][row_col][i];  // xsi direction, sps at row i
				local_vort[i][1] = vort_in[el][i][row_col];  //eta direction, sps at column i
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
						calc_boundary_comm_vals_meth2(el, el_B, ijp, ibnd, bndr_vort[el][ijp], bndr_grad_vort[el][ijp],
							comm_vort[el][ijp], comm_grad_vort[el][ijp], BC_switch_diffusion[el_B], BC_diffusion[el_B]);

						for (int m=0; m<Knod; ++m)
							BC_vorticity[el_B][(elem_side / 2) * (Knod - 2 * m - 1) + m] = comm_vort[el][ijp][m]; //get the BC_vorticity from the comm_vort values

						//for (int k = 0; k < Knod; ++k) velocity_jump[el_B][k] = comm_vort[el][ijp][k];
					}
					else { //if the element boundary ijp is internal
						const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
						ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];

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
					const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
					ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
					//comm_grad_vort[el][ijp][Knod] is Average of bndr_grad_vort[el][ijp][Knod] + [eln][ijpm][Knod]; same works for Left, South, and North
					double sign = 2. * ((ijp + ijpm) % 2) - 1.;  //to see if the normal fluxes are in different orientations
					int revert_coeff = (ijp + ijpm - ijp / 2 - ijpm / 2 + 1) % 2; //equal to (d+t+dn+tn+1)%2 = (ijp-d+ijpm-dn+1)%2
					for (int k = 0; k < Knod; ++k) comm_grad_vort[el][ijp][k] = 0.5 * (bndr_grad_vort[el][ijp][k] + sign * bndr_grad_vort[eln][ijpm][revert_coeff * (Knod - 2 * k - 1) + k]);
					/*
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
					*/
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
				for (int dummy = 0; dummy < Knod; ++dummy) dw_dxsi += vort_in[el][ky][dummy] * sps_sps_grad_basis[dummy][kx];
				dw_dxsi += (comm_vort[el][0][ky] - bndr_vort[el][0][ky]) * sps_grad_radau[kx][0];
				dw_dxsi += (comm_vort[el][1][ky] - bndr_vort[el][1][ky]) * sps_grad_radau[kx][1];

				for (int dummy = 0; dummy < Knod; ++dummy) dw_deta += vort_in[el][dummy][kx] * sps_sps_grad_basis[dummy][ky];
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
	free_array<double>(local_vort);
	free_array<double>(bndr_vort);
	free_array<double>(bndr_grad_vort);
	free_array<double>(comm_vort);
	free_array<double>(comm_grad_vort);
	free_array<double>(f_tilda);
	free_array<double>(f_tilda_B);
	free_array<double>(g_tilda);
	free_array<double>(g_tilda_B);

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
	/*
	bool revert = (ijp == ijpm) || (ijp + ijpm == 3);// orientation of the two faces are the same or not.
	if (revert)
		for (int k = 0; k < Knod; ++k) com_vort[k] = 0.50 * (vort[k] + vort_nbr[Knod - k-1]); //The common value of Phi at all flux points on one edge of element(eq . 3.1 in hyun diffusion paper))
	else
		for (int k = 0; k < Knod; ++k) com_vort[k] = 0.50 * (vort[k] + vort_nbr[k]); //The common value of Phi at all flux points on one edge of element(eq . 3.1 in hyun diffusion paper))
	*/
	int revert_coeff = (ijp + ijpm - ijp / 2 - ijpm / 2 + 1) % 2; //equal to (d+t+dn+tn+1)%2 = (ijp-d+ijpm-dn+1)%2
	for (int k = 0; k < Knod; ++k) com_vort[k] = 0.50 * (vort[k] + vort_nbr[revert_coeff*(Knod - 2*k - 1)+k]); //The common value of Phi at all flux points on one edge of element(eq . 3.1 in hyun diffusion paper))

	//calc the corrected values of grad(Phi) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG' = -gRB' (=-Knod^2/2)

	for (int k = 0; k < Knod; ++k) Dvort[k] += (com_vort[k] - vort[k]) * g_prime[ibnd]; //derivative of eq 3.7 to use in 4.3 (common derivative))

	//calc the values of grad(Phi) along the mesh interface;

	for (int k = 0; k < Knod; ++k) {
		cross_Dvort[k] = 0.;
		for (int dummy = 0; dummy < Knod; ++dummy)
			cross_Dvort[k] += com_vort[dummy] * sps_sps_grad_basis[dummy][k]; //derivative of common value(at the element edge) in the other direction at flux point k, i.e.derivative w.r.t eta if ijp=0,1 and derivative w.r.t xsi if ijp=2,3
	}

	//calculate the common flux
	for (int k = 0; k < Knod; ++k)
		Dvort[k] = (Dvort[k] * face_Acoef[el][ijp][k] - cross_Dvort[k] * face_Bcoef[el][ijp][k]) / face_jac[el][ijp][k]; //the common flux normal to the global boundary

	delete[] cross_Dvort;
}

void HO_2D::calc_boundary_comm_vals_meth2(int el, int el_b, int ijp, int ibnd, double* vort, double* Dvort, double* com_vort, double* com_grad_vort, const unsigned char BC_type, const double* BC_vals) {
	/*
	el_b: global boundary index
	vort = bndr_vort[el][ijp][0:Knod-1]: extrapolated vorticity at the ijp side of an element
	Dvort = bndr_grad_vort[el][ijp][0:Knod-1]; in: the extrapolated grad of vorticity on the boundary, out: the corrected grad of vorticity on the boundary
	com_vort = comm_vort[el][ijp][0:Knod-1; out: common vorticity value at the ijp side (which is on global boundary)
	com_grad_vort = comm_grad_vort[el][ijp]; out: the common diffusive flux on the global boundary
	*/
	cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir
	double* cross_Dvort = new double[Knod];   //cross derivatives using common values com_vort
	double** NeuMatrix = new double* [Knod], * NeuRHS = new double[Knod];

	for (int m = 0; m < Knod; ++m) NeuMatrix[m] = new double[Knod];

	if (BC_type == DirichletBC) {  //the vorticity (velocity) is provided on the global boundary as BC (Haji bug: it is assumed here that the BC_values have the same orientation as the local elements edge directions: fix it later)
		for (int m = 0; m < Knod; ++m) com_vort[m] = BC_vals[(elem_side / 2) * (Knod - 2 * m - 1) + m]; //the vorticity at the global boundary, this is taken as the common value

		//Get the corrected values of grad(vort) at the mesh interfaces;  Eq. (7.17) in Huynh w gLB'=gDG' = -gRB' (=-Knod^2/2)
		for (int k = 0; k < Knod; ++k)
			Dvort[k] += (com_vort[k] - vort[k]) * g_prime[ibnd]; //corrected grad vorticity at the ibnd global boundary: L/R sided derivative = Derivative(vort^{L/R})

		//calc the corrected values of grad(Phi) along the mesh interface (parallel to the global boundary);
		for (int k = 0; k < Knod; ++k) {
			cross_Dvort[k] = 0.;
			for (int m = 0; m < Knod; ++m)
				cross_Dvort[k] += com_vort[m] * sps_sps_grad_basis[m][k];  // parallel to the boundary component of grad(vorticity) in the xsi^(1-idir) direction where idir=ijp/2
		}

		//calculate the common flux
		for (int k = 0; k < Knod; ++k)
			com_grad_vort[k] = (Dvort[k] * face_Acoef[el][ijp][k] - cross_Dvort[k] * face_Bcoef[el][ijp][k]) / face_jac[el][ijp][k]; // the common flux (diffusive flux) normal to the global boundary at the flux point k

	}
	else if (BC_type == NeumannBC) { //the normal gradient of vorticity is provided on the global boundary as BC: grad(w).n^hat
		for (int k = 0; k < Knod; ++k)
			com_grad_vort[k] = BC_vals[(elem_side / 2) * (Knod - 2 * k - 1) + k] * face_Anorm[el][ijp][k];  //flux created from the vorticity normal gradient: serves as common flux (haji bug: BC_values and the normal derivative should be fixed)

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
	const int N_el = mesh.N_el;
	const int Ksq = Knod * Knod;
	//const int K4 = Ksq * Ksq;;
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
				const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
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

void HO_2D::Poisson_solver_AMGCL_setup(double*** laplacian_center, double**** laplacian_neighbor) {
	// sets up the LHS of the Poisson equation via the Eigen library
	const int N_el = mesh.N_el;
	const int Ksq = Knod * Knod;
	//const int K4 = Ksq * Ksq;

	std::vector<int>    ptr, col; //col_start = ptr[row]; col_end = ptr[row + 1], col has all the column indices
	std::vector<double> val; //val has all the matrix coeffients

	std::vector<std::vector<int>> columns(N_el*Ksq); // a 2D vecor to store the columns indices that have nonzero values for all solution points
	std::vector<std::vector<double>> values(N_el * Ksq);

	nnz = 0;
	for (int el = 0; el < N_el; ++el)
		for (int ij = 0; ij < Ksq; ++ij) {
			const int sps_index = el * Ksq + ij;
			columns[sps_index].clear();
			values[sps_index].clear();
			for (int rs = 0; rs < Ksq; ++rs)
				if (std::fabs(laplacian_center[el][ij][rs]) > 1.e-14) {
					nnz++;
					columns[sps_index].push_back(el * Ksq + rs);
					values[sps_index].push_back(laplacian_center[el][ij][rs]);
				}
		}

	for (int el = 0; el < N_el; ++el)
		for (int elem_side = south; elem_side <= west; elem_side++)
			if (!(mesh.elem_neighbor[el].is_on_boundary[elem_side])) {
				const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
				for (int ij = 0; ij < Ksq; ++ij) {
					const int sps_index = el * Ksq + ij;
					for (int rs = 0; rs < Ksq; ++rs)
						if (std::fabs(laplacian_neighbor[el][elem_side][ij][rs]) > 1.e-14) {
							nnz++;
							columns[sps_index].push_back(eln * Ksq + rs);
							values[sps_index].push_back(laplacian_neighbor[el][elem_side][ij][rs]);
						}
				}
			}

	// now form the required vectors
	ptr.clear(); ptr.resize(N_el * Ksq+1); ptr.push_back(0);
	col.clear(); col.resize(nnz); // We use 5-point stencil, so the matrix
	val.clear(); val.resize(nnz); // will have at most n2 * 5 nonzero elements.
	for (int el = 0; el < N_el; ++el) {
		for (int ij = 0; ij < Ksq; ++ij) {
			int sps_index = el * Ksq + ij;
			std::copy(std::begin(columns[sps_index]), std::end(columns[sps_index]), std::begin(col) + ptr[sps_index]);
			std::copy(std::begin(values[sps_index]), std::end(values[sps_index]), std::begin(val) + ptr[sps_index]);
			ptr[sps_index + 1] = ptr[sps_index ] + columns[sps_index].size();
		}
	}

	//AMGCL_solver.prm.solver.maxiter = 500;
	//AMGCL_solver.prm.solver.tol = 1.e-6;
	AMGCL_Solver::params prm;
	double tol = 1e-15;
	prm.solver.tol = tol;
	prm.solver.maxiter = 150;

	A_AMGCL = std::make_tuple(N_el*Ksq, ptr, col, val);
	AMGCL_solver = new AMGCL_Solver(A_AMGCL, prm);
	//AMGCL_solver = &AMGCL_solver_tmp;
	RHS_AMGCL.resize(N_el * Ksq);
	//(*AMGCL_solver).system_matrix
	//AMGCL_solver.system_matrix= A_amgcl;

/*
	std::fill(RHS_AMGCL.begin(), RHS_AMGCL.end(), 1.);
	int iters;
	double error;
	std::vector<double> AMGCL_sol(N_el * Ksq, 0.);
	std::tie(iters, error) = (*AMGCL_solver)(RHS_AMGCL, AMGCL_sol);
	printf("#Iters=%d Error=%lf\n", iters, error);
*/
}

void HO_2D::save_output(int n) {
	/*  commented for now, no need for it. It was used to test the code for the diffusion or the Poisson equation separetl.
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
						int node_index = mesh.elements[el].nodes[mesh.tensor2FEM(nx, ny)];
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

				Linf_Norm = std::max(Linf_Norm, fabs(u_exact - stream_function[el][j][i])); //can replace streamfunction with vorticity too
				L1_Norm = L1_Norm + fabs(u_exact - stream_function[el][j][i]); //can replace streamfunction with vorticity too
				L2_Norm = L2_Norm + std::pow(u_exact - stream_function[el][j][i], 2);//can replace streamfunction with vorticity too

				file_handle << el << " \t " << j << " \t " << i << " \t " << coor.x << " \t " << coor.y << " \t " << stream_function[el][j][i] << " \t " << u_exact << " \t " << 100. * (u_exact - stream_function[el][j][i]) / u_exact << std::endl;
			}
		}
	}
	L1_Norm /= ((double)mesh.N_el * Knod * Knod);
	L2_Norm = std::sqrt(L2_Norm / ((double)mesh.N_el * Knod * Knod));
	file_handle_norm << L1_Norm << " \t " << L2_Norm << " \t " << Linf_Norm;
	std::cout << "Done writing to file" << std::endl;
	file_handle.close();
	file_handle_norm.close();
	*/
}

void HO_2D::calc_velocity_vector_from_streamfunction() {
	//calculate u and v from the streamfunction field
	/*
	bndr_psi[el][ijp][row_col]: interpolated value for streamfunction from the sps of element el on the ijp side on the flux point row_col
	bndr_grad_psi[el][ijp][row_col]: interpolated value for streamfunction derivative from the sps of element el on the ijp side on the flux point row_col. the derivative is in xsi direction for ijp=0,1 and in eta direction for ijp=2,3;
			out: corrected derivative of psi (in the xsi^{idir} direction) at the row_col flux point of ijp side of element el
	comm_psi: common psi value at the flux points
	comm_grad_psi: the common flux normal to the global boundary at the flux point k
	out: velocity_cart[el][j][i].x,.y
	*/

	// ****************** definition, memory allocation and initialization ****************
	double u_dpsi_dx = 0., u_dpsi_dy = 1., v_dpsi_dx = -1., v_dpsi_dy = 0.; //coeffients to calculate u, v velocity components multiplying d(psi)/dx and d(psi)/dy (by default for vortical velocity)
	if (!get_curl) {u_dpsi_dy = v_dpsi_dx = 0.;  u_dpsi_dx = v_dpsi_dy = 1.; } //for potential velocity
	//double du_dxsi, du_deta;
	unsigned char ijp, ijpm; //left boundary(xsi): ijp=0, right_boundary (xsi): ijp=1; south boundary (eta): ijp=2; north boundary (eta): ijp = 3
	//local array to hold the vorticity along a row of csi [0], and row of eta [1] direction
	double** local_psi = allocate_array<double>(Knod, 2);

	//streamfunction and grad of streamfunction at sps of the L/R (0:1) and S/N (2:3) boundaries of all elements
	double*** bndr_psi = allocate_array<double>(mesh.N_el, 4, Knod);
	double*** bndr_grad_psi = allocate_array<double>(mesh.N_el, 4, Knod);
	//common streamfunction and common grad streamfunction per element, per face (west, east, south, north) at flux points
	double*** comm_psi = allocate_array<double>(mesh.N_el, 4, Knod);
	double*** comm_grad_psi = allocate_array<double>(mesh.N_el, 4, Knod);

	//f_tilda is the derivative of the continuous streamfunction at all sps, i.e. d/dxsi(psi^C)
	double** f_tilda = allocate_array<double>(Knod, Knod);
	//g_tilda is the derivative of the continuous streamfunction at all sps, i.e. d/deta(psi^C)
	double** g_tilda = allocate_array<double>(Knod, Knod);

	// // f_tilda_B is defined on W, E side, g_tilda_B is defined on S, N side
	double** f_tilda_B = allocate_array<double>(2, Knod);
	double** g_tilda_B = allocate_array<double>(2, Knod);

	//**********************************************************************

	//Extrapolate the unknown vorticity and its derivative to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) {
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (int side = 0; side < 4; ++side) bndr_psi[el][side][row_col] = bndr_grad_psi[el][side][row_col] = comm_psi[el][side][row_col] = comm_grad_psi[el][side][row_col] = 0.;
			for (int i = 0; i < Knod; i++) {
				local_psi[i][0] = stream_function[el][row_col][i];  //xsi direction (eta=constant), sps at row i
				local_psi[i][1] = stream_function[el][i][row_col];  //eta direction, sps at column i
			}

			ijp = 0;

			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					for (int m = 0; m < Knod; m++) bndr_psi[el][ijp][row_col] += local_psi[m][idir] * sps_boundary_basis[m][ibnd]; //stream function on the ijp side at row_col flux point
					for (int m = 0; m < Knod; ++m) bndr_grad_psi[el][ijp][row_col] += local_psi[m][idir] * sps_boundary_grad_basis[m][ibnd]; //derivative of streamfunction in xsi/eta direction (xsi^idir) on the ijp side on row_col flux point: west/east:xsi derivative; north/south: eta derivative
					ijp++;
				}
			}
		}
	}
    double* tmp_psi_deriv = new double [Knod]; // to convert the dpsi/dn into dpsi/dh when Neumann BC exists for psi
	if (HuynhSolver_type == 2) { //Use g = g_DG, so presumably DG method
		for (int el = 0; el < mesh.N_el; el++) {
			ijp = 0;
			for (int idir = 0; idir < 2; ++idir) { //idir=0 (xsi); idir=1: eta direction
				for (int ibnd = 0; ibnd < 2; ++ibnd) { //ibnd=0: left/south); ibnd=1: right/north
					cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to combination of ibnd and idir, gives the SENW side
					if (mesh.elem_neighbor[el].is_on_boundary[elem_side]) {
						int el_B = mesh.elem_neighbor[el].boundary_index[elem_side];

						if (BC_switch_Poisson[el_B] == DirichletBC) { //comGradPsi(1:Knod, ijP, el) / Face_Norm(1:Knod, ijP, el)
                            calc_boundary_comm_vals_meth2(el, el_B, ijp, ibnd, bndr_psi[el][ijp], bndr_grad_psi[el][ijp], comm_psi[el][ijp], comm_grad_psi[el][ijp], BC_switch_Poisson[el_B], BC_Poisson[el_B]);
							for (int k = 0; k < Knod; ++k) velocity_jump[el_B][(elem_side / 2) * (Knod - 2 * k - 1) + k] = comm_grad_psi[el][ijp][k] / face_Anorm[el][ijp][k];// partial(psi)/partial h
                        }
						else if (BC_switch_Poisson[el_B] == NeumannBC) {
                            double sign = 2.*(ijp%2)-1.;
                            for (int k=0; k<Knod; ++k) tmp_psi_deriv[k] = sign * BC_Poisson[el_B][k];
                            calc_boundary_comm_vals_meth2(el, el_B, ijp, ibnd, bndr_psi[el][ijp], bndr_grad_psi[el][ijp], comm_psi[el][ijp], comm_grad_psi[el][ijp], BC_switch_Poisson[el_B], tmp_psi_deriv);
							for (int k = 0; k < Knod; ++k) {
								int boundary_flux_point_index = (elem_side / 2) * (Knod - 2 * k - 1) + k;
								BC_normal_vel[el_B][boundary_flux_point_index] = 0.;
								double sign = 1. - 2. * (ijp / 2);// gives +1 for west and east faces, -1 for south and north faces
								for (int dummy = 0; dummy < Knod; ++dummy)
									BC_normal_vel[el_B][boundary_flux_point_index] += sign*comm_psi[el][ijp][dummy] * sps_sps_grad_basis[dummy][k]; //the normal component of velocity. The normal is in the local h direction (g1, g2) of element
                                BC_normal_vel[el_B][boundary_flux_point_index] /= face_Anorm[el][ijp][k]; //added on 14 May, to fix the bug
							}
                        }
					}

					else { //if the element boundary ijp is internal
						const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
						ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];

						calc_internal_comm_vals_meth2(el, ijp, ijpm, ibnd, bndr_psi[el][ijp], bndr_psi[eln][ijpm], bndr_grad_psi[el][ijp], comm_psi[el][ijp]);
					}
					ijp++;
				}
			}
		}

		for (int el = 0; el < mesh.N_el; el++) {
			for (ijp = 0; ijp < 4; ++ijp) {
				cell_sides elem_side = i2f[ijp]; //the side of current element corresponding to ijp
				if (!(mesh.elem_neighbor[el].is_on_boundary[elem_side])) { //if it is internal local boundary
					const int eln = mesh.elem_neighbor[el].neighbor[elem_side]; //element number of the neighbor
					ijpm = f2i[mesh.elem_neighbor[el].neighbor_common_side[elem_side]];
					//comm_grad_psi[el][ijp][Knod] is Average of bndr_grad_psi[el][ijp][Knod] + [eln][ijpm][Knod]; same works for Left, South, and North
					double sign = 2. * ((ijp + ijpm) % 2) - 1.;  //to see if the normal fluxes are in different orientations
					int revert_coeff = (ijp + ijpm - ijp / 2 - ijpm / 2 + 1) % 2; //equal to (d+t+dn+tn+1)%2 = (ijp-d+ijpm-dn+1)%2
					for (int k = 0; k < Knod; ++k) comm_grad_psi[el][ijp][k] = 0.5 * (bndr_grad_psi[el][ijp][k] + sign * bndr_grad_psi[eln][ijpm][revert_coeff * (Knod - 2 * k - 1) + k]);

				}
			}
		}
	}
	delete[] tmp_psi_deriv;

	double dpsi_dxsi, dpsi_deta;
	for (int el = 0; el < mesh.N_el; ++el) {
		//calculate first derivative of the variable
		for (int ky = 0; ky < Knod; ++ky) {
			for (int kx = 0; kx < Knod; ++kx) {
				//calc grads of unknown psi along xsi(for eta = const) and along eta(for xsi = const); du_dxsi and du_deta are the du / dxsi and du / deta at sps[ky][kx]
				dpsi_dxsi = dpsi_deta = 0.;
				for (int dummy = 0; dummy < Knod; ++dummy) dpsi_dxsi += stream_function[el][ky][dummy] * sps_sps_grad_basis[dummy][kx];
				dpsi_dxsi += (comm_psi[el][0][ky] - bndr_psi[el][0][ky]) * sps_grad_radau[kx][0];
				dpsi_dxsi += (comm_psi[el][1][ky] - bndr_psi[el][1][ky]) * sps_grad_radau[kx][1];

				for (int dummy = 0; dummy < Knod; ++dummy) dpsi_deta += stream_function[el][dummy][kx] * sps_sps_grad_basis[dummy][ky];
				dpsi_deta += (comm_psi[el][2][kx] - bndr_psi[el][2][kx]) * sps_grad_radau[ky][0];
				dpsi_deta += (comm_psi[el][3][kx] - bndr_psi[el][3][kx]) * sps_grad_radau[ky][1];

				// now calculate f^tilda g^tilda as per Huynh's paper, i.e. eq 7.16 and a line below in Hyun diff paper
				f_tilda[ky][kx] = (dpsi_dxsi * (vol_Dx_Dxsi[el][ky][kx].y * vol_Dx_Dxsi[el][ky][kx].y + vol_Dy_Dxsi[el][ky][kx].y * vol_Dy_Dxsi[el][ky][kx].y)
					- dpsi_deta * (vol_Dx_Dxsi[el][ky][kx].x * vol_Dx_Dxsi[el][ky][kx].y + vol_Dy_Dxsi[el][ky][kx].x * vol_Dy_Dxsi[el][ky][kx].y)) / vol_jac[el][ky][kx];

				g_tilda[ky][kx] = (dpsi_deta * (vol_Dx_Dxsi[el][ky][kx].x * vol_Dx_Dxsi[el][ky][kx].x + vol_Dy_Dxsi[el][ky][kx].x * vol_Dy_Dxsi[el][ky][kx].x)
					- dpsi_dxsi * (vol_Dx_Dxsi[el][ky][kx].x * vol_Dx_Dxsi[el][ky][kx].y + vol_Dy_Dxsi[el][ky][kx].x * vol_Dy_Dxsi[el][ky][kx].y)) / vol_jac[el][ky][kx];
			}
		}
		// calculate tailda values at the mesh boundaries
		for (int ibnd = 0; ibnd < 2; ++ibnd) {
			for (int k = 0; k < Knod; k++) {
				f_tilda_B[ibnd][k] = comm_grad_psi[el][ibnd][k];
				for (int dummy = 0; dummy < Knod; ++dummy) f_tilda_B[ibnd][k] -= f_tilda[k][dummy] * sps_boundary_basis[dummy][ibnd]; //common derivative - extrapolated derivative of continuous function at left/right of mesh
				g_tilda_B[ibnd][k] = comm_grad_psi[el][ibnd + 2][k];
				for (int dummy = 0; dummy < Knod; ++dummy) g_tilda_B[ibnd][k] -= g_tilda[dummy][k] * sps_boundary_basis[dummy][ibnd]; //common derivative - extrapolated derivative of continuous function at left/right of mesh
			}
		}

		for (int ky = 0; ky < Knod; ++ky) {
			for (int kx = 0; kx < Knod; ++kx) {
				double tempx = f_tilda[ky][kx] + f_tilda_B[0][ky] * sps_radau[kx][0] + f_tilda_B[1][ky] * sps_radau[kx][1];
				double tempy = g_tilda[ky][kx] + g_tilda_B[0][kx] * sps_radau[ky][0] + g_tilda_B[1][kx] * sps_radau[ky][1];
				double dpsi_dx = (tempx * vol_Dx_Dxsi[el][ky][kx].x + tempy * vol_Dx_Dxsi[el][ky][kx].y) / vol_jac[el][ky][kx];//=partial(psi)/partial(x)
				double dpsi_dy = (tempx * vol_Dy_Dxsi[el][ky][kx].x + tempy * vol_Dy_Dxsi[el][ky][kx].y) / vol_jac[el][ky][kx];//=partial(psi)/partial(y)
				velocity_cart[el][ky][kx].x = u_dpsi_dx * dpsi_dx + u_dpsi_dy * dpsi_dy; //combine both potential and vortical velocity cases in one statement to avoid if statements
				velocity_cart[el][ky][kx].y = v_dpsi_dx * dpsi_dx + v_dpsi_dy * dpsi_dy;
			}
		}
	}

	// ******************* free memory on heap ***************
	free_array<double>(local_psi);

	free_array<double>(bndr_psi);
	free_array<double>(bndr_grad_psi);
	free_array<double>(comm_psi);
	free_array<double>(comm_grad_psi);

	free_array<double>(f_tilda);
	free_array<double>(g_tilda);
	free_array<double>(f_tilda_B);
	free_array<double>(g_tilda_B);
}

void HO_2D::save_output_vtk(int indx) {
// The older version of save data to file which uses the boundary values
/*
	//writes the data in VTK format for vorticity and velocity
	Solution points sps (Gauss Legendre) projected onto ParaView nodes(gps), and extrapolated to bounday, edges and corners
	//

	int Lnod = mesh.Lnod;
	int N_el = mesh.N_el;
	int L2 = Lnod * Lnod; //number of geometric nodes per element
	double** om = new double* [N_el], **uv = new double* [N_el], ** vv = new double* [N_el]; //the vorticity and u, v velocity components
	std::vector<std::vector<bool>> used_edge( N_el, std::vector<bool>(4, false)); //2D vector to keep track of the values found on the edges of internal cells (to make sure they ra ecalculated once)
	std::vector<double> nodes_vorticity(mesh.N_nodes); // to save the vorticity field for each node (save data from elements to nodes), so that we dont have to write the data for all nodes for all elements
	std::vector<Cmpnts2> nodes_velocity(mesh.N_nodes); // to save the velocity vector for each node (save data from elements to nodes), so that we dont have to write the data for all nodes for all elements
	for (int i = 0; i < N_el; ++i) {
		om[i] = new double[L2];
		uv[i] = new double[L2];
		vv[i] = new double[L2];
	}

	std::vector<int> loc_ID, loc_ID_inv; //sweeping the cell indices row by row in the eta direction
	if (Lnod == 4) {
		loc_ID = {0, 4, 5, 1, 11, 12, 13, 6, 10, 15, 14, 7, 3, 9, 8, 2};
		//su2pv_ID = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		loc_ID_inv = {0, 3, 15, 12, 1, 2, 7, 11, 14, 13, 8, 4, 5, 6, 10, 9};
	}
	else if  (Lnod == 3) {
		loc_ID = {0, 4, 1, 7, 8, 5, 3, 6, 2};
		//su2pv_ID = {0, 1, 2, 3, 4, 5, 6, 7, 8};
		loc_ID_inv = {0, 2, 8, 6, 1, 5, 7, 3, 4};
	}
	else if (Lnod == 2) {
		loc_ID = {0, 1, 3, 2};
		//su2pv_ID = {0, 1, 2, 3};
		loc_ID_inv = {0, 1, 3, 2};
	}

	// ************ get the u,v field corresponding to vorticity at current time step *************
	double current_time = (ti + 1) * dt;
	update_BCs(current_time);
	solve_Poisson(vorticity); //calculate the streamfunction corresponding to the vort_in, i.e. Laplacian(psi)=-vort_in. solves for stream_function. vorticity is at ti+1 time step
	calc_velocity_vector_from_streamfunction(); //calculate the Cartesian velocity field velocity_cart from psi corresponding to ti+1 time step
	// *********************************************************************************************
	for (int el = 0; el < N_el; ++el) {
		for (int Ly = 0; Ly < Lnod; ++Ly)  //gps in the eta direction
			for (int Lx = 0; Lx < Lnod; ++Lx) {  //gps in the xsi direction
				int gps_index = loc_ID[Ly * Lnod + Lx];
				om[el][gps_index] = uv[el][gps_index] = vv[el][gps_index] = 0.;
				for (int ky = 0; ky < Knod; ++ky)
					for (int kx = 0; kx < Knod; ++kx) {
						double coeff = sps_gps_basis[ky][Ly] * sps_gps_basis[kx][Lx];
						om[el][gps_index] += coeff * vorticity[el][ky][kx];
						uv[el][gps_index] += coeff * velocity_cart[el][ky][kx].x;
						vv[el][gps_index] += coeff * velocity_cart[el][ky][kx].y;
					}
			}
	}

	calc_RHS_diffusion(vorticity);  //to get the updated values for BC_vorticity from the vorticity[el][j][i] field (corresponding to ti+1 time step)
	// ************** Now for the element faces: 1) the inter-mesh values, take the average of values on the left and right (average the contribution from both elements), 2) the boundary edges, just consider the BC values ***************

	for (int el = 0; el < N_el; ++el) {
		for (int elem_side = 0; elem_side < 4; elem_side++) {  //loop over the SENW sides of every element
			cell_sides ijp = f2i[elem_side]; //convert SENW to ijp index
			int d = ijp / 2;
			int t = ijp % 2;
            if (!(used_edge[el][elem_side]) && mesh.elem_neighbor[el].is_on_boundary[elem_side]==false){  //if elem_side is inter-mesh face and has not been computed before
				int eln = mesh.elem_neighbor[el].neighbor[elem_side];
				int elem_siden = mesh.elem_neighbor[el].neighbor_common_side[elem_side];
				int ijpm = f2i[elem_siden]; // the ijp of the neighbor element
				int dn = ijpm / 2;
				int tn = ijpm % 2;
				used_edge[eln][elem_siden] = true;
				for (int L = 0; L < Lnod; ++L) { //the loop on gps of element face (in the direction of element numering)
					int su2pv = L * ((1 - d) * Lnod + d) + t * (Lnod - 1) * (d * (Lnod - 1) + 1);
					int su2pvn = (((t+tn+d+dn+1)%2)*(Lnod-2*L-1)+L) * ((1 - dn) * Lnod + dn) + tn * (Lnod - 1) * (dn * (Lnod - 1) + 1);
					int gps_index = loc_ID[su2pv];
					int gps_indexn = loc_ID[su2pvn]; //the gps index for the neigbor element eln
					om[el][gps_index] = om[eln][gps_indexn] = 0.5 * (om[el][gps_index] + om[eln][gps_indexn]);
					uv[el][gps_index] = uv[eln][gps_indexn] = 0.5 * (uv[el][gps_index] + uv[eln][gps_indexn]);
					vv[el][gps_index] = vv[eln][gps_indexn] = 0.5 * (vv[el][gps_index] + vv[eln][gps_indexn]);
				}
			}
		}
	}

	// ************* Now assign boundary conditions; First initialize the BC values because it was written over above, then fill in with BC values *************

	for (int el_B = 0; el_B < mesh.N_edges_boundary; ++el_B) {
		int el = mesh.boundary_elem_ID[el_B].element_index;  // index of the element whose face is on the global boundary
		cell_sides elem_side = mesh.boundary_elem_ID[el_B].side; //SENW index of the element face that is located on the global boundary
		cell_sides ijp = f2i[elem_side];
		int d = ijp / 2;
		int t = ijp % 2;
		for (int L = 0; L < Lnod; ++L) { //the loop on each boundary edge (in the direction of element numering, NOT local boundary edge numbering)
			int su2pv = L * ((1 - d) * Lnod + d) + t * (Lnod - 1) * (d * (Lnod - 1) + 1);
			int gps_index = loc_ID[su2pv];
			om[el][gps_index] = 0.;
			uv[el][gps_index] = 0.;
			vv[el][gps_index] = 0.;
			for (int m = 0; m < Knod; ++m) {
				double coeff = sps_gps_basis[m][L];
				om[el][gps_index] += coeff * BC_vorticity[el_B][(elem_side / 2) * (Knod - 2 * m - 1) + m];
				uv[el][gps_index] += coeff * BC_cart_vel[el_B][(elem_side / 2) * (Knod - 2 * m - 1) + m].x;
				vv[el][gps_index] += coeff * BC_cart_vel[el_B][(elem_side / 2) * (Knod - 2 * m - 1) + m].y;
			}
		}
	}

	// *********************************************************************************************************************************************************

	// ********************************************************************************************************************************************************************************************
	// ***************************** write omega, and velocity components into a file ****************************
	std::string file_name = "problem";
	file_name.append(std::to_string(problem_type));
	file_name.append("_K");
	file_name.append(std::to_string(Knod));
	file_name.append("_timestep");
	file_name.append(std::to_string(ti+1));
	file_name.append(".vtk");

	std::cout << "     Writing the results after " << ti+1 << "  timesteps into the file:  " << file_name << std::endl;

	double time = (ti+1) * dt;
	std::ofstream file_handle(file_name);

	file_handle << "# vtk DataFile Version 3.0" << std::endl;
	file_handle << "2D Unstructured Grid of Quads" << std::endl;
	file_handle << "ASCII" << std::endl << std::endl;
	file_handle << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file_handle << "POINTS " << N_el * L2 << " double" << std::endl;

	//for (int node_index = 0; node_index < mesh.N_nodes; node_index++)
	//	file_handle << mesh.nodes[].coor.x << " " << mesh.nodes[node_index].coor.y << " 0.0" << std::endl;

	for (int el = 0; el < N_el; el++)
		for (int node_index = 0; node_index < L2; node_index++) {
			int index = mesh.elements[el].nodes[node_index];
			file_handle << mesh.nodes[index].coor.x << " " << mesh.nodes[index].coor.y<< " 0.0" << std::endl;
		}

	file_handle << "CELLS " << N_el << " " << N_el * (L2 + 1) << std::endl;

	//for (int el = 0; el < N_el; ++el) {
	//	file_handle << L2;
	//	for (int node_index = 0; node_index < Lnod; node_index++)
	//		file_handle << " " << mesh.elements[el].nodes[node_index];
	//	file_handle << std::endl;
	//}

	for (int el = 0; el < N_el; ++el) {
		file_handle << L2;
		for (int node_index = 0; node_index < L2; node_index++)
			file_handle << " " << el * L2 + node_index;
		file_handle << std::endl;
	}

	file_handle << "CELL_TYPES " << N_el << std::endl;
	for (int el = 0; el < N_el; el++) file_handle << 70 << std::endl;

	file_handle << "POINT_DATA " << N_el*L2 << std::endl;
	file_handle << "SCALARS vorticity double 1 " << std::endl;
	file_handle << "LOOKUP_TABLE default " << std::endl;

	for (int el = 0; el < N_el; ++el) {
		for (int gps_index = 0; gps_index < L2; gps_index++)
			file_handle << " " << om[el][gps_index];
		file_handle << std::endl;
	}

	file_handle << "VECTORS velocity double " << std::endl;

	for (int el = 0; el < N_el; ++el) {
		for (int gps_index = 0; gps_index < L2; gps_index++)
			file_handle << " " << uv[el][gps_index] << " " << vv[el][gps_index] << " 0.0";
		file_handle << std::endl;
	}

	std::cout << "Done writing to file" << std::endl;
	file_handle.close();


	for (int i = 0; i < N_el; ++i) {
		delete[] om[i];
		delete[] uv[i];
		delete[] vv[i];
	}
	delete[] om;
	delete[] uv;
	delete[] vv;
*/
}

void HO_2D::save_smooth_vtk(int indx) {

	/*
	writes the data in VTK format for vorticity and velocity with properly averaging the values from all coincident elements
	Solution points sps (Gauss Legendre) projected onto ParaView nodes(gps), and extrapolated to bounday, edges and corners
	*/

	int Lnod = mesh.Lnod;
	int N_el = mesh.N_el;
	int L2 = Lnod * Lnod; //number of geometric nodes per element
	double** om = new double* [N_el]; //the vorticity for elements in Gmsh numbering
	Cmpnts2** vel = new Cmpnts2* [N_el]; //the velocity vector for elements in Gmsh numbering
	for (int i = 0; i < N_el; ++i) {
		om[i] = new double[L2];
		vel[i] = new Cmpnts2[L2];
	}

	int gps_index; // the local index of a node based on Gmsh numbering
	int node_index; //the global node index

	std::vector<double> nodes_vorticity(mesh.N_nodes, 0.); // to save the vorticity field for each node (save data from elements to nodes), so that we dont have to write the data for all nodes for all elements
	std::vector<Cmpnts2> nodes_velocity(mesh.N_nodes); // to save the velocity vector for each node (save data from elements to nodes), so that we dont have to write the data for all nodes for all elements
	std::vector<int> nodes_repetition(mesh.N_nodes, 0.); //keeps the number of contributions added on the nodes_vorticity and nodes_velocity


	std::vector<int> loc_ID, /*su2pv_ID,*/ loc_ID_inv; //sweeping the cell indices row by row in the eta direction
	if (Lnod == 4) {
		loc_ID = { 0, 4, 5, 1, 11, 12, 13, 6, 10, 15, 14, 7, 3, 9, 8, 2 };

	}
	else if (Lnod == 3) {
		loc_ID = { 0, 4, 1, 7, 8, 5, 3, 6, 2 };

	}
	else if (Lnod == 2) {
		loc_ID = { 0, 1, 3, 2 };
	}

	// ************ get the u,v field corresponding to vorticity at current time step *************
	//double current_time = (indx + 1) * dt;
	//update_BCs(current_time);
	//solve_Poisson(vorticity); //calculate the streamfunction corresponding to the vort_in, i.e. Laplacian(psi)=-vort_in. solves for stream_function. vorticity is at ti+1 time step
	//calc_velocity_vector_from_streamfunction(); //calculate the Cartesian velocity field velocity_cart from psi corresponding to ti+1 time step
	// *********************************************************************************************
	for (int el = 0; el < N_el; ++el) {
		for (int Ly = 0; Ly < Lnod; ++Ly)  //gps in the eta direction
			for (int Lx = 0; Lx < Lnod; ++Lx) {  //gps in the xsi direction
				gps_index = loc_ID[Ly * Lnod + Lx];
				om[el][gps_index] = 0.; //the vel array is by default initiated to zero in the cmpnts2 constructor, so no need to do it here
				for (int ky = 0; ky < Knod; ++ky)
					for (int kx = 0; kx < Knod; ++kx) {
						double coeff = sps_gps_basis[ky][Ly] * sps_gps_basis[kx][Lx];
						om[el][gps_index] += coeff * vorticity[el][ky][kx];
						vel[el][gps_index].plus(coeff, velocity_cart[el][ky][kx]);
					}
			}
	}

    // ************* calculate the velocity vector on the global boundary based on normal and tangential components ************
    //Cmpnts2** BC_cart_vel = new Cmpnts2* [mesh.N_edges_boundary];
    //for (int i=0; i<mesh.N_edges_boundary; ++i) BC_cart_vel[i] = new Cmpnts2 [Knod];

    for (int el_b=0; el_b<mesh.N_edges_boundary; ++el_b) { //the normal vel is in the local g direction but parl vel is CCW.
        const int elem_index = mesh.boundary_elem_ID[el_b].element_index;
        const cell_sides SENW_side = mesh.boundary_elem_ID[el_b].side;
        const int ijp = f2i[SENW_side];
        const int dir = ijp/2;
        //int bnd = ijp%2;
        const double parl_vel_sign = 1. - 2. * (SENW_side/2); //the sign of tangent velocity, if it is opposite to g direction
        for (int k=0; k<Knod; ++k) {
            double dx_dxsi = face_Dx_Dxsi[elem_index][ijp][k].x;
            double dx_deta = face_Dx_Dxsi[elem_index][ijp][k].y;
            double dy_dxsi = face_Dy_Dxsi[elem_index][ijp][k].x;
            double dy_deta = face_Dy_Dxsi[elem_index][ijp][k].y;

            double g_1_norm = std::sqrt(dx_dxsi*dx_dxsi + dy_dxsi*dy_dxsi);
            double g_2_norm = std::sqrt(dx_deta*dx_deta + dy_deta*dy_deta);

            if (dir==0) {
                double u_g1 = BC_normal_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k];// the normal velocity
                double u_g_2 = parl_vel_sign * BC_parl_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k];// the tangential in the local g direction
                BC_cart_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k].x = (1./g_2_norm) * (u_g1*dy_deta + u_g_2*dx_deta);
                BC_cart_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k].y = (1./g_2_norm) * (-u_g1*dx_deta + u_g_2*dy_deta);
            }

            else {
                double u_g2 = BC_normal_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k]; // the normal velocity
                double u_g_1 = parl_vel_sign * BC_parl_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k]; // the tangential in the local g direction
                BC_cart_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k].x = (1./g_1_norm) * (-u_g2*dy_dxsi + u_g_1*dx_dxsi);
                BC_cart_vel[el_b][(SENW_side / 2) * (Knod - 2 * k - 1) + k].y = (1./g_1_norm) * (u_g2*dx_dxsi + u_g_1*dy_dxsi);
            }
        }
    }



	//calc_RHS_diffusion(vorticity);  //to get the updated values for BC_vorticity from the vorticity[el][j][i] field (corresponding to ti+1 time step)
	// ************** Now  update the elements adjacent to the global boundary ***************
	for (int el_B = 0; el_B < mesh.N_edges_boundary; ++el_B) {
		int el = mesh.boundary_elem_ID[el_B].element_index;  // index of the element whose face is on the global boundary
		cell_sides elem_side = mesh.boundary_elem_ID[el_B].side; //SENW index of the element face that is located on the global boundary
		cell_sides ijp = f2i[elem_side];
		int d = ijp / 2;
		int t = ijp % 2;
		for (int L = 0; L < Lnod; ++L) { //the loop on each boundary edge (in the direction of element numering, NOT local boundary edge numbering)
			int su2pv = L * ((1 - d) * Lnod + d) + t * (Lnod - 1) * (d * (Lnod - 1) + 1); //index in the row by row format (su2pv index)
			gps_index = loc_ID[su2pv]; //Gmsh-style local index
			om[el][gps_index] = 0.;
			vel[el][gps_index].set_coor(0., 0.); //for now just use the extrapolation field from the u, v field in the element to calculate the velocity vector on the boundary.
			for (int m = 0; m < Knod; ++m) {
				double coeff = sps_gps_basis[m][L];
				om[el][gps_index] += coeff * BC_vorticity[el_B][(elem_side / 2) * (Knod - 2 * m - 1) + m];
				vel[el][gps_index].plus(coeff, BC_cart_vel[el_B][(elem_side / 2) * (Knod - 2 * m - 1) + m]);
			}
		}
	}


	// *****************************************************************************************
	// ******** Now accumulate the contributions from all the elements ********
	for (int el = 0; el < N_el; ++el) {
		for (int L = 0; L < L2; ++L) {
			node_index = mesh.elements[el].nodes[L];
			nodes_velocity[node_index].plus(1., vel[el][L]);
			nodes_vorticity[node_index] += om[el][L];
			nodes_repetition[node_index]++;
		}
	}
	// **************************************************************************
	// ****** Now calculate the average values *******
	for (int node_id = 0; node_id < mesh.N_nodes; node_id++) {
		nodes_vorticity[node_id]  /= nodes_repetition[node_id];
		nodes_velocity[node_id].scale(1. / nodes_repetition[node_id]);
	}
	// ***********************************************
	// ***************************** write omega, and velocity components into a file ****************************
	std::string file_name = "problem";
	file_name.append(std::to_string(problem_type));
	file_name.append("_K");
	file_name.append(std::to_string(Knod));
	file_name.append("_timestep");
	file_name.append(std::to_string(indx));
	file_name.append(".vtk");

	std::cout << "     Writing the results after " << indx << "  timesteps into the file:  " << file_name << std::endl;

	//double time = (indx) * dt;
	std::ofstream file_handle(file_name);

	file_handle << "# vtk DataFile Version 3.0" << std::endl;
	file_handle << "2D Unstructured Grid of Quads" << std::endl;
	file_handle << "ASCII" << std::endl << std::endl;
	file_handle << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file_handle << "POINTS " << mesh.N_nodes << " double" << std::endl;

	for (int node_id = 0; node_id < mesh.N_nodes; node_id++)
		file_handle << mesh.nodes[node_id].coor.x << " " << mesh.nodes[node_id].coor.y << " 0.0" << std::endl; //0.0 is the z component

	file_handle << "CELLS " << N_el << " " << N_el * (L2 + 1) << std::endl;

	for (int el = 0; el < N_el; ++el) {
		file_handle << L2;
		for (int L = 0; L < L2; L++)
			file_handle << " " << mesh.elements[el].nodes[L];
		file_handle << std::endl;
	}

	file_handle << "CELL_TYPES " << N_el << std::endl;
	for (int el = 0; el < N_el; el++) file_handle << 70 << std::endl;

	file_handle << "POINT_DATA " << mesh.N_nodes << std::endl;
	file_handle << "SCALARS vorticity double 1 " << std::endl;
	file_handle << "LOOKUP_TABLE default " << std::endl;

	for (int node_id=0; node_id<mesh.N_nodes; node_id++)
		file_handle << nodes_vorticity[node_id] <<std::endl;

	file_handle << "VECTORS velocity double " << std::endl;

	for (int node_id = 0; node_id < mesh.N_nodes; node_id++)
		file_handle << nodes_velocity[node_id].x << " " << nodes_velocity[node_id].y << " 0.0" << std::endl;

	std::cout << "Done writing to file" << std::endl;
	file_handle.close();


	for (int i = 0; i < N_el; ++i) {
		delete[] om[i];
		delete[] vel[i];
	}
	delete[] om;
	delete[] vel;

	//***********************************************************************
	//******** Now write the values for the sample points **********
	if (!N_sample_points) return;
	file_name = "sampling_points";
	file_name.append("_timestep");
	file_name.append(std::to_string(indx));
	file_name.append("_L");
	file_name.append(std::to_string(Lnod));
	file_name.append("_K");
	file_name.append(std::to_string(Knod));
	file_name.append(".dat");

	file_handle.open(file_name);
	for (int sp = 0; sp < N_sample_points; ++sp) {
		double vort = 0.;
		Cmpnts2 vel;
		int el = sample_points_containing_elements[sp];
		for (int L = 0; L < L2; ++L) {
			node_index = mesh.elements[el].nodes[L];
			vort += gps_shapefunc_on_sample_points[L][sp] * nodes_vorticity[node_index];
			vel.plus(gps_shapefunc_on_sample_points[L][sp],	nodes_velocity[node_index]);
		}
		file_handle << sample_points_coor[sp].x << " " << sample_points_coor[sp].y << " " << vort << " " << vel.x << " " << vel.y << std::endl;
	}

	//for (int i=0; i<mesh.N_edges_boundary; ++i) delete[] BC_cart_vel[i];
    //delete[] BC_cart_vel;
}

void HO_2D::save_vorticity_vtk(int indx) {
	const int N_el = mesh.N_el;
	std::string file_name = "vorticity";
	file_name.append(std::to_string(problem_type));
	file_name.append("_K");
	file_name.append(std::to_string(Knod));
	file_name.append("_timestep");
	file_name.append(std::to_string(indx));
	file_name.append(".vtk");

	std::cout << "     Writing the results after " << indx << "  timesteps into the file:  " << file_name << std::endl;

	//double time = indx * dt;
	std::ofstream file_handle(file_name);

	file_handle << "# vtk DataFile Version 3.0" << std::endl;
	file_handle << "2D Unstructured Grid of Quads" << std::endl;
	file_handle << "ASCII" << std::endl << std::endl;
	file_handle << "DATASET POLYDATA" << std::endl;
	file_handle << "POINTS " << Knod * Knod * N_el << " double" << std::endl;

	for (int el = 0; el < N_el; el++)
		for (int j = 0; j < Knod; j++)
			for (int i = 0; i < Knod; i++) {
				Cmpnts2 coor(0., 0.); //coordinate of sps[j][i]
				for (int ny = 0; ny < mesh.Lnod; ++ny)
					for (int nx = 0; nx < mesh.Lnod; ++nx) {
						int node_index = mesh.elements[el].nodes[mesh.tensor2FEM(nx, ny)];
						coor.plus(gps_sps_basis[ny][j] * gps_sps_basis[nx][i], mesh.nodes[node_index].coor);
					}
				file_handle << coor.x << " " << coor.y << " " << " 0." << std::endl;
			}

	file_handle << "POINT_DATA " << Knod * Knod * N_el << std::endl;
	file_handle << "SCALARS vorticity double 1 " << std::endl;
	file_handle << "LOOKUP_TABLE default " << std::endl;

	for (int el = 0; el < N_el; el++) {
		for (int j = 0; j < Knod; j++)
			for (int i = 0; i < Knod; i++)
				file_handle << vorticity[el][j][i] << " ";
		file_handle << std::endl;
	}

	file_handle << "VECTORS velocity double " << std::endl;
	for (int el = 0; el < N_el; el++)
		for (int j = 0; j < Knod; j++)
			for (int i = 0; i < Knod; i++)
				file_handle << velocity_cart[el][j][i].x << " " << velocity_cart[el][j][i].y << " 0." << std::endl;

	std::cout << "Done writing to file" << std::endl;
	file_handle.close();
}

//
// Methods to support hybrid operation
//

void HO_2D::set_defaults() {
	Knod = 3;
	using_hybrid = false;
	Reyn_inv = 1.0;
	HuynhSolver_type = 2;
	time_integration_type = 2;
	advection_type = 1;
	problem_type = 11;
	num_time_steps = 0;
	dt = 1.0;
	ti = 0;
	time_end = 0.0;
	dump_frequency = 100;
	fast = 0;
	LHS_type = 3;	// 1=Eigen, 3=amgcl
}

void HO_2D::enable_hybrid() {
	using_hybrid = true;
}

void HO_2D::set_elemorder(const int32_t _eorder) {
	Knod = (unsigned int)_eorder;
}

void HO_2D::allocate_hybrid_arrays(const size_t _nel, const size_t _neb, const size_t _knod) {

	// first, deallocate if these point to anything
	clean_up();

	// allocate anew
	BC_VelNorm_start = allocate_array<double>(_neb, _knod);
	BC_VelNorm_end   = allocate_array<double>(_neb, _knod);
	BC_VelParl_start = allocate_array<double>(_neb, _knod);
	BC_VelParl_end   = allocate_array<double>(_neb, _knod);
	BC_Vort_start    = allocate_array<double>(_neb, _knod);
	BC_Vort_end      = allocate_array<double>(_neb, _knod);

	Vort_start = allocate_array<double>(_nel, _knod, _knod);
	Vort_end   = allocate_array<double>(_nel, _knod, _knod);
	Vort_wgt   = allocate_array<double>(_nel, _knod, _knod);
}

// method to accept indices and values and return a boundary
boundary HO_2D::arrays_to_boundary(
		const std::string _name,
		const int32_t _n, const int32_t* _idx,
		const int32_t _nval, const double* _bc) {

	boundary thisbdry = boundary();

	thisbdry.name = _name;
	thisbdry.N_edges = (unsigned int)_n / mesh.Lnod;
	const unsigned int dtype = mesh.edge_type_node_number_inv.at(mesh.Lnod);
	thisbdry.edges.resize(thisbdry.N_edges);

	for (size_t i=0; i<thisbdry.N_edges; ++i) {
		// make a new edge
		mesh.edges.emplace_back(edge());
		edge& thisedge = mesh.edges.back();

		// now fill in the information
		thisedge.edge_type = dtype;
		thisedge.N_nodes = mesh.Lnod;
		thisedge.nodes.resize(mesh.Lnod);
		for (size_t j=0; j<thisedge.N_nodes; ++j) {
			thisedge.nodes[j] = (unsigned int)_idx[mesh.Lnod*i+j];
			//std::cout << _name << " bdry edge " << i << " node " << j << " is " << thisedge.nodes[j] << std::endl;
		}
		
		// add the index to the boundary
		thisbdry.edges[i] = (unsigned int)mesh.edges.size() - 1;
		//std::cout << "  and is index " << thisbdry.edges[i] << std::endl;

		// finally, search for the owning element
		for (int el=0; el<mesh.N_elements; ++el) {
			// check each side of the element
			for (int s = south; s <= west; s++) {
				// match first and second node idx
				const unsigned int n0 = mesh.elements[el].nodes[s];
				const unsigned int n1 = mesh.elements[el].nodes[(s + 1) % 4];
				if (n0 == thisedge.nodes[0] && n1 == thisedge.nodes[1]) {
					mesh.elements[el].edges[s] = thisbdry.edges[i];
					//std::cout << _name << " edge " << thisbdry.edges[i] << " matches element " << el << std::endl;
				}
			}
		}
	}

	// deal with the BCs now
    if (_nval > 0) {
		std::cout << "boundary " << _name << " has BC values" << std::endl;
		// this could be tangential (slip wall) or normal (inlet/outlet/bleed surface)
		std::cout << "  " << _nval << " values over " << thisbdry.N_edges << " elements" << std::endl;

		// how many vels per element?
		const size_t nvalper = _nval / thisbdry.N_edges;
		assert(_nval % thisbdry.N_edges == 0 && "Value array is not an even multiple of element count");
		// if 2, then tangential is first, then normal
		thisbdry.parlvel.resize(thisbdry.N_edges);
		thisbdry.normvel.resize(thisbdry.N_edges);
		std::fill(thisbdry.normvel.begin(), thisbdry.normvel.end(), 0.0);

		// HOW DO WE ACCOUNT FOR HIGHER-ORDER GEOM ELEMS? HOW MANY VALS PER ELEM?
		// note that BC_parl_vel, BC_normal_vel are not yet allocated!
		// save one or two numbers in these - we will adjust other BCs in setup_IC_BC_SRC

		if (nvalper == 1) {
			for (size_t i=0; i<thisbdry.N_edges; ++i) {
				thisbdry.parlvel[i] = _bc[i];
			}
		} else if (nvalper == 2) {
			for (size_t i=0; i<thisbdry.N_edges; ++i) {
				thisbdry.parlvel[i] = _bc[i*2];
				thisbdry.normvel[i] = _bc[i*2+1];
			}
		} else {
			assert(false && "Unexpected number of BCs per element in HO_2D::arrays_to_boundary");
		}
	}

	return thisbdry;
}

// load in all the nodes, elements, wall, and open BCs
void HO_2D::load_mesh_arrays_d(const int32_t _iorder,
		const int32_t _nnodes, const double* _xynodes,
		const int32_t _nelems, const int32_t* _idxelems,
		const int32_t _nwall,  const int32_t* _idxwall,
		const int32_t _nopen,  const int32_t* _idxopen) {

	// set some defaults
	mesh.OCC_format = 0;
	mesh.dx_ratio = 1.0;
	mesh.fac = 0.0;

	// Lnod is in Mesh
	mesh.Lnod_in = (unsigned int)_iorder;
	mesh.Lnod = mesh.Lnod_in + 1;

	//std::cout << "in HO_2D::load_mesh_arrays_d with " << _iorder << " " << _nnodes << " " << _nelems << " " << _nwall << " " << _nopen << std::endl;

	// load in the nodes
	mesh.N_nodes = _nnodes/2;
	mesh.nodes.resize(mesh.N_nodes);
	for (int i=0; i<mesh.N_nodes; ++i) {
		mesh.nodes[i].coor.x = _xynodes[2*i];
		mesh.nodes[i].coor.y = _xynodes[2*i+1];
		mesh.nodes[i].node_type = 9;	// do we care?
	}

	// load in the elements - all should be the same type
	const unsigned int nnodeper = mesh.Lnod*mesh.Lnod;
	mesh.N_elements = _nelems / nnodeper;
	const unsigned int etype = mesh.face_type_node_number_inv.at(nnodeper);
	mesh.elements.resize(mesh.N_elements);
	for (int i=0; i<mesh.N_elements; ++i) {
		mesh.elements[i].element_type = etype;
		mesh.elements[i].N_nodes = nnodeper;
		mesh.elements[i].nodes.resize(nnodeper);
		for (size_t j=0; j<nnodeper; ++j) {
			mesh.elements[i].nodes[j] = (unsigned int)_idxelems[nnodeper*i+j];
		}
		for (size_t j=0; j<4; ++j) {
			mesh.elements[i].edges[j] = 0;	// these must be assigned here!!!
		}
	}

	// load in the boundaries
	mesh.edges.clear();
	mesh.boundaries.clear();

	// first, the wall boundary
	if (_nwall > 0) {
		mesh.boundaries.emplace_back(arrays_to_boundary("wall",_nwall,_idxwall,0,nullptr));
	}

	// then the open boundary
	if (_nopen > 0) {
		mesh.boundaries.emplace_back(arrays_to_boundary("open",_nopen,_idxopen,0,nullptr));
	}

	// read in optional inlets and outlets after this but before processing mesh
}

// load in any normal vels and indices for inlets and outlets
void HO_2D::load_inout_arrays_d(
		const int32_t _ninval, const double* _velin,
		const int32_t _nin, const int32_t* _idxin,
		const int32_t _noutval, const double* _velout,
		const int32_t _nout, const int32_t* _idxout) {

	//std::cout << "in HO_2D::load_inout_arrays_d with " << _ninval << " " << _nin << " " << _noutval << " " << _nout << std::endl;

	if (_nin > 0) {
		mesh.boundaries.emplace_back(arrays_to_boundary("inlet",_nin,_idxin,_ninval,_velin));
	}
	if (_nout > 0) {
		mesh.boundaries.emplace_back(arrays_to_boundary("outlet",_nout,_idxout,_noutval,_velout));
	}

	// testing - stop here
	//exit(1);
	std::cout << "GOT FAR ENOUGH" << std::endl;
}

// all input is done, process the mesh and boundaries
void HO_2D::process_mesh_input() {

	mesh.N_edges_boundary = mesh.edges.size();
	mesh.N_Gboundary = mesh.boundaries.size();

	// complete processing of the mesh
	mesh.process_mesh();

	// allocate the proper space for the calculation
	allocate_arrays();
	allocate_hybrid_arrays(mesh.N_el, mesh.N_edges_boundary, Knod);

	// now that the mesh is loaded, prepare everything else
	setup_sps_gps();
	//read_process_sample_points();
	form_bases();
	form_metrics();
	setup_IC_BC_SRC();
	form_Laplace_operator_matrix();

	// dump some debug info
	if (false) for (int i=0; i<mesh.N_el; i+=2000) {
		std::cout << "Elem " << i << " has vol_jac" << std::endl;
		for (int j=0; j<Knod; ++j) {
			for (int k=0; k<Knod; ++k) {
				std::cout << " " << vol_jac[i][j][k];
			}
			std::cout << std::endl;
		}
		std::cout << "  and face_jac" << std::endl;
		for (int j=0; j<4; ++j) {
			for (int k=0; k<Knod; ++k) {
				std::cout << " " << face_jac[i][j][k];
			}
			std::cout << std::endl;
		}
	}
}


// get data from this Eulerian solver

int32_t HO_2D::getsolnptlen() {
	return 2 * mesh.N_el * Knod * Knod;
}

void HO_2D::getsolnpts_d(const int32_t _ptlen, double* _xypts) {

	double** xloc = allocate_array<double>(mesh.Lnod,mesh.Lnod);
	double** yloc = allocate_array<double>(mesh.Lnod,mesh.Lnod);

	// first get the mesh nodes
	for (int i=0; i<mesh.N_el; ++i) {
		for (int j=0; j<mesh.Lnod; ++j) for (int k=0; k<mesh.Lnod; ++k) {
			// xloc(jx,jy) = xcoord(nodeID(t2f(jx,jy),i))
			const size_t nid = mesh.elements[i].nodes[mesh.tensor2FEM(k,j)];
			xloc[j][k] = mesh.nodes[nid].coor.x;
			yloc[j][k] = mesh.nodes[nid].coor.y;
		}

		// then convert these into the solution nodes
		size_t idx = 2*i*Knod*Knod;
		for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
			double x = 0.0;
			double y = 0.0;

			for (int l=0; l<mesh.Lnod; ++l) for (int m=0; m<mesh.Lnod; ++m) {
				// x = x + GeomNodesLgrangeBasis(jx,kx) * GeomNodesLgrangeBasis(jy,ky) * xloc(jx,jy)
				x += gps_sps_basis[l][j] * gps_sps_basis[m][k] * xloc[l][m];
				y += gps_sps_basis[l][j] * gps_sps_basis[m][k] * yloc[l][m];
			}

			_xypts[idx+0] = x;
			_xypts[idx+1] = y;
			idx += 2;
		}
	}

	free_array<double>(xloc);
	free_array<double>(yloc);
}

void HO_2D::getsolnareas_d(const int32_t _veclen, double* _areas) {
	for (int i=0; i<mesh.N_el; ++i) {
		size_t idx = i*Knod*Knod;
		for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
			// areas(ptnum) = wgt(i) * wgt(j) * Vol_Jac(i,j,el)
			//_areas[idx] = vol_jac[i][j][k];	// wrong?
			_areas[idx] = sps_weight[j] * sps_weight[k] * vol_jac[i][j][k];
			idx++;
		}
	}
}

int32_t HO_2D::getopenptlen() {
	// loop over boundary and look for "open"
	for (auto &bdry : mesh.boundaries) {
		if (bdry.name == "open") {
			return 2 * bdry.N_edges * Knod;
		}
	}
	return 0;
}

void HO_2D::getopenpts_d(const int32_t _nopen, double* _xyopen) {

	for (auto &bdry : mesh.boundaries) if (bdry.name == "open") {
		double* xloc = new double [mesh.Lnod];
		double* yloc = new double [mesh.Lnod];

		// now get the solution points on this boundary
		for (unsigned int i=0; i<bdry.N_edges; ++i) {

			for (int j=0; j<mesh.Lnod; ++j) {
				// xloc1D(jx) = xcoord(boundarynodeID(t2f(jx),i))
				const size_t nid = mesh.edges[bdry.edges[i]].nodes[mesh.tensor2FEM(j)];
				//if (i<1000) std::cout << "open pt " << i << " " << j << " nodes " << bdry.edges[i] << " " << nid << std::endl;
				xloc[j] = mesh.nodes[nid].coor.x;
				yloc[j] = mesh.nodes[nid].coor.y;
			}

			size_t idx = 2*i*Knod;
			for (int j=0; j<Knod; ++j) {
				double x = 0.0;
				double y = 0.0;

				for (int l=0; l<mesh.Lnod; ++l) {
					// x = x + GeomNodesLgrangeBasis(jx,kx) * xloc1D(jx)
					x += gps_sps_basis[l][j] * xloc[l];
					y += gps_sps_basis[l][j] * yloc[l];
				}
				//if (i<1000) std::cout << "open pt " << i << " " << j << " is " << x << " " << y << std::endl;

				_xyopen[idx+0] = x;
				_xyopen[idx+1] = y;
				idx += 2;
			}
		}
		delete[] xloc;
		delete[] yloc;
	}
	return;
}

// send vels and vorts to this solver

// update the x,y velocity at the open boundary solution nodes
void HO_2D::setopenvels_d(const int32_t _veclen, double* _xyvel) {
	std::cout << "  In setopenvels_d, incoming size is " << _veclen << std::endl;

	// move data from end to start
	//BC_VelNorm_start = BC_VelNorm_end
	//BC_VelParl_start = BC_VelParl_end
	for (int i=0; i<mesh.N_edges_boundary; ++i) for (int j=0; j<Knod; ++j) {
		BC_VelNorm_start[i][j] = BC_VelNorm_end[i][j];
	}
	for (int i=0; i<mesh.N_edges_boundary; ++i) for (int j=0; j<Knod; ++j) {
		BC_VelParl_start[i][j] = BC_VelParl_end[i][j];
	}

	// read new BC vels into "end"
	unsigned int eoffset = 0;
	for (auto &bdry : mesh.boundaries) {
	if (bdry.name == "open") {
		assert(2*(int)bdry.N_edges*Knod == _veclen && "setopenvels input vector mismatch");

		for (unsigned int i=0; i<bdry.N_edges; ++i) {

			// see line 1948 for code like this
			// find the enclosing element
			const int edge_index = eoffset + i;	// edge index in the master list
			const int el = mesh.boundary_elem_ID[edge_index].element_index;
			// and the side of that element
			const int elem_side = mesh.boundary_elem_ID[edge_index].side;
			// and others
			const int ijp = f2i[elem_side];
			const int t = ijp%2;  //Left/right (South/north face of element)
			const int d = ijp/2; //the direction of the face (0: xsi side, 1: eta side)
			const double coeff = 2. * t -1.; //convert the local h direction to outward normal direction

			// loop over Knod compute nodes along this edge
			for (int j=0; j<Knod; ++j) {

				// pull x and y vels into temps
				const double xvel = _xyvel[2*(i*Knod+j)+0];
				const double yvel = _xyvel[2*(i*Knod+j)+1];

				// find unit normal vector here
				double nx = 0, ny = 0;
				if (d == 1) {
					nx = -face_Dy_Dxsi[el][ijp][j].x;
					ny = face_Dx_Dxsi[el][ijp][j].x;
				} else {
					nx = face_Dy_Dxsi[el][ijp][j].y;
					ny = -face_Dx_Dxsi[el][ijp][j].y;
				}

				nx *= coeff/face_Anorm[el][ijp][j];
				ny *= coeff/face_Anorm[el][ijp][j];

				// convert into Norm and Parl components
				// edge_flux_point_index
				//const int efpi = (elem_side / 2) * (Knod - 2 * j - 1) + j;
				//if (i==0) std::cout << "open pt " << i << " " << j << " in elem " << el << " side " << elem_side << std::endl;

				// BC_VelNorm_end(l,bel) = (Srf_Dx_iDxsi_j(l,2,2,bel) * xvel(l) - Srf_Dx_iDxsi_j(l,1,2,bel) * yvel(l))
				//BC_VelNorm_end[i][j] = face_Dx_Dxsi[el][ijp][efpi].y * xvel - face_Dx_Dxsi[el][ijp][efpi].x * yvel;
				BC_VelNorm_end[i][j] = nx*xvel + ny*yvel;

				// BC_VelParl_end(l,bel) = -(Srf_Dx_iDxsi_j(l,1,2,bel) * xvel(l) + Srf_Dx_iDxsi_j(l,2,2,bel) * yvel(l)) / Face_Norm(l,k,j)
				//BC_VelParl_end[i][j] = -(face_Dx_Dxsi[el][ijp][efpi].x * xvel + face_Dx_Dxsi[el][ijp][efpi].y * yvel) / face_Anorm[el][ijp][j];
				BC_VelParl_end[i][j] = nx*yvel - ny*xvel;
				//if (i%1 == 0) std::cout << "new edge vel bc at " << i << " " << j << " is " << BC_VelNorm_end[i][j] << " " << BC_VelParl_end[i][j]<< std::endl;
			}
		}
	}
	// keep track of where in the full list of boundary edges this boundary's edges begin
	eoffset += bdry.N_edges;
	}

	// now we have "start" and "end" values for the upcoming time step
	return;
}

// update the vorticity at the open boundary solution nodes
void HO_2D::setopenvort_d(const int32_t _veclen, double* _invort) {
	std::cout << "  In setopenvort_d, incoming size is " << _veclen << std::endl;

	for (auto &bdry : mesh.boundaries) if (bdry.name == "open") {
		assert((int)bdry.N_edges*Knod == _veclen && "setopenvort input vector mismatch");

		for (unsigned int i=0; i<bdry.N_edges; ++i) {
			for (int j=0; j<Knod; ++j) {
				// first - move what's in "end" to "start"
				BC_Vort_start[i][j] = BC_Vort_end[i][j];
				// then - read the new BC vorts into "end"
				BC_Vort_end[i][j] = _invort[Knod*i+j];
				//if (i%1 == 0) std::cout << "edge vort bc new and old at " << i << " " << j << " are " << BC_Vort_start[i][j] << " " << BC_Vort_end[i][j] << std::endl;
			}
		}
	}

	// we're going to need this somewhere, using vorticity[el][eta][xi]

	return;
}

// set/update the vorticity at all solution nodes (useful for initialization)
void HO_2D::setsolnvort_d(const int32_t _veclen, double* _invort) {
	assert(mesh.N_el*Knod*Knod == _veclen && "setsolnvort input vector mismatch");

	// first - move what's in "end" to "start"
	for (int i=0; i<mesh.N_el; ++i) {
		for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
			Vort_start[i][j][k] = Vort_end[i][j][k];
		}
	}

	// then - read the new vorts into "end"
	for (int i=0; i<mesh.N_el; ++i) {
		size_t idx = i*Knod*Knod;
		for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
			Vort_end[i][j][k] = _invort[idx+j*Knod+k];
		}
		//if (i%100 == 0) std::cout << "new vort bc at " << i << " is " << Vort_end[i][0][0] << std::endl;
	}

	// this will eventually affect vorticity[el][eta][xi]

	return;
}

// set the particle-to-grid weights for all solution nodes
//   (should be 0.0 most places, 1.0 at the open boundary, falling off quickly)
void HO_2D::setptogweights_d(const int32_t _veclen, double* _inwgt) {
	assert(mesh.N_el*Knod*Knod == _veclen && "setsolnvort input vector mismatch");

	// replace the weight with the incoming weight
	for (int i=0; i<mesh.N_el; ++i) {
		size_t idx = i*Knod*Knod;
		for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
			Vort_wgt[i][j][k] = _inwgt[idx+j*Knod+k];
		}
		//if (i%100 == 0) std::cout << "vwght at " << i << " is " << Vort_wgt[i][0][0] << std::endl;
	}

	return;
}

// march forward

void HO_2D::solveto_d(const double _outerdt, const int32_t _numstep,
		const int32_t _integtype, const double _reyn) {

	// dt is a global (class) variable
	dt = _outerdt / (double)_numstep;
	std::cout << "In solveto_d, marching forward " << _outerdt << " in " << _numstep << " steps" << std::endl;

	// set other class variables
	time_integration_type = (int)_integtype;
	Reyn_inv = 1.0 / _reyn;

	// need times for start and end of outer time step
	// copy from old end
	time_start = time_end;
	time_end += _outerdt;

	// perform time integration, note ti is global step count
	for (int32_t iter = 0; iter < _numstep; ++iter) {

		// refer to the global time step count
		std::cout << "timestep  " << ti << std::endl;

		// take a time step
		solve_advection_diffusion();

		// increment ti, the global time counter
		++ti;
	}
}

// retrieve results

// caller asks for vorticity on all solution nodes
void HO_2D::getallvorts_d(const int32_t _veclen, double* _outvort) {
	assert(mesh.N_el*Knod*Knod == _veclen && "getallvorts vector length mismatch");

	// copy from master vorticity
	for (int i=0; i<mesh.N_el; ++i) {
		size_t idx = i*Knod*Knod;
		for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
			_outvort[idx+j*Knod+k] = vorticity[i][j][k];
		}
		//if (i%100 == 0) std::cout << "vort at " << i << " is " << vorticity[i][0][0] << std::endl;
	}
}

// caller needs array of weights for a solution element
void HO_2D::get_hoquad_weights_d(const int32_t _veclen, double* _outvec) {
	assert(Knod*Knod == _veclen && "get_hoquad_weights vector length mismatch");

	// fraction of the element's area taken up by this solution node's domain
	for (int j=0; j<Knod; ++j) for (int k=0; k<Knod; ++k) {
		_outvec[j*Knod+k] = 0.25 * sps_weight[j] * sps_weight[k];
	}
}

// 
void HO_2D::trigger_write(const int32_t _indx) {
	save_vorticity_vtk((int)_indx);
	save_smooth_vtk((int)_indx);
}

// close, finish

void HO_2D::clean_up() {

	// call the existing memory cleanup routine
	//release_memory();

	// and then clean up hybrid arrays
	free_array<double>(BC_VelNorm_start);
	free_array<double>(BC_VelNorm_end);
	free_array<double>(BC_VelParl_start);
	free_array<double>(BC_VelParl_end);
	free_array<double>(BC_Vort_start);
	free_array<double>(BC_Vort_end);

	free_array<double>(Vort_start);
	free_array<double>(Vort_end);
	free_array<double>(Vort_wgt);
}

