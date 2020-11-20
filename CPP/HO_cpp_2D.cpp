#include "calculation.h"
#include "misc.hpp"
#include <cmath>

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

	sps_radau = new LR_boundary[Knod];
	sps_grad_radau = new LR_boundary[Knod];

	vol_Dx_Dxsi = new Cmpnts2** [N_el];
	vol_Dy_Dxsi = new Cmpnts2** [N_el];
	for (int i = 0; i < N_el; ++i) {
		vol_Dx_Dxsi[i] = new Cmpnts2* [Knod];
		vol_Dy_Dxsi[i] = new Cmpnts2* [Knod];
		for (int j = 0; j < Knod; ++j) {
			vol_Dx_Dxsi[i][j] = new Cmpnts2 [Knod];
			vol_Dy_Dxsi[i][j] = new Cmpnts2 [Knod];
		}
	}

	face_Acoef = new double** [N_el];
	face_Bcoef = new double** [N_el];
	face_jac = new double** [N_el];
	face_norm = new double** [N_el];
	for (int i = 0; i < N_el; ++i) {
		face_Acoef [i] = new double* [4];
		face_Bcoef [i] = new double* [4];
		face_jac[i] = new double* [4];
		face_norm[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			face_Acoef[i][j] = new double [Knod];
			face_Bcoef[i][j] = new double[Knod];
			face_jac[i][j] = new double[Knod];
			face_norm[i][j] = new double[Knod];
		}

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
			gps_boundary_basis[k].right *= 1. - gps_local_coor[j];
			gps_boundary_basis[k].left *= -1. - gps_local_coor[j];
			for (int m = 0; m < Knod; ++m) gps_sps_basis[k][m] *= sps_local_coor[m] - gps_local_coor[j];
			grad_numerator.left = grad_numerator.right = 1.;
			for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] = 1.;

			for (int i = 0; i < Knod; ++i) {
				if (i == k || i == j) continue;
				grad_numerator.left *= -1.0 - gps_local_coor[i];
				grad_numerator.right *= 1.0 - gps_local_coor[i];
				for (int m = 0; m < Knod; ++m) sps_grad_numerator[m] *= sps_local_coor[m] - gps_local_coor[i];
			}
			gps_boundary_grad_basis[k].left += grad_numerator.left;
			gps_boundary_grad_basis[k].right += grad_numerator.right;

			for (int m = 0; m < Knod; ++m) gps_sps_grad_basis[k][m] += sps_grad_numerator[m];
		}

		gps_boundary_basis[k].left /= denominator;
		gps_boundary_basis[k].right /= denominator;

		gps_boundary_grad_basis[k].left /= denominator;
		gps_boundary_grad_basis[k].right /= denominator;

		for (int m = 0; m < Knod; ++m) gps_sps_grad_basis[k][m] /= denominator;
		for (int m = 0; m < Knod; ++m) gps_sps_basis[k][m] /= denominator;
	}

	// *************************************** In this sub section ***************************************
	// calculates and stores the value of right(sps_radau[k].right) and left(sps_radau[k].left) radau functions on the Knod sps and stores in sps_radau[k].left,right
	// calculates and stores the derivative of right(sps_grad_radau[k].right) and left(sps_grad_radau[k].left) radau functions on the Knod sps and stores in sps_grad_radau[k].left,right

	double coef, coefd;

		//		R_k(x) = (-1) ^ k / 2 * (P(x)_k - P(x)_(k - 1)), where P(x)_k is Legendre polynomial of order k, 3.15 in Hyun diffusion paper
		//		derivative :D[R_k(x), x] = (k / 2) * (-1) ^ k * (P(x)_k + P(x)_(k - 1)) / (x + 1)
	if (!Knod) {
		std::cout << "ERROR: Order 0 Radau is Undefined!" << std::endl;
		exit(1);
	}
	else if (Knod == 1) {
		sps_radau[0].left = sps_radau[0].right = 0.5;
		sps_grad_radau[0].right = -0.5;
		sps_grad_radau[0].left = 0.5;
	}
	else {
		coef =  0.5 * minus_one_to_power(Knod);
		coefd = 0.5 * Knod * minus_one_to_power(Knod);
		for (int i = 0; i < Knod; ++i) {
			//	remember right radau(xsi) = left radau(-xsi) : R_k(x) | Right = R_k(-x) | Left
			sps_radau[i].right = coef * (Legendre(Knod, sps_local_coor[i]) - Legendre(Knod-1, sps_local_coor[i])); //value of right Radau function(g_LB) on the sps
			sps_radau[i].left  = coef * (Legendre(Knod, -sps_local_coor[i]) - Legendre(Knod - 1, -sps_local_coor[i])); //value of left Radau function(g_LB) on the sps; 

			// D[R_k(x), x] | Right = -D[R_k(-x), x] | Left
			sps_grad_radau[i].right = coefd * (Legendre(Knod, sps_local_coor[i]) + Legendre(Knod - 1, sps_local_coor[i])) / (1. + sps_local_coor[i]);
			sps_grad_radau[i].left = -coefd * (Legendre(Knod, -sps_local_coor[i]) + Legendre(Knod - 1, -sps_local_coor[i])) / (1. - sps_local_coor[i]);
		}		
	}
	delete[] sps_grad_numerator;
}

void HO_2D::setup_IC_BC_SRC() {
	//setup initial condition, boundary condition and the source/sink terms
	for (int el = 0; el < mesh.N_el; ++el)
		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i)
				initial_vorticity[el][j][i] = 0.;

}

void HO_2D::form_metrics() {
	/* forms
	vol_Dx_Dxsi; //the derivative of x wrt to xsi_s(s is 0, 1 to make dx / dxsi, dx / deta) at the sps(jx, jy) in tensor form on element el.This forms Vol_Dx_iDxsi_j(el, jy, jx).x, .y
	vol_Dy_Dxsi; //the derivative of y wrt to xsi_s(s is 0, 1 to make dy / dxsi, dy / deta) at the sps(jx, jy) in tensor form on element el.This forms Vol_Dy_iDxsi_j(el, jy, jx).x, .y		
	
	!Haji : Vol_Jac(jx, jy, el) : is the cross product of vectors g_1 = (dx / dcsi, dy / dcsi) x g_2 = (dx / deta, dy / deta) at the sps(jx, jy) of element el.g_1 and g_2 are based on Fotis notes(covariant bases).So, the cross product is G = det[dx / dcsi, dy / dxsi; dx / deta, dy / deta] = ratio of volume of original element / volume of transformed element(del_csi(2) * del_eta(2) = 4)
		!Haji : Face_Jac(i(1:Knod), r(0:3), el(1:Nel)) is the G = dx / dxsi * dy / deta - dx / deta * dy / dxsi on the boundary face at the r = 0(left) or r = 1 (right) or r = 2(south) or r = 3(north)face on the i sps of the element el
		!Haji : Face_Acoef(i(1:Knod), r(0:3), el(1:Nel)) is the Case 1) squared length of g_2 = (dx / deta, dy / deta), i.e.g_2(dot)g_2, on the leftand right boundaries(r = 0, 1); Case 2) squared length of g_1 = (dx / dcsi, dy / dcsi), i.e.g_1(dot)g_1, on the bottomand top boundaries(r = 2, 3)
		!Haji: Face_Bcoef(i(1:Knod), r(0:3), el(1:Nel)) is the dot product of g_1 with g_2 = g_1(dot)g_2 = dx / dxsi * dx / deta + dy / dxsi * dy / deta) on the 4 boundaries(left, right, bottom, top(r = 0, 1, 2, 3)
			!Haji: Face_Norm(i(1:Knod), r(0:3), el(1:Nel)) is the norm of Face_Acoef, i.e.Face_Norm(jy, i, el) = Sqrt(Face_Acoef(jy, i, el))
			!Haji : The Face_Acoef and Face_Bcoef, forms the covariant metric tensor g_ij = [g11 g12; g12, g22]] (symmetric matrix))
	*/
	unsigned int Lnod = mesh.Lnod;
	double dx_dxsi, dy_dxsi, dx_deta, dy_deta;
	Cmpnts2** local_coor = new Cmpnts2 * [Lnod]; //local array tp store the coor of the sps in an element, dont really need it just for convenience
	for (int i = 0; i < Lnod; ++i) local_coor[i] = new Cmpnts2 [Lnod];

	for (int el = 0; el < mesh.N_el; ++el) {
		for (int j = 0; j < Lnod; ++j)
			for (int i = 0; i < Lnod; ++i)
				local_coor[j][i] = mesh.nodes[mesh.node_ID[el][tensor2FEM(i, j)]].coor;

		for (int j = 0; j < Knod; ++j)
			for (int i = 0; i < Knod; ++i)
				for (int n = 0; n < Lnod; ++n)
					for (int m = 0; m < Lnod; ++m) {
						vol_Dx_Dxsi[el][j][i].x += gps_sps_grad_basis[m][i] * gps_sps_basis[n][j] * local_coor[n][m].x;	//grad at x - dir  * no grad at y - dir  * x / y - coord of geom
						vol_Dy_Dxsi[el][j][i].x += gps_sps_grad_basis[m][i] * gps_sps_basis[n][j] * local_coor[n][m].y;	//grad at x - dir  * no grad at y - dir  * x / y - coord of geom

						vol_Dx_Dxsi[el][j][i].y += gps_sps_grad_basis[n][j] * gps_sps_basis[m][i] * local_coor[n][m].x;	//grad at y - dir  * no grad at x - dir  * x / y - coord of geom
						vol_Dy_Dxsi[el][j][i].y += gps_sps_grad_basis[n][j] * gps_sps_basis[m][i] * local_coor[n][m].y;	//grad at y - dir  * no grad at x - dir  * x / y - coord of geom
					}

		for (int k = 0; k < Knod; ++k) { //loop on all sps on the faces to calculate metrics
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			// ****** sps on the left (west) boundary (xsi=-1) *********
			for (int j=0; j<Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) { 
					dx_dxsi += gps_boundary_grad_basis[i].left * gps_sps_basis[j][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_grad_basis[i].left * gps_sps_basis[j][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_basis[i].left * gps_sps_grad_basis[j][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_basis[i].left * gps_sps_grad_basis[j][k] * local_coor[j][i].y;
				}
			// note: On lines of constant xsi, it is proven in page 6 of my curvilinear notes that dS = g_2 d_eta, so ||dS|| = ||g_2|| * d_eta, dS is element length
			face_jac[el][3][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][3][k] = dx_deta * dx_deta + dy_deta * dy_deta; //g_2 dot g_2 in fotis notes = g_{22}
			face_Bcoef[el][3][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_norm[el][3][el] = sqrt(face_Acoef[el][3][k]); //||g_2||
			// **********************************************************
			// ****** sps on the right (east) boundary (xsi=+1) *********
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			for (int j = 0; j < Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) {
					dx_dxsi += gps_boundary_grad_basis[i].right * gps_sps_basis[j][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_grad_basis[i].right * gps_sps_basis[j][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_basis[i].right * gps_sps_grad_basis[j][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_basis[i].right * gps_sps_grad_basis[j][k] * local_coor[j][i].y;
				}
			face_jac[el][1][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][1][k] = dx_deta * dx_deta + dy_deta * dy_deta; //g_2 dot g_2 in fotis notes = g_{22}
			face_Bcoef[el][1][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_norm[el][1][el] = sqrt(face_Acoef[el][1][k]); //||g_2||
			// ************************************************************
			// ****** sps on the bottom (south) boundary (eta=-1) *********
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			for (int j = 0; j < Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) {
					dx_dxsi += gps_boundary_basis[j].left * gps_sps_grad_basis[i][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_basis[j].left * gps_sps_grad_basis[i][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_grad_basis[j].left * gps_sps_basis[i][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_grad_basis[j].left * gps_sps_basis[i][k] * local_coor[j][i].y;
				}
			face_jac[el][0][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][0][k] = dx_dxsi * dx_dxsi + dy_dxsi * dy_dxsi; //g_1 dot g_1 in fotis notes = g_{11}
			face_Bcoef[el][0][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_norm[el][0][el] = sqrt(face_Acoef[el][0][k]); //||g_1||
			// ************************************************************
			// ****** sps on the top (north) boundary (eta=+1) *********
			dx_dxsi = dx_deta = dy_dxsi = dy_deta = 0.;
			for (int j = 0; j < Lnod; ++j)
				for (int i = 0; i < Lnod; ++i) {
					dx_dxsi += gps_boundary_basis[j].right * gps_sps_grad_basis[i][k] * local_coor[j][i].x;
					dy_dxsi += gps_boundary_basis[j].right * gps_sps_grad_basis[i][k] * local_coor[j][i].y;

					dx_deta += gps_boundary_grad_basis[j].right * gps_sps_basis[i][k] * local_coor[j][i].x;
					dy_deta += gps_boundary_grad_basis[j].right * gps_sps_basis[i][k] * local_coor[j][i].y;
				}
			face_jac[el][2][k] = dx_dxsi * dy_deta - dx_deta * dy_dxsi;  //local area: G = 1/J in fotis notes
			face_Acoef[el][2][k] = dx_dxsi * dx_dxsi + dy_dxsi * dy_dxsi; //g_1 dot g_1 in fotis notes = g_{11}
			face_Bcoef[el][2][k] = dx_dxsi * dx_deta + dy_dxsi * dy_deta; //g_1 dot g_2 in fotis notes = g_{12}
			face_norm[el][2][el] = sqrt(face_Acoef[el][2][k]); //||g_1||
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


}

char HO_2D::solve_advection_diffusion() {
	//solves the advection diffusion eq. for the vorticity field
	if (time_integration_type == 1) { // first order Euler time integration method
		calc_advection_flux();
		Euler_time_integrate(dt, Reyn_inv, k1, HuynhSolver_type);
			Vort = Vort + k1;
	}




}

char HO_2D::solve_Poisson() {
	//solves the Poisson's equation for the stream function field

}

char HO_2D::calc_advection_flux() {
	// ********* This subroutine calculates the - div(VW), V is Cartesian velocity vector, W is vorticity *********
	
	Cmpnts2* local_vort = new Cmpnts2[Knod]; //local array to hold the vorticity along a row of csi (.x), and row of eta (.y) direction
	Cmpnts2* local_vel  = new Cmpnts2[Knod]; //local array to hold the velocity along a row of csi (.x), and row of eta (.y) direction
	LR_boundary **bndr_vort, **bndr_flx, **upwnd_flx;
	bndr_vort = new LR_boundary* [mesh.N_el];
	bndr_flx = new LR_boundary * [mesh.N_el];
	upwnd_flx = new LR_boundary * [mesh.N_el];

	for (int i = 0; i < mesh.N_el; ++i) {
		bndr_vort[i] = new LR_boundary [Knod];
		bndr_flx[i] = new LR_boundary [Knod];
		upwnd_flx[i] = new LR_boundary [Knod];
		for (int j = 0; j < Knod; ++j)
			for (int bnd_side=0; bnd_side<1; ++bnd_side)
				bndr_vort[i][j][bnd_side] = bndr_vort[i][j][bnd_side] = bndr_flx[i][j][bnd_side] = bndr_flx[i][j][bnd_side] = upwnd_flx[i][j].left = upwnd_flx[i][j].right = 0.;

	}




	//Extrapolate the unknown, Phi and the Flux to the mesh boundaries using Lagrange polynomials of order Knod - 1
	for (int el = 0; el < mesh.N_el; el++) {
		for (int row_col = 0; row_col < Knod; ++row_col) {
			for (int i = 0; i < Knod; i++) {
				local_vort[i].x = vorticity[el][row_col][i];
				local_vort[i].y = vorticity[el][i][row_col];
				local_vel [i].x =  velocity_cart[el][row_col][i].x * vol_Dy_Dxsi[el][row_col][i].y - velocity_cart[el][row_col][i].y * vol_Dx_Dxsi[el][row_col][i].y;  //contravariant xsi FLUX component along xsi direction, based on eq. 2.11 in fotis class notes
				local_vel [i].y = -velocity_cart[el][i][row_col].x * vol_Dy_Dxsi[el][i][row_col].x + velocity_cart[el][i][row_col].y * vol_Dx_Dxsi[el][i][row_col].x;  //contravariant eta FLUX component along eta direction,
			}
		}

		// *************** XSI direction fluxes ***************
		// ****************************************************
		for (int row_col = 0; row_col < Knod; ++row_col) {
			bndr_vort[el][row_col][bnd_side] += local_vort[row_col].x * sps_boundary_basis[row_col].left;
			bndr_vort[el][row_col].right += local_vort[row_col].x * sps_boundary_basis[row_col].right;
		}

		





	}

	

/*
		ijP = 0
		!Extrapolation operations in x(= 1) and y(= 2) directions
		DO idir = 1, 2

		!Extraploated boundary values of Phi and Velocity * Phi to left / south(= 0) and right / north(= 1)
		DO ibnd = 0, 1
		bndrPhi(j, ijP, el) = dot_product(loclPhi(1:Knod, idir), SolnBndryLgrangeBasis(1:Knod, ibnd))

		!mesh to the left of left face or right of right face
		eln = elemID(i2f(ijP), el)
		IF(eln.lt. 0) THEN
		!IF(BC_Switch(-eln).eq.DirichletBC) THEN
		bndrFlx(j, ijP, el) = bndrPhi(j, ijP, el) * BC_Values(j, -eln)
		!ELSEIF(BC_Switch(-eln).eq.NeumannBC) THEN
		!ENDIF
		ELSE
		!Haji: value of u * w; ijp = 0: on left boundary at the j_th row of sps, ijp = 1 : on right boundary at the j_th row of sps, ijp = 2 : on bottom boundary at the j_th column of sps, ijp = 3, top at j_th col
		bndrFlx(j, ijP, el) = bndrPhi(j, ijP, el) * dot_product(loclVel(1:Knod, idir), SolnBndryLgrangeBasis(1:Knod, ibnd))
		ENDIF
		ijP = ijP + 1
		ENDDO

		!Discontinuous flux values on(1:Knod, j) solution points for horizontal direction
		!and on(j, 1:Knod) solution points for vertical direction
		DO i = 1, Knod
		discFlx(i, j, idir, el) = loclVel(i, idir) * loclPhi(i, idir)
		ENDDO

		ENDDO

		ENDDO

		ENDDO
*/





		delete[] local_vort;
		delete[] local_vel;
}

