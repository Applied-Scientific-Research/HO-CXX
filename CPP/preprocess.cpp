#include "preprocess.h"

char Mesh::read_msh_file() {
	// reads a msh file output from the Gmsh software. The msh file is in ASCII 4.1 version of the Gmsh output
	const char* filename = input_msh_file;
	std::cout << "     Gmsh file   ***** " << filename << " *****   opened for reading ..." << std::endl << std::endl;
	int retval = 1, tmp, tmp1, tmp2;
	unsigned int index;
	unsigned int nodes_min_index, nodes_max_index, raw_N_nodes, tag_N_nodes, nodes_total_entities, group_tag, entity_dim, unorganized_node_index = 0;
	unsigned int elements_min_index, elements_max_index, tag_N_elements, elements_total_entities, element_type;
	unsigned int N_boundary, node_index;
	unsigned int entities_N_points, entities_N_curves, entities_N_faces;
	bool found = false;
	double double_field, coorX, coorY, coorZ;
	std::string temp_string;
	std::vector<unsigned int>::iterator its;
	bool check_start;
	edge _edge;
	element2d _face;
	std::vector<Cmpnts2> unorganized_nodes_coor; //store all the raw read coordinartes of all nodes
	std::vector<unsigned int> unorganized_node_mapping; // maps the node numbers read from the msh file into the unorganized indices (helps to remove the indices gap)
	std::vector<char> unorganized_node_type; //0 for corner, 1 for edge and 2 for nodes on the faces
	std::ifstream mshfile(filename);
	if (mshfile.fail()) {
		std::cout << "Input file opening failed.\n";
		exit(1);
	}

	// *************** Now read in the Boundaries field: locate the keyword $PhysicalNames *****************
	//-----------------------------------------------------------------------------------------
	tmp = locate_in_file(mshfile, "$PhysicalNames");
	mshfile >> N_boundary;
	boundaries.resize(N_boundary);
	std::vector<unsigned int> tmp_boundary_tag(N_boundary); //temporary vector to store the boundary indices in MSH file (it is often irregular)
	for (int i = 0; i < N_boundary; ++i) {
		mshfile >> tmp >> tmp_boundary_tag[i] >> boundaries[i].name;
		boundaries[i].name.erase(boundaries[i].name.size() - 1); //get rid of the "" that are in the name field read from MSH file
		boundaries[i].name.erase(0, 1);
	}

	// *************** Now read in the entities field: locate the keyword $Entities *****************
	//-----------------------------------------------------------------------------------------
	tmp = locate_in_file(mshfile, "$Entities");
	mshfile >> entities_N_points >> entities_N_curves >> entities_N_faces >> tmp; //tmp is number of vols which is 0 for 2D problem
	assert(!tmp); //double check that there is no volume in the domain, so 2D problem

	for (int i = 0; i < entities_N_points; ++i) //skip the points
		do
			getline(mshfile, temp_string);
	while (temp_string.length() == 0);


	std::vector<unsigned int> tmp_curve_entity_tag(entities_N_curves); //temporary vector to store the curves entity tag in MSH file
	std::vector<unsigned int> tmp_curves_boundary_tag(entities_N_curves); //writes the boundry tag of curves
	for (int i = 0; i < entities_N_curves; ++i) { //read in the curves entity tag and the boundry tag
		mshfile >> tmp_curve_entity_tag[i] >> double_field >> double_field >> double_field >> double_field >> double_field >> double_field
			>> tmp >> tmp_curves_boundary_tag[i];
		getline(mshfile, temp_string); //read in the rst of data which are not needed at this point
	}

	// *************** Now read in the NODES field: locate the keyword $Nodes *****************
	//-----------------------------------------------------------------------------------------
	tmp = locate_in_file(mshfile, "$Nodes");

	mshfile >> nodes_total_entities >> raw_N_nodes >> nodes_min_index >> nodes_max_index;
	unorganized_nodes_coor.resize(raw_N_nodes);
	unorganized_node_type.resize(raw_N_nodes);
	tmp = nodes_max_index - nodes_min_index + 1;
	unorganized_node_mapping.resize(tmp);

	// read in all the coordinates
	for (unsigned int node_entity = 0; node_entity < nodes_total_entities; ++node_entity) {
		mshfile >> entity_dim >> group_tag; //entity_dim=0(corners), 1(edge nodes) , 2(surface nodes); group_tag: tag of group for each entity
		mshfile >> tmp >> tag_N_nodes; //tmp=0,1 means No parametric or with parametric; NNodes: number of nodes in this tag
		std::fill(unorganized_node_type.begin() + unorganized_node_index, unorganized_node_type.begin() + unorganized_node_index + tag_N_nodes, entity_dim);
		for (unsigned int node = 0; node < tag_N_nodes; ++node) { //store the indices
			mshfile >> index;
			index -= nodes_min_index; //shift all indices s.t. the indices start from zero
			unorganized_node_mapping[index] = unorganized_node_index++;	//shifting all indices, such that min_index becomes zero index
		}
		unorganized_node_index -= tag_N_nodes; //restore to read the coordinates
		for (unsigned int node = 0; node < tag_N_nodes; ++node) { //now store the coordinates
			mshfile >> unorganized_nodes_coor[unorganized_node_index].x;
			mshfile >> unorganized_nodes_coor[unorganized_node_index++].y >> coorZ;
		}
	}

	// *************** Now read in the ELEMENTS field: locate the keyword $Elements *****************
	//-----------------------------------------------------------------------------------------------
	unsigned int edge_index = 0; //to store in the boundaries, the edges that form each boundary curve
	//std::vector<std::vector<unsigned int>> tmp_curves_edges_tag(entities_N_curves); //store the tags for the edges forming each curve (curve is an entity, edgeis discretized form of curves)
	tmp = locate_in_file(mshfile, "$Elements");

	mshfile >> elements_total_entities >> N_elements >> elements_min_index >> elements_max_index; //N_element= total number of nodes, edges and 2d elements, ignore the 0d elements
	assert(elements_total_entities == entities_N_points + entities_N_curves + entities_N_faces);
	for (unsigned int element_entity = 0; element_entity < elements_total_entities; ++element_entity) {
		mshfile >> entity_dim >> group_tag; ////entity_dim=0(0d), 1(1d) , 2(2d) features; group_tag: tag of entity
		mshfile >> element_type >> tag_N_elements; //1,8,26,27,28: 2-node, 3-node, 4-node, 5-node and 6-node lines; 3,10,16,36,37: 4-node, 9-node, 8-node, 16-node, 25-node 2D elements
		if (entity_dim == 0 /*element_type==15*/)  //single-node point
			for (unsigned int element = 0; element < tag_N_elements; ++element) //skip the nodes definitions
				mshfile >> tmp1 >> tmp2;

		else if (entity_dim == 1) { // element_type corresponds to edge
			its = std::find(edge_type_node_number[0].begin(), edge_type_node_number[0].end(), element_type);
			check_start = its != edge_type_node_number[0].end(); //true means the element is of edge type
			if (!check_start) {
				std::cout << "This edge type " << element_type << " is not supported, remesh" << std::endl;
				return 2;
			}

			index = its - edge_type_node_number[0].begin();
			_edge.N_nodes = edge_type_node_number[1][index]; //number of nodes for this edge type element
			_edge.edge_type = element_type; //type of edge based on the value written in the msh file

			// in general group_tag (curve_tag here) forming edges can be written NOT in the same order as in the entity section. So, find the tmp_curve_entity_tag index that has group_tag
			its = std::find(tmp_curve_entity_tag.begin(), tmp_curve_entity_tag.end(), group_tag);
			index = its - tmp_curve_entity_tag.begin();  //index of the curve_tag (group_tag) in the vector tmp_curve_entity_tag (to find the corresponding boundary tag)


			int curve_boundary_tag = tmp_curves_boundary_tag[index];
			its = std::find(tmp_boundary_tag.begin(), tmp_boundary_tag.end(), curve_boundary_tag);
			index = its - tmp_boundary_tag.begin();  //index of the boundary_tag in the vector tmp_boundary_tag which corresponds to the curve_boundary_tag. it is used to extract the boundary name below
			boundaries[index].N_edges += tag_N_elements; //number of edges that from the boundary
			boundaries[index].edges.resize(boundaries[index].N_edges); //open up more space at the end of the vector, in case a boundary is made of multiple curves (and each curve has multiple edges)

			//tmp_curves_edges_tag[index].resize(tag_N_elements); //to store the edges tags for the curves (which are defined in entities)

			for (unsigned int i = boundaries[index].N_edges - tag_N_elements; i < boundaries[index].N_edges; ++i) {
				boundaries[index].edges[i] = edge_index++; // the index corresponding to edges vector
				//boundaries[index].edges.push_back(edge_index++); // the index corresponding to edges vector

				mshfile >> tmp;
				_edge.nodes.clear();
				for (unsigned int j = 0; j < _edge.N_nodes; ++j) {
					mshfile >> node_index;
					_edge.nodes.push_back(node_index - nodes_min_index);
				}
				edges.push_back(_edge);
			}
		}
		else if (entity_dim == 2) { // element_type corresponds to face
			its = std::find(face_type_node_number[0].begin(), face_type_node_number[0].end(), element_type);
			check_start = its != face_type_node_number[0].end(); //true means the element is of face type
			if (!check_start) {
				std::cout << "The element type " << element_type << " is not supported, remesh" << std::endl;
				return 3;
			}
			index = its - face_type_node_number[0].begin();
			_face.N_nodes = face_type_node_number[1][index]; //number of nodes for this face type element
			_face.element_type = element_type; //type of face based on the value written in the msh file
			for (unsigned int i = 0; i < tag_N_elements; ++i) {
				mshfile >> tmp;
				_face.nodes.clear();
				for (unsigned int j = 0; j < _face.N_nodes; ++j) {
					mshfile >> index;
					_face.nodes.push_back(index - nodes_min_index);
				}
				elements.push_back(_face);
			}
		}
	}

	// Now that the edge and elements are stored, the nodes constituting them should be renumbered
	// use only the nodes that are used in the edges and elements vectors
	node _node;
	node_index = 0;
	nodes.reserve(raw_N_nodes); //maximum possible number of nodes, some may not be used in the mesh
	std::vector<bool> used_node(raw_N_nodes, false); //to keep track of the nodes in the edges and elements that are used once
	std::vector<unsigned int> second_mapping(raw_N_nodes);
	for (unsigned int i = 0; i < edges.size(); ++i) {
		for (unsigned int j = 0; j < edges[i].N_nodes; ++j) {
			tmp = edges[i].nodes[j]; //tmp should map into node_index now
			tmp1 = unorganized_node_mapping[tmp];
			if (!used_node[tmp1]) {
				used_node[tmp1] = true;
				second_mapping[tmp1] = node_index++;
				_node.node_type = unorganized_node_type[tmp1];
				_node.coor = unorganized_nodes_coor[tmp1];
				nodes.emplace_back(_node);
			}
			edges[i].nodes[j] = second_mapping[tmp1];  //redistributing the indices

		}
	}
	for (unsigned int i = 0; i < elements.size(); ++i) {
		for (unsigned int j = 0; j < elements[i].N_nodes; ++j) {
			tmp = elements[i].nodes[j]; //tmp should map into index now
			tmp1 = unorganized_node_mapping[tmp];

			if (!used_node[tmp1]) {
				used_node[tmp1] = true;
				second_mapping[tmp1] = node_index++;
				_node.node_type = unorganized_node_type[tmp1];
				_node.coor = unorganized_nodes_coor[tmp1];
				nodes.emplace_back(_node);
			}
			elements[i].nodes[j] = second_mapping[tmp1];  //redistributing the indices
		}
	}

	std::cout << "         Done" << std::endl;
	mshfile.close();

	for (int i = 0; boundaries.size(); ++i)
		N_edges_boundary += boundaries[i].N_edges; //total number of edges on boundaries

	return retval;
}

int Mesh::locate_in_file(std::ifstream& filestream, const std::string& searchname) {
	//to find a specific keyword in the MSH file and return the file stream
	std::string temp;
	bool found = false;
	while (!filestream.eof()) {
		getline(filestream, temp);
		found = false;
		if (temp == searchname) {
			found = true;
			break;
		}
	}
	if (!found) {
		std::cout << "The  " << searchname << "  Field Not Found! " << std::endl;
		return 10;
	}
	return 0;
}

char Mesh::setup_mesh_problem(unsigned int prob_type) {
	//setup the mesh for problem types other than problem_type==10
	double xL = 1., yL = 1.; //domain size
	double dx, dy; //mesh spacings between nodes when they are uniformly distributed
	bool periodic=false;
	unsigned char grid_type = 1; //hard-coded for now
	unsigned int nd = 0;
	Cmpnts2 average_spacing;;

	if (fabs(dx_ratio - 1.0) > 1.e-4) grid_type = 4;
	if (prob_type <= 2) periodic = true;
	if (prob_type > 6) grid_type = 3;

	N_nodes = (N_el_i * Lnod_in + 1) * (N_el_j * Lnod_in + 1); //Haji: total number of geometry nodes in the domain
	if (periodic || grid_type==3) N_nodes = (N_el_i * Lnod_in) * (N_el_j * Lnod_in + 1);
	nodes.resize(N_nodes);
	N_edges_boundary = 2 * (N_el_i + N_el_j);  //number of edges on the boundary
	N_el = N_el_i * N_el_j; //total number of elements
	elem_neighbor = new element_neighbor [N_el]; //resize the elem_neighbors arrays
	//for (int i = 0; i < N_el; ++i) elem_neighbor[i] = new int[4]; //4 neighbors for an element
	boundary_elem_ID = new boundary_element[N_edges_boundary];
	boundary_node_ID = new unsigned int* [N_edges_boundary];
	for (int i = 0; i < N_edges_boundary; ++i) boundary_node_ID[i] = new unsigned int[Lnod]; //Lnod nodes per boundary edge
	node_ID = new unsigned int* [(int)N_el]; //resize the node_ID arrays
	for (int i = 0; i < N_el; ++i) node_ID[i] = new unsigned int[Lnod * Lnod]; //Lnod*Lnod nodes per element


	/*
	A few grid types
		square box w / random perturbation = 1
		square box w / sinusoidal perturbation = 2
		a pair of concentric circles = 3
		square boxes 1 and 2 for periodic BC
		*/

	if (prob_type != 1 && prob_type != 2 && prob_type != 6) {
		//std::cout << "The problem under consideration doesnt exist!" << std::endl;
		//exit(2);
	}

	if (prob_type == 1) {
		xL = 0.2;
		yL = 1.0;
	}

	if (grid_type == 1) { //a regular cartesian grid in a 1*1 domain but with perturbed coordinates with random distribution with magnitude "fac"
		if (Lnod > 2) {
			std::cout << "Only bilinear elements allowed for grid_type 1! " << std::endl;
			exit(3);
		}

		dx = xL / ((double)N_el_i * Lnod_in);
		dy = yL / ((double)N_el_j * Lnod_in);
		average_spacing.set_coor(dx, dy);
		//Set geometry nodes(not the same as solution nodes)
		nodes[nd++].coor.set_coor(0., 0.);  //Haji: SW corner
		for (int i = 1; i < N_el_i; ++i) //South edge
			nodes[nd++].coor.set_coor(i * dx, 0.0);

		nodes[nd++].coor.set_coor(xL, 0.0); // SE corner

		for (int j = 1; j < N_el_j; ++j) {
			nodes[nd++].coor.set_coor(0.0, j * dy); //west edge
			for (int i = 1; i < N_el_i; ++i) {
				nodes[nd++].coor.set_coor(i * dx, j * dy); //inside the domain
				if (fac < 1.0)
					nodes[nd - 1].coor.plus(fac * (0.5 - rand() / RAND_MAX), average_spacing);
			}
			nodes[nd++].coor.set_coor(xL, j * dy); //east edge
		}
		nodes[nd++].coor.set_coor(0.0, yL); //NW corner
		for (int i = 1; i < N_el_i; ++i) //north edge
			nodes[nd++].coor.set_coor(i * dx, yL);
		nodes[nd++].coor.set_coor(xL, yL);  //NE corner
	}

	else if (grid_type == 2) {  //grid extending in a sinusoidal way
		dx = xL / ((double)N_el_i * Lnod_in);
		dy = yL / ((double)N_el_j * Lnod_in);

		for (int j = 0; j < Lnod_in * N_el_j; ++j) {
			for (int i = 0; i < Lnod_in * N_el_i; ++i)
				nodes[nd++].coor.set_coor(i * dx + sin(i * 2. * pi * dx) * sin(j * 2. * pi * dy) * fac, j * dy + sin(i * 2. * pi * dx) * sin(j * 2. * pi * dy) * fac);
			//already makes a straight boundaries for the west and south boundaries only when i=j=0
			nodes[nd++].coor.set_coor(xL, j * dy); //east edge
		}

		for (int i = 0; i <= Lnod_in * N_el_i; ++i)  //North boundary
			nodes[nd++].coor.set_coor(i * dx, yL);

		//nodes[nd++].coor.set_coor(xL, yL);	//NE corner
	}

	else if (grid_type == 3) { //makes the grid in the polar coordinate (r,teta). i dir=teta, j dir= radial
		double R_in = 0.5, R_out = 1.0, R_j, percent, theta_in, theta_out, theta_i;
		double dr = (R_out - R_in) / ((double)N_el_j * Lnod_in);
		std::cout << "enter percent rotation (0 --> 1): ";
		std::cin >> percent;
		theta_in = 0.0 + percent * 2.0 * pi;
		theta_out = 2.0 * pi * (1.0 + percent); //make it a full 360 rotation
		double dtheta = (theta_out - theta_in) / ((double)N_el_i * Lnod_in);

		for (int j = 0; j <= N_el_j * Lnod_in; ++j)
			for (int i = 0; i < N_el_i * Lnod_in; ++i) {
				theta_i = theta_in + i * dtheta;
				R_j = R_in + j * dr;
				nodes[nd++].coor.set_coor(R_j * cos(-theta_i), R_j * sin(-theta_i)); //clockwise rotation of theta gridlines
			}
	}

	else if (grid_type == 4) { //shrinking and then expanding grid in x, y direction by a factor dxrat
		if (Lnod_in != 1) {
			std::cout << "only linear elements allowed for grid_type 4! " << std::endl;
			exit(4);
		}
		if (N_el_i % 2 || N_el_j % 2) {
			std::cout << "only even element counts allowed for grid type 4" << std::endl;
			exit(5);
		}

		double* geometric_X = new double[(int)N_el_i + 1]; //a 1D array to temporarily store the node locations in x direction
		double* geometric_Y = new double[(int)N_el_j + 1]; //a 1D array to temporarily store the node locations in y direction
		dx = 0.5 * xL * (1.0 - dx_ratio) / (1.0 - pow(dx_ratio, N_el_i / 2.)); // first element in geometrical series (t1)
		dy = 0.5 * yL * (1.0 - dx_ratio) / (1.0 - pow(dx_ratio, N_el_j / 2.)); //t1 in y direction
		std::cout << "Smallest and largest dx: " << dx << "  " << dx * pow(dx_ratio, N_el_i / 2. - 1.) << std::endl;
		geometric_X[0] = 0.; geometric_Y[0] = 0.; geometric_X[1] = dx; geometric_Y[1] = dy;
		geometric_X[N_el_i] = xL; geometric_Y[N_el_j] = yL; geometric_X[N_el_i - 1] = xL - dx; geometric_Y[N_el_j - 1] = yL - dy;

		for (int i = 2; i <= N_el_i / 2; ++i) {
			geometric_X[i] = geometric_X[i - 1] + (geometric_X[i - 1] - geometric_X[i - 2]) * dx_ratio;
			geometric_X[N_el_i - i] = xL - geometric_X[i];
		}

		for (int j = 2; j <= N_el_j / 2; ++j) {
			geometric_Y[j] = geometric_Y[j - 1] + (geometric_Y[j - 1] - geometric_Y[j - 2]) * dx_ratio;
			geometric_Y[N_el_j - j] = yL - geometric_Y[j];
		}

		for (int j = 0; j <= N_el_j; ++j)
			for (int i = 0; i <= N_el_i; ++i)
				nodes[nd++].coor.set_coor(geometric_X[i], geometric_Y[j]);

		delete[] geometric_X;
		delete[] geometric_Y;
	}

	else {
		std::cout << "selected grid_type is not supported!" << std::endl;
		exit(6);
	}

	// ********** set the node_ID that makes elements ***********
	unsigned int nc = 0;  //cell counter
	nd = 0; //node index counter
	for (int j = 0; j < N_el_j; ++j) {
		for (int i = 0; i < N_el_i; ++i) {
			node_ID[nc][0] = nd; //SW corner
			node_ID[nc][1] = nd + Lnod_in;  //SE corner
			node_ID[nc][2] = nd + Lnod_in * (Lnod_in * N_el_i + 2); // NE corner (indices go row by row)
			node_ID[nc][3] = nd + Lnod_in * (Lnod_in * N_el_i + 1); //NW corner
			
			if (grid_type == 3) { //it has one less node in the theta direction, so adjust
				node_ID[nc][2] -= Lnod_in;
				node_ID[nc][3] -= Lnod_in;
				if (i == N_el_i - 1) {
					node_ID[nc][1] -= Lnod_in * N_el_i;  //or equivalently node_ID[nc-N_el_i+1][0]
					node_ID[nc][2] -= Lnod_in * N_el_i;
				}
			}

			if (Lnod_in == 2) { //set the rest of nodes or higher order elements
				node_ID[nc][4] = node_ID[nc][0] + 1;
				node_ID[nc][6] = node_ID[nc][3] + 1;
				node_ID[nc][7] = (node_ID[nc][0] + node_ID[nc][3]) / 2; //average of local node 0,3. This way no need to add special case for periodic condition
				node_ID[nc][5] = (node_ID[nc][1] + node_ID[nc][2]) / 2; //average of local node 1,2
				node_ID[nc][8] = node_ID[nc][7] + 1;
			}
			else if (Lnod_in == 3) {
				node_ID[nc][4] = node_ID[nc][0] + 1;
				node_ID[nc][5] = node_ID[nc][0] + 2;
				node_ID[nc][8] = node_ID[nc][3] + 2;
				node_ID[nc][9] = node_ID[nc][3] + 1;
				node_ID[nc][11] = (2*node_ID[nc][0] + node_ID[nc][3]) / 3; //get rid of special case of periodic condition
				node_ID[nc][10] = (node_ID[nc][0] + 2*node_ID[nc][3]) / 3;
				node_ID[nc][6] = (2*node_ID[nc][1] + node_ID[nc][2]) / 3; //get rid of special case of periodic condition
				node_ID[nc][7] = (node_ID[nc][1] + 2*node_ID[nc][2]) / 3; //get rid of special case of periodic condition
				node_ID[nc][12] = node_ID[nc][11] + 1;
				node_ID[nc][13] = node_ID[nc][11] + 2;
				node_ID[nc][14] = node_ID[nc][10] + 2;
				node_ID[nc][15] = node_ID[nc][10] + 1;
			}
			else if (Lnod_in > 3) {
				std::cout << "ERROR: Up to bicubic elements supported! Set Lnod = 3, 2 or 1" << std::endl;
				exit(7);
			}

			nd += Lnod_in;
			nc++;
		}
		nd += (Lnod_in - 1) * (Lnod_in * N_el_i + 1) + 1;
		if (grid_type == 3) nd -= Lnod_in;
	}
	// **************************************************************

	// *********************** set the boundary elements and nodes of cells **********************
	/*
	Set mesh connectivity to nodes(nodeID); also get connectivity to neighboring meshes(elemID)
	nodeID starts from SW cornersand goes CCW; elemID starts with S faceand goes CCW
	NOTE: This is of course a poor man's (structured grid) approach to the general problem
	Also, the - ve numbers refer to boundary faces, which are assigned in the BC_IC setup routine
	*/
	int nb = 0;  //counter for the edges locaed on the boundaries of the domain
	nc = 0;  //cell counter
	nd = 0; //node index
	int bool_boundary = false;
	if (Lnod_in > 3) {
		std::cout << "ERROR: Up to bicubic elements supported! Set Lnod = 3, 2 or 1" << std::endl;
		exit(7);
	}
	for (int j = 0; j < N_el_j; ++j) 
		for (int i = 0; i < N_el_i; ++i) {
			bool_boundary = false;
			if (!j) { //the cells attached to the south boundary of the domain
				elem_neighbor[nc].is_on_boundary[south]=true;
				elem_neighbor[nc].boundary_index[south] = nb;
				boundary_elem_ID[nb].element_index = nc;
				boundary_elem_ID[nb].side = south;
				boundary_node_ID[nb][0] = node_ID[nc][0];
				boundary_node_ID[nb][1] = node_ID[nc][1];
				if (Lnod_in == 2) boundary_node_ID[nb][2] = node_ID[nc][4];
				else if (Lnod_in == 3) {
					boundary_node_ID[nb][2] = node_ID[nc][4];
					boundary_node_ID[nb][3] = node_ID[nc][5];
				}
				nb++;
			}
			else elem_neighbor[nc].neighbor[south] = nc - N_el_i; //index of the neighboring element in the South

			if (j == N_el_j - 1) { //the cells attached to the north boundary of the domain
				elem_neighbor[nc].is_on_boundary[north] = true;
				elem_neighbor[nc].boundary_index[north] = nb;  //fix this later, it is wrong as of now
				boundary_elem_ID[nb].element_index = nc;
				boundary_elem_ID[nb].side = north;
				boundary_node_ID[nb][0] = node_ID[nc][3];
				boundary_node_ID[nb][1] = node_ID[nc][2];
				if (Lnod_in == 2) boundary_node_ID[nb][2] = node_ID[nc][6];
				else if (Lnod_in == 3) {
					boundary_node_ID[nb][2] = node_ID[nc][9];
					boundary_node_ID[nb][3] = node_ID[nc][8];
				}
				nb++;
			}
			else elem_neighbor[nc].neighbor[north] = nc + N_el_i; //index of the neighboring element in the South

			if (!i) { //the cells attached to the west boundary of the domain
				if (grid_type == 3 || periodic)
					elem_neighbor[nc].neighbor[west] = nc + N_el_i - 1;
				else {
					elem_neighbor[nc].is_on_boundary[west] = true;
					elem_neighbor[nc].boundary_index[west] = nb;
					boundary_elem_ID[nb].element_index = nc;
					boundary_elem_ID[nb].side = west;
					//CCW rotation for meshing, For now, I'm assuming all sides are in positive x/y direction. NOTE: need to account for direction later
					boundary_node_ID[nb][0] = node_ID[nc][0];
					boundary_node_ID[nb][1] = node_ID[nc][3];
					if (Lnod_in == 2) boundary_node_ID[nb][2] = node_ID[nc][7];
					else if (Lnod_in == 3) {
						boundary_node_ID[nb][2] = node_ID[nc][11];
						boundary_node_ID[nb][3] = node_ID[nc][10];
					}
					nb++;
				}
			}
			else elem_neighbor[nc].neighbor[west] = nc - 1; //index of the neighboring element in the South

			if (i == N_el_i - 1) { //the cells attached to the east boundary of the domain
				if (grid_type == 3 || periodic)
					elem_neighbor[nc].neighbor[east] = nc - N_el_i +1;
				else {
					elem_neighbor[nc].is_on_boundary[east] = true;
					elem_neighbor[nc].boundary_index[east] = nb;
					boundary_elem_ID[nb].element_index = nc;
					boundary_elem_ID[nb].side = east;
					boundary_node_ID[nb][0] = node_ID[nc][1];
					boundary_node_ID[nb][1] = node_ID[nc][2]; 
					if (Lnod_in == 2) boundary_node_ID[nb][2] = node_ID[nc][5];
					else if (Lnod_in == 3) {
						boundary_node_ID[nb][2] = node_ID[nc][6];
						boundary_node_ID[nb][3] = node_ID[nc][7];
					}
					nb++;
				}
			}
			else elem_neighbor[nc].neighbor[east] = nc + 1; //index of the neighboring element in the South

			nc++;
			//if (std::any_of(std::begin(elem_neighbor[nc].is_on_boundary), std::end(elem_neighbor[nc].is_on_boundary), [](bool i) {return i; } ) ) nb++; //update the 
		}  //for i < N_el_i , for j < N_el_j

	N_edges_boundary = nb; //the corrected number of edges on the boundary (it exclude the edge woth periodic BC)

}

