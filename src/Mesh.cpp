/*
 * Mesh.cpp - Unstructured mesh class implementation
 *
 * (c)2020-1 Applied Scientific Research, Inc.
 *           Mohammad Hajit <mhajit@gmail.com>
 */

#include "Mesh.hpp"

#include <cmath>


int Mesh::tensor2FEM(int i) {
	// converts the tensor index to the FEM node ordering for 1D case
	if (!i) return 0;
	else if (i == Lnod_in) return 1;
	else return (i + 1);
}

int Mesh::tensor2FEM(int i, int j) {
	// converts the tensor index to the FEM node ordering for 2D case
	int t2f = -1;

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

char Mesh::read_msh_file() {
	// reads a msh file output from the Gmsh software. The msh file is in ASCII 4.1 version of the Gmsh output
	std::string filename = input_msh_file;
	std::cout << "     Gmsh file   ***** " << filename << " *****   opened for reading ..." << std::endl << std::endl;
	int retval = 1;
	int index;
	int nodes_min_index, nodes_max_index, raw_N_nodes, tag_N_nodes, nodes_total_entities, group_tag, entity_dim, unorganized_node_index = 0;
	int elements_min_index, elements_max_index, tag_N_elements, elements_total_entities, element_type;
	int node_index;
	int entities_N_points, entities_N_curves, entities_N_faces;
	double double_field;
	std::string temp_string;
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
	int tmp = locate_in_file(mshfile, "$PhysicalNames");
	if (tmp==10) {std::cout<< "could not find the physicalNames field, check mesh \n"; exit(1);}
	int N_physical_names, dim_physical_name;
	mshfile >> N_physical_names;  //number of global boundaries with boundary conditions
	std::getline(mshfile, temp_string);
	std::vector<unsigned int> tmp_boundary_tag; //temporary vector to store the boundary indices in MSH file (it is often irregular)
	for (int i = 0; i < N_physical_names; ++i) {
        mshfile >> dim_physical_name;
        if (dim_physical_name!=1) continue; //we are not interested in 2D or 0D physical names
        tmp_boundary_tag.resize(N_Gboundary+1);
        boundaries.resize(N_Gboundary+1);
		mshfile >> tmp_boundary_tag[N_Gboundary] >> boundaries[N_Gboundary].name; std::getline(mshfile, temp_string);
		boundaries[N_Gboundary].name.erase(boundaries[N_Gboundary].name.size() - 1); //get rid of the "" that are in the name field read from MSH file
		boundaries[N_Gboundary].name.erase(0, 1);
		N_Gboundary++;
	}

	// *************** Now read in the entities field: locate the keyword $Entities *****************
	//-----------------------------------------------------------------------------------------
	tmp = locate_in_file(mshfile, "$Entities");
	if (tmp==10) {std::cout<< "could not find the Entities field, check mesh \n"; exit(1);}
	mshfile >> entities_N_points >> entities_N_curves >> entities_N_faces >> tmp; //tmp is number of vols which is 0 for 2D problem
	std::getline(mshfile, temp_string);
	assert(!tmp); //double check that there is no volume in the domain, so 2D problem

	for (int i = 0; i < entities_N_points; ++i) std::getline(mshfile, temp_string); //skip the points

	std::vector<unsigned int> tmp_curve_entity_tag(entities_N_curves); //temporary vector to store the curves entity tag in MSH file
	std::vector<unsigned int> tmp_curves_boundary_tag(entities_N_curves); //writes the boundry tag of curves
	std::vector<int> tmp_curves_num_physical_tag(entities_N_curves); //this curve entity belongs to how many physical boundaries. If zero, then it is an interface between 2 subdomains.
	for (int i = 0; i < entities_N_curves; ++i) { //read in the curves entity tag and the boundry tag
	  mshfile >> tmp_curve_entity_tag[i] >> double_field >> double_field >> double_field >> double_field >> double_field >> double_field >> tmp_curves_num_physical_tag[i] >> tmp_curves_boundary_tag[i];
	  std::getline(mshfile, temp_string); //read in the rst of data which are not needed at this point
	}

	// *************** Now read in the NODES field: locate the keyword $Nodes *****************
	//-----------------------------------------------------------------------------------------
	tmp = locate_in_file(mshfile, "$Nodes");
	if (tmp==10) {std::cout<< "could not find the nodes field, check mesh \n"; exit(1);}
	mshfile >> nodes_total_entities >> raw_N_nodes >> nodes_min_index >> nodes_max_index;
	std::getline(mshfile, temp_string);
	unorganized_nodes_coor.resize(raw_N_nodes);
	unorganized_node_type.resize(raw_N_nodes);
	tmp = nodes_max_index - nodes_min_index + 1;
	unorganized_node_mapping.resize(tmp);

	// read in all the coordinates
	for (int node_entity = 0; node_entity < nodes_total_entities; ++node_entity) {
		mshfile >> entity_dim >> group_tag; //entity_dim=0(corners), 1(edge nodes) , 2(surface nodes); group_tag: tag of group for each entity
		mshfile >> tmp >> tag_N_nodes; //tmp=0,1 means No parametric or with parametric; NNodes: number of nodes in this tag
		std::getline(mshfile, temp_string);
		std::fill(unorganized_node_type.begin() + unorganized_node_index, unorganized_node_type.begin() + unorganized_node_index + tag_N_nodes, entity_dim);
		for (int node = 0; node < tag_N_nodes; ++node) { //store the indices
		  mshfile >> index; std::getline(mshfile, temp_string);
			index -= nodes_min_index; //shift all indices s.t. the indices start from zero
			unorganized_node_mapping[index] = unorganized_node_index++;	//shifting all indices, such that min_index becomes zero index
		}
		unorganized_node_index -= tag_N_nodes; //restore to read the coordinates
		for (int node = 0; node < tag_N_nodes; ++node) { //now store the coordinates
		  double a,b,c;
		  mshfile >> a >>b >>c; std::getline(mshfile, temp_string); // The x,y,z coordinates of the nodes
		  unorganized_nodes_coor[unorganized_node_index++].set_coor(a,b);
		}
	}

	// *************** Now read in the ELEMENTS field: locate the keyword $Elements *****************
	//-----------------------------------------------------------------------------------------------
	unsigned int edge_index = 0; //to store in the boundaries, the edges that form each boundary curve
	//std::vector<std::vector<unsigned int>> tmp_curves_edges_tag(entities_N_curves); //store the tags for the edges forming each curve (curve is an entity, edgeis discretized form of curves)
	tmp = locate_in_file(mshfile, "$Elements");
	if (tmp==10) {std::cout<< "could not find the Elements field, check mesh \n"; exit(1);}
	mshfile >> elements_total_entities >> N_elements >> elements_min_index >> elements_max_index; //N_element= total number of nodes, edges and 2d elements, ignore the 0d elements
	std::getline(mshfile, temp_string);

	assert(elements_total_entities == entities_N_points + entities_N_curves + entities_N_faces);
	for (int element_entity = 0; element_entity < elements_total_entities; ++element_entity) {
		mshfile >> entity_dim >> group_tag; ////entity_dim=0(0d), 1(1d) , 2(2d) features; group_tag: tag of entity
		mshfile >> element_type >> tag_N_elements; //1,8,26,27,28: 2-node, 3-node, 4-node, 5-node and 6-node lines; 3,10,16,36,37: 4-node, 9-node, 8-node, 16-node, 25-node 2D elements
		std::getline(mshfile, temp_string);
		if (entity_dim == 0 /*element_type==15*/) { //single-node point
			for (int element = 0; element < tag_N_elements; ++element) { //skip the nodes definitions
				std::getline(mshfile, temp_string);
			}

		} else if (entity_dim == 1) { // element_type corresponds to edge
			//number of nodes for this edge type element
			try {
				_edge.N_nodes = edge_type_node_number.at(element_type);
			} catch (...) {
				std::cout << "The element type " << element_type << " is not supported, remesh" << std::endl;
				return 2;
			}

			_edge.edge_type = element_type; //type of edge based on the value written in the msh file

			// in general group_tag (curve_tag here) forming edges can be written NOT in the same order as in the entity section. So, find the tmp_curve_entity_tag index that has group_tag
			std::vector<unsigned int>::iterator its =
				std::find(tmp_curve_entity_tag.begin(), tmp_curve_entity_tag.end(), group_tag);
			index = its - tmp_curve_entity_tag.begin();  //index of the curve_tag (group_tag) in the vector tmp_curve_entity_tag (to find the corresponding boundary tag)

            if (!tmp_curves_num_physical_tag[index]) {  //exclude/skip the internal interfaces
                for (int i=0; i<tag_N_elements; ++i) std::getline(mshfile, temp_string);
                    continue;
            }

			int curve_boundary_tag = tmp_curves_boundary_tag[index]; // the physical tag of this
			its = std::find(tmp_boundary_tag.begin(), tmp_boundary_tag.end(), curve_boundary_tag);
			index = its - tmp_boundary_tag.begin();  //index of the boundary_tag in the vector tmp_boundary_tag which corresponds to the curve_boundary_tag. it is used to extract the boundary name below
			boundaries[index].N_edges += tag_N_elements; //number of edges that form the boundary
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
				std::getline(mshfile, temp_string);
				edges.push_back(_edge);
			}
		}
		else if (entity_dim == 2) { // element_type corresponds to face
			//number of nodes for this face type element, throws exception if it does not exist
			try {
				_face.N_nodes = face_type_node_number.at(element_type);
			} catch (...) {
				std::cout << "The element type " << element_type << " is not supported, remesh" << std::endl;
				return 3;
			}

			_face.element_type = element_type; //type of face based on the value written in the msh file
			for (int i = 0; i < tag_N_elements; ++i) {
				mshfile >> tmp;
				_face.nodes.clear();
				for (unsigned int j = 0; j < _face.N_nodes; ++j) {
					mshfile >> index;
					_face.nodes.push_back(index - nodes_min_index);
				}
				std::getline(mshfile, temp_string);
				elements.push_back(_face);
			}
		}
	}
	mshfile.close();

	Lnod = _edge.N_nodes;
	Lnod_in = Lnod - 1;

	// ******** Reorient in the CCW direction if the mesh file is in OCC format ************
	std::vector<unsigned int> tmp_node_ids(Lnod*Lnod);
    if (OCC_format) {
        for (int el=0; el<(int)elements.size(); ++el) {
            tmp_node_ids = elements[el].nodes;
            for (int j_old=0; j_old<Lnod; ++j_old)
                for (int i_old=0; i_old<Lnod; ++i_old) {
                    int j_new = Lnod_in - i_old;
                    int i_new = Lnod_in - j_old;
                    elements[el].nodes[tensor2FEM(i_new,j_new)] = tmp_node_ids[tensor2FEM(i_old,j_old)];
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
			const int tmp = edges[i].nodes[j]; //tmp should map into node_index now
			const int tmp1 = unorganized_node_mapping[tmp];
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
			const int tmp = elements[i].nodes[j]; //tmp should map into index now
			const int tmp1 = unorganized_node_mapping[tmp];

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

	//std::cout << "         Done" << std::endl;
	//	std::cout << elements.size();
	//std::cout << "out of the loop" << elements[10].nodes[5] << " "<<elements[20].nodes[6] << " "<< elements[30].nodes[4] << " ";

	process_mesh();

	return retval;
}

void Mesh::process_mesh() {
	/*
	processes the mesh that is read from file. It finds the elements neighbors and stores it in elem_neighbors,
	*/
	N_nodes = nodes.size(); //total number of nodes in the domain
	N_edges_boundary = 0;
	for (int i = 0; i < (int)boundaries.size(); ++i) {
		std::cout << "adding " << boundaries[i].N_edges << " edges from boundary " << i << std::endl;
		N_edges_boundary += boundaries[i].N_edges; //total number of edges on global boundaries
	}
	N_el = elements.size();  //total number of elements (2d faces)
	elem_neighbor = new element_neighbor[N_el];
	boundary_elem_ID = new boundary_element[N_edges_boundary];
	std::cout << "  edges vector begins with " << edges.size() << " elements" << std::endl;

	// *************** detect the total interior edges in the domain and add them to the edges vector (which already has boundary edges) and add them to the elements vector *********************
	edge _edge;
	bool found;
	unsigned int edge_index = edges.size();
	for (int el = 0; el < N_el; ++el) {  //loop over all 2D elements
		const unsigned int element_type = elements[el].element_type;
		//number of nodes for this edge type on the edges of the element el
		const unsigned int N_edge_nodes = element_edge_node_number.at(element_type);
		//type of edge on the boundary of the element el
		const unsigned int edge_type = edge_type_node_number_inv.at(N_edge_nodes);

		_edge.N_nodes = N_edge_nodes;
		_edge.edge_type = edge_type;
		_edge.nodes.resize(N_edge_nodes);
		for (int s = south; s <= west; s++) { //loop over all 4 sides of each element to find the neighbors
			_edge.nodes[0] = elements[el].nodes[s];
			_edge.nodes[1] = elements[el].nodes[(s + 1) % 4];
			//_edge.nodes.insert(_edge.nodes.begin() + 2, elements[el].nodes.begin() + s * (N_edge_nodes - 2) + 4, elements[el].nodes.begin() + (s + 1) * (N_edge_nodes - 2) + 4);
			std::copy(elements[el].nodes.begin() + s * (N_edge_nodes - 2) + 4, elements[el].nodes.begin() + (s + 1) * (N_edge_nodes - 2) + 4, _edge.nodes.begin() + 2);

			//search edge _edge in the edges vector to check if it doesnt already exist
			found = false; //if the _edge exists in the edges list already
			for (int edg = 0; edg < (int)edges.size(); ++edg) {
				if ((edges[edg].nodes == _edge.nodes) ||
					(edges[edg].nodes[1] == _edge.nodes[0] && edges[edg].nodes[0] == _edge.nodes[1])) { //this edge already exists in the list
					elements[el].edges[s] = edg; // the index of the edge that exists
					// This is added on April 8, 2021: noticed the boundary edges may rotate in CW or CCW depending on how they are created. Lets adopt a standard. As we march on the boundary edges, the domain should be on our left.
					edges[edg].nodes = _edge.nodes; //making sure that the edges are oriented in the same CCW orientation of the element attach to the boundary edge.
					found = true;
					break;
				}
			}

			if (!found) {  // is the _edge does not exist in the list of edges
				edges.push_back(_edge);
				elements[el].edges[s] = edge_index++;
			}
		} //for s=south, ...
	} //for el
	std::cout << "  after adding interior edges, vector now has " << edges.size() << " elements" << std::endl;

	//*********************************************************************************************************************************************************
	// ************************************* now detect the neighbors of each element **************************************
	for (int el = 0; el < N_el; ++el) {
		for (int s = south; s <= west; s++) { //loop over all 4 sides of each element to find the neighbors
			found = false;
			for (int eln = 0; eln < N_el; ++eln) { //loop over all elements to see if the edge[s] belongs to them too
				if (eln == el) continue;
				for (int sn = south; sn <= west; sn++) { //loop over all 4 sides of neighboring element to check if any of them are common to s
					if (elements[el].edges[s] == elements[eln].edges[sn]) {
						elem_neighbor[el].neighbor[s] = eln;
						elem_neighbor[el].neighbor_common_side[s] = sn;
						found = true;
						break;
					}
				}
				if (found) break;
			}
			if (!found) { // if the side s of element el does not have a neighboring element, then side s SHOULD be on the boundary, so find the edge index
				elem_neighbor[el].is_on_boundary[s] = true;
				elem_neighbor[el].boundary_index[s] = elements[el].edges[s];
			}
		}
	}
	//*********************************************************************************************************************
	// ****************************** form the boundary_elem_ID which has the details of all edges on the global boundary *****************************
	for (int el = 0; el < N_el; ++el) {
		for (int s = south; s <= west; s++) {
			if (elem_neighbor[el].is_on_boundary[s]) { //if the side s of element el is on global boundary
				edge_index = elem_neighbor[el].boundary_index[s];
				boundary_elem_ID[edge_index].element_index = el;
				boundary_elem_ID[edge_index].side = (cell_sides)s;
			}
		}
	}
	/*
	// *********************************************************************************************************************
	// ************* rotate the edges of elements so that north of this element is neighboring the south of neighboring element, ... ***************
	//std::vector<bool> fixed_elements(N_el, false); //the elements that are fixed by rotating the edges are turned into true
	//fixed_elements[0] = true; //I assume element 0 is fixed, so orient the rest based on element 0
	std::vector<bool> listed_elements_to_orient(N_el, false); //a list of the elements that are added to the vector tet_guess_for_hex_vertices
	std::vector<boundary_element> to_be_oriented_list; //a list of the elements indices to be oriented and the side it should be fixed based upon: to be completed gradually as it runs
	to_be_oriented_list.reserve(N_el);
	listed_elements_to_orient[0] = true;
	boundary_element tmp_BE;
	//tmp_BE.element_index = 0; tmp_BE.side = south;
	to_be_oriented_list.emplace_back(0,south); //the element=0 is added and its side south should be oriented in south, sop no work on element 0

	unsigned int elem_pos = 0; //the index of element in the to_be_oriented_list vector
	unsigned int orig_edges[4];  //to make a copy of the elements[el].edges and then reorient them
	element_neighbor orig_elem_neighbor; //to make a copy of the elem_neighbor and then reorient the components
	while (elem_pos < N_el) {
		int el = to_be_oriented_list[elem_pos].element_index;
		cell_sides desired_side = to_be_oriented_list[elem_pos].side; //the side direction that we want to move the south edge on it (south edge should be oriented to this side)
		// ************ reorient such that the local south edge orients in the to_be_oriented_list[el].side
		int delta_orient = 4 - (desired_side - south); //old south index should move into this new side: old 0 index of edge should be this new index: new_edge[0] = old_edge[delta_orient], ...
		if (delta_orient <4) {
			std::copy(elements[el].edges, elements[el].edges+4, orig_edges);  // copy the edges array to shift them by desired_side amount
			orig_elem_neighbor = elem_neighbor[el];
			for (int s = south; s <= west; ++s) { //reorient the order
				elements[el].edges[s] = orig_edges[(s + delta_orient) % 4];
				elem_neighbor[el].boundary_index[s] = orig_elem_neighbor.boundary_index[(s + delta_orient) % 4];
				elem_neighbor[el].neighbor[s] = orig_elem_neighbor.neighbor[(s + delta_orient) % 4];
				elem_neighbor[el].is_on_boundary[s] = orig_elem_neighbor.is_on_boundary[(s + delta_orient) % 4];
				elem_neighbor[el].neighbor_common_side[s] = orig_elem_neighbor.neighbor_common_side[(s + delta_orient) % 4];
			}
		}
		// *******************************************************************************************
		// ************** add the 4 neighbors elements if they are not on boundary ****************
		for (int s = south; s <= west; ++s)
			if (!elem_neighbor[el].is_on_boundary[s]) {
				int eln = elem_neighbor[el].neighbor[s];
				if (!listed_elements_to_orient[eln]) {
					listed_elements_to_orient[eln] = true;
					int sn = elem_neighbor[el].neighbor_common_side[s];
					to_be_oriented_list.emplace_back(eln, (cell_sides)((6 + s - sn) % 4)); //the element index and what the south side should get projected to
				}
			}
		elem_pos++;
	}
	*/
}

int Mesh::locate_in_file(std::ifstream& filestream, const std::string& searchname) {
	//to find a specific keyword in the MSH file and return the file stream
	std::string temp;
	bool found = false;
	while (filestream.good()) {
	  std::getline(filestream, temp);
	  //temp.erase(temp.size() - 1); //get rid of the "" that are in the name field read from MSH file (is only the case in windows)
	  // temp.erase(0, 1); // get rid of the first character
	  std::cout << temp << std::endl;
	  found = false;
	  if (temp==searchname) {
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
	//setup the mesh manually for problem types other than problem_type==10
	double xL = 1., yL = 1.; //domain size
	double dx, dy; //mesh spacings between nodes when they are uniformly distributed
	bool periodic=false;
	unsigned char grid_type = 1; //hard-coded for now
	unsigned int nd = 0;
	Cmpnts2 average_spacing;

	if (fabs(dx_ratio - 1.0) > 1.e-4) grid_type = 4;
	if (prob_type > 6) periodic = true;
	if (prob_type > 6) grid_type = 3;  //for peridic boundary

	N_Gboundary = (grid_type == 3) ? 2 : 4; // for the cylinder case, the only boundaries are the innermost circleand outermost circle
	boundaries.resize(N_Gboundary);
	N_nodes = (periodic || grid_type == 3) ? (N_el_i * Lnod_in) * (N_el_j * Lnod_in + 1) : (N_el_i * Lnod_in + 1) * (N_el_j * Lnod_in + 1); //Haji: total number of geometry nodes in the domain
	nodes.resize(N_nodes);
	N_edges_boundary = (periodic || grid_type == 3) ? 2*N_el_i : 2 * (N_el_i + N_el_j);  //number of edges on the boundary
	N_el = N_el_i * N_el_j; //total number of elements
	elem_neighbor = new element_neighbor [N_el]; //resize the elem_neighbors arrays
	boundary_elem_ID = new boundary_element[N_edges_boundary];
	edges.resize(N_edges_boundary);
	//boundary_node_ID = new unsigned int* [N_edges_boundary];
	for (int i = 0; i < N_edges_boundary; ++i) {
		/*boundary_node_ID[i] = new unsigned int[Lnod];*/ //Lnod nodes per boundary edge
		edges[i].nodes.resize(Lnod);
		edges[i].N_nodes = Lnod;
	}
	elements.resize(N_el); //resize the elements to store the nodes forming each element
	for (int i = 0; i < N_el; ++i) elements[i].nodes.resize(Lnod * Lnod); //Lnod*Lnod nodes per element

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
				nodes[nd++].coor.set_coor(i * dx + sin(i * 2. * M_PI * dx) * sin(j * 2. * M_PI * dy) * fac*dx, j * dy + sin(i * 2. * M_PI * dx) * sin(j * 2. * M_PI * dy) * fac*dy);
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
		theta_in = 0.0 + percent * 2.0 * M_PI;
		theta_out = 2.0 * M_PI * (1.0 + percent); //make it a full 360 rotation
		double dtheta = (theta_out - theta_in) / ((double)N_el_i * Lnod_in);

		for (int j = 0; j < N_el_j * Lnod_in+1; ++j)
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

	// ********** set the node_ID that makes elements, replaced by elements vector to be consisent with read mesh file data ***********
	unsigned int nc = 0;  //cell counter
	nd = 0; //node index counter
	for (int j = 0; j < N_el_j; ++j) {
		for (int i = 0; i < N_el_i; ++i) {
			elements[nc].nodes[0] = nd; //SW corner
			elements[nc].nodes[1] = nd + Lnod_in;  //SE corner
			elements[nc].nodes[2] = nd + Lnod_in * (Lnod_in * N_el_i + 2); // NE corner (indices go row by row)
			elements[nc].nodes[3] = nd + Lnod_in * (Lnod_in * N_el_i + 1); //NW corner

			if (grid_type == 3) { //it has one less node in the theta direction, so adjust
				elements[nc].nodes[2] -= Lnod_in;
				elements[nc].nodes[3] -= Lnod_in;
				if (i == N_el_i - 1) {
					elements[nc].nodes[1] -= Lnod_in * N_el_i;  //or equivalently node_ID[nc-N_el_i+1][0]
					elements[nc].nodes[2] -= Lnod_in * N_el_i;
				}
			}

			if (Lnod_in == 2) { //set the rest of nodes or higher order elements
				elements[nc].nodes[4] = elements[nc].nodes[0] + 1;
				elements[nc].nodes[6] = elements[nc].nodes[3] + 1;
				elements[nc].nodes[7] = (elements[nc].nodes[0] + elements[nc].nodes[3]) / 2; //average of local node 0,3. This way no need to add special case for periodic condition
				elements[nc].nodes[5] = (elements[nc].nodes[1] + elements[nc].nodes[2]) / 2; //average of local node 1,2
				elements[nc].nodes[8] = elements[nc].nodes[7] + 1;
			}
			else if (Lnod_in == 3) {
				elements[nc].nodes[4] = elements[nc].nodes[0] + 1;
				elements[nc].nodes[5] = elements[nc].nodes[0] + 2;
				elements[nc].nodes[8] = elements[nc].nodes[3] + 2;
				elements[nc].nodes[9] = elements[nc].nodes[3] + 1;
				elements[nc].nodes[11] = (2*elements[nc].nodes[0] + elements[nc].nodes[3]) / 3; //get rid of special case of periodic condition
				elements[nc].nodes[10] = (elements[nc].nodes[0] + 2*elements[nc].nodes[3]) / 3;
				elements[nc].nodes[6] = (2*elements[nc].nodes[1] + elements[nc].nodes[2]) / 3; //get rid of special case of periodic condition
				elements[nc].nodes[7] = (elements[nc].nodes[1] + 2*elements[nc].nodes[2]) / 3; //get rid of special case of periodic condition
				elements[nc].nodes[12] = elements[nc].nodes[11] + 1;
				elements[nc].nodes[13] = elements[nc].nodes[11] + 2;
				elements[nc].nodes[14] = elements[nc].nodes[10] + 2;
				elements[nc].nodes[15] = elements[nc].nodes[10] + 1;
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

	// ************************* set the boundary elements and nodes of cells **********************
	/*
	Set mesh connectivity to nodes(nodeID); also get connectivity to neighboring meshes(elemID)
	nodeID starts from SW cornersand goes CCW; elemID starts with S faceand goes CCW
	NOTE: This is of course a poor man's (structured grid) approach to the general problem
	Also, the - ve numbers refer to boundary faces, which are assigned in the BC_IC setup routine
	*/
	int nb = 0;  //counter for the edges locaed on the boundaries of the domain
	nc = 0;  //cell counter
	nd = 0; //node index
	int top_B = -1;
	if (Lnod_in > 3) {
		std::cout << "ERROR: Up to bicubic elements supported! Set Lnod = 3, 2 or 1" << std::endl;
		exit(7);
	}
	if (grid_type == 3 || periodic) {  //there is no boundary in the theta direction (periodic)
		boundaries[0].name = "bottom";
		boundaries[1].name = "top";
		top_B = 1;
	}
	else {
		boundaries[0].name = "bottom";
		boundaries[1].name = "right";
		boundaries[2].name = "top";
		boundaries[3].name = "left";
		top_B = 2;
	}
	for (int j = 0; j < N_el_j; ++j)
		for (int i = 0; i < N_el_i; ++i) {
			if (!j) { //the cells attached to the south boundary of the domain
				boundaries[0].edges.push_back(nb);
				boundaries[0].N_edges++;
				elem_neighbor[nc].is_on_boundary[south]=true;
				elem_neighbor[nc].boundary_index[south] = nb;
				boundary_elem_ID[nb].element_index = nc;
				boundary_elem_ID[nb].side = south;
				edges[nb].nodes[0] = elements[nc].nodes[0];
				edges[nb].nodes[1] = elements[nc].nodes[1];
				if (Lnod_in == 2) edges[nb].nodes[2] = elements[nc].nodes[4];
				else if (Lnod_in == 3) {
					edges[nb].nodes[2] = elements[nc].nodes[4];
					edges[nb].nodes[3] = elements[nc].nodes[5];
				}
				nb++;
			}
			else elem_neighbor[nc].neighbor[south] = nc - N_el_i; //index of the neighboring element in the South

			if (j == N_el_j - 1) { //the cells attached to the north boundary of the domain
				boundaries[top_B].edges.push_back(nb);
				boundaries[top_B].N_edges++;
				elem_neighbor[nc].is_on_boundary[north] = true;
				elem_neighbor[nc].boundary_index[north] = nb;  //fix this later, it is wrong as of now
				boundary_elem_ID[nb].element_index = nc;
				boundary_elem_ID[nb].side = north;
				edges[nb].nodes[0] = elements[nc].nodes[3];
				edges[nb].nodes[1] = elements[nc].nodes[2];
				if (Lnod_in == 2) edges[nb].nodes[2] = elements[nc].nodes[6];
				else if (Lnod_in == 3) {
					edges[nb].nodes[2] = elements[nc].nodes[9];
					edges[nb].nodes[3] = elements[nc].nodes[8];
				}
				nb++;
			}
			else elem_neighbor[nc].neighbor[north] = nc + N_el_i; //index of the neighboring element in the South

			if (!i) { //the cells attached to the west boundary of the domain
				if (grid_type == 3 || periodic)
					elem_neighbor[nc].neighbor[west] = nc + N_el_i - 1;
				else {
					boundaries[3].edges.push_back(nb);
					boundaries[3].N_edges++;
					elem_neighbor[nc].is_on_boundary[west] = true;
					elem_neighbor[nc].boundary_index[west] = nb;
					boundary_elem_ID[nb].element_index = nc;
					boundary_elem_ID[nb].side = west;
					//CCW rotation for meshing, For now, I'm assuming all sides are in positive x/y direction. NOTE: need to account for direction later
					edges[nb].nodes[0] = elements[nc].nodes[0];
					edges[nb].nodes[1] = elements[nc].nodes[3];
					if (Lnod_in == 2) edges[nb].nodes[2] = elements[nc].nodes[7];
					else if (Lnod_in == 3) {
						edges[nb].nodes[2] = elements[nc].nodes[11];
						edges[nb].nodes[3] = elements[nc].nodes[10];
					}
					nb++;
				}
			}
			else elem_neighbor[nc].neighbor[west] = nc - 1; //index of the neighboring element in the South

			if (i == N_el_i - 1) { //the cells attached to the east boundary of the domain
				if (grid_type == 3 || periodic)
					elem_neighbor[nc].neighbor[east] = nc - N_el_i +1;
				else {
					boundaries[1].edges.push_back(nb);
					boundaries[1].N_edges++;
					elem_neighbor[nc].is_on_boundary[east] = true;
					elem_neighbor[nc].boundary_index[east] = nb;
					boundary_elem_ID[nb].element_index = nc;
					boundary_elem_ID[nb].side = east;
					edges[nb].nodes[0] = elements[nc].nodes[1];
					edges[nb].nodes[1] = elements[nc].nodes[2];
					if (Lnod_in == 2) edges[nb].nodes[2] = elements[nc].nodes[5];
					else if (Lnod_in == 3) {
						edges[nb].nodes[2] = elements[nc].nodes[6];
						edges[nb].nodes[3] = elements[nc].nodes[7];
					}
					nb++;
				}
			}
			else elem_neighbor[nc].neighbor[east] = nc + 1; //index of the neighboring element in the South

			nc++;
			//if (std::any_of(std::begin(elem_neighbor[nc].is_on_boundary), std::end(elem_neighbor[nc].is_on_boundary), [](bool i) {return i; } ) ) nb++; //update the
		}  //for i < N_el_i , for j < N_el_j

	N_edges_boundary = nb; //the corrected number of edges on the boundary (it exclude the edge woth periodic BC)
	return 1;
}

