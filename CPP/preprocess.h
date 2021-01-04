#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>


//enum edges_types {_2Node_edge = 1, _3Node_edge = 8, _4Node_edge = 26, _5Node_edge = 27, _6Node_edge = 28};
//std::vector< std::vector<unsigned int>> edge_type_node_number{ {1,8,26,27,28}, {2,3,4,5,6} }; //holds the number of nodes on each edge type
//std::vector< std::vector<unsigned int>> face_type_node_number{ {3,10,16,36,37}, {4,9,8,16,25} }; //holds the number of nodes on each 2D element type

// class for 2D coordinate system

enum cell_sides { south = 0, east = 1, north = 2, west = 3 }; //the 4 sides of a cell

struct Cmpnts2 {
	double x, y; //the components for 2D coordinate system

	Cmpnts2(double _x, double _y) : x(_x), y(_y) {};
	Cmpnts2(const Cmpnts2& cmp) :
		x(cmp.x), y(cmp.y)
	{}
	Cmpnts2() {
		x = y = 0.; //default constructor
	}

	Cmpnts2& operator=(const Cmpnts2& rhs) // compound assignment (does not need to be a member,
	{     // addition of rhs to *this takes place here
		x = rhs.x;
		y = rhs.y;
		return *this; // return the result by reference
	}

	~Cmpnts2() {}	//destructor define here later

	void add(Cmpnts2 a, Cmpnts2 b) {
		x = a.x + b.x;
		y = a.y + b.y;
	}
	void subtract(Cmpnts2 a, Cmpnts2 b) {
		x = a.x - b.x;
		y = a.y - b.y;
	}
	void multiply(double r, Cmpnts2 a) {
		x = r * a.x;
		y = r * a.y;
	}
	void scale(double r) {
		x *= r;
		y *= r;
	}
	double norm2() {
		return std::sqrt(x * x + y * y);
	}
	void set_coor(double a, double b) {
		x = a; y = b;
	}

	void plus(double b, Cmpnts2 Y) {
		x += b*Y.x; y += b*Y.y;
	}
};
struct element2d { // the 2D element
	unsigned int element_type;
	unsigned int N_nodes;
	std::vector<unsigned int> nodes;  //all nodes used in the formation of an element. The details of these nodes can be found in vector nodes
	unsigned int edges[4];  //all 4 edges used in the formation of an element. The details of these edges can be found in vector edges
};
struct edge {  //the side edges of the 2D elements
	unsigned int edge_type;
	unsigned int N_nodes;
	std::vector<unsigned int> nodes;
};
struct node {  //the nodes
	unsigned int node_type;  //0 for corner, 1 for nodes on the edges and 2 for nodes on the faces
	Cmpnts2 coor;
};
struct boundary {  //the boundary for 2D mesh
	std::string name;  //name of the boundary
	unsigned int N_edges=0; //number of the 1D edges constituting the boundary
	std::vector<unsigned int> edges;  //the index of the edges that form the boundary
	boundary() : N_edges(0) {} //default constructor
};
struct side {
	bool west, east, north, south;
};
struct boundary_element {  //to store the element index of the edges that are located on the boundary and the side of the element
	unsigned int element_index;
	cell_sides side;
};
struct element_neighbor {
	unsigned int neighbor[4]; //4 sides, south, east, north, west
	bool is_on_boundary[4] = { false,false,false,false}; //by default the 4 side are not located on the boundary
	unsigned int boundary_index[4]; //the boundary index for the element sides that are located on the boundary
	int neighbor_common_side[4]; //the side of the neighboring element that is common
};

class Mesh {
private:
	int N_nodes = 0;     //# of nodes read from msh file: total number of geometry nodes in the domain
	int N_elements = 0;  //# of elements read from msh file
	int N_edges_boundary;  // total number of edges on the boundary
	int N_Gboundary; //number of global boundaries with boundary conditions: number of boundary conditions (number of different boundaries)
	std::vector<node> nodes; //coordinates of the nodes
	std::vector<edge> edges;		 // types and nodes constituting edges of elements
	std::vector<element2d> elements; // types and nodes constituting elements
	std::vector<boundary> boundaries; // the boundaries of the problem, including the names and the edges that form each boundary
	std::string input_msh_file; //default name of the mesh file. It can be read from the file too
	double dx_ratio; //multiplier to expand the grid geometrically by factor dx_ratio
	double fac; //multiplier to make grids randomly non-uniform; uniform: fac=0.
	unsigned int Lnod; //number of geometry nodes on each edge (=Lnod_in + 1)
	unsigned int N_el_i, N_el_j; //number of quad elements in i (csi direction), j (eta direction) if the mesh is structured quad only
	unsigned int N_el;  //total number of quad elements
	unsigned int Lnod_in; //geometry parametric order per elem (must generalize to include x, y dirs). Lnod_in is the degree of the polynomial made from the geometry nodes, so the number of nodes on the element geometry is one more than Lnod_in.e.g. 5 nodes on edge is Lnod_in = 4
	element_neighbor* elem_neighbor; //index of the neighboring elements on the 4 sides and the boundary index (if located on the boundary). south neighbor is index0, east=1, north=2, west=3 (so CCW)
	boundary_element* boundary_elem_ID;  //index pointing to mesh element number containing that boundary element (size is number of edges on the boundaries)
	//boundaryPointsElemID is a 1D array that keeps the indices of the elements that have an edge on the boundary. The location of the boundary is known by elemID(1:4,nc) when it stores negative values
	//unsigned int** node_ID; //an array to hold the indices of the nodes of the elements, going from SW CCW, redundant as elements[] have the same info
	//unsigned int** boundary_node_ID; //The nodes that form each boundary element. redundant, as edges[] and boundaries[] have the information

	//enum edges_types {_2Node_edge = 1, _3Node_edge = 8, _4Node_edge = 26, _5Node_edge = 27, _6Node_edge = 28};
	std::vector< std::vector<unsigned int>> edge_type_node_number{ {1,8,26,27,28}, {2,3,4,5,6} }; //holds the number of nodes on each edge type
	std::vector< std::vector<unsigned int>> face_type_node_number{ {3,10,16,36,37}, {4,9,8,16,25} }; //holds the number of nodes on each 2D element type
	std::vector<std::vector<unsigned int>> element_edge_node_number{ {3,10,16,36,37}, {2,3,3,4,5} }; //holds the number of nodes on each edge of a 2D element type


public:
	Mesh() { //default constructor
		input_msh_file = "mesh_file.msh"; // the default name of the mesh file exported from Gmsh
		fac = 0.; //uniform
		dx_ratio = 1; //no expanding grid
		N_edges_boundary = 0;
		N_Gboundary = 4;
	} 

	~Mesh() {}	//destructor define here later

	char read_msh_file();
	char setup_mesh_problem(unsigned int problem_type);

	int locate_in_file(std::ifstream& filestream, const std::string& searchname); //to find a searchname in the MSH file and return the file stream

	void process_mesh(); //processes the mesh that is read from file, finding the elements neighbors, ...

	friend class HO_2D;

};

