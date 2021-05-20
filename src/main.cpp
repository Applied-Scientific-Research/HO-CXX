/*
 * main.cpp - Driver code for HO-CXX, a high-order solver in vorticity variables
 *
 * (c)2018-21 Applied Scientific Research, Inc.
 *            Mohammad Haji <mhajit@gmail.com>
 *            Mark J Stock <markjstock@gmail.com>
 *            Based on work by
 *            Adrin Gharakhani <adrin@applied-scientific.com>
 *            Christoper Stone <cpstone@gmail.com>
 */

#include "HO_2D.hpp"

#include <iostream>

int main(int argc, char const *argv[])
{
    //file listing all the settings, default name
    std::string input_file_name = "input.dat";

    // grab the mesh file from the command line
    if (argc == 2) {
        input_file_name = argv[1];
    } else if (argc > 2) {
        std::cout << std::endl << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << " [input.dat]" << std::endl << std::endl;
        return -1;
    }

    // the kitchen sink class
    HO_2D ho_2d;

    //true means gmsh file made with opencascade kernel, zero means the built-in format of GMSH
    ho_2d.read_input_file(input_file_name);

    //if problem type is 10 then read from file, otherwise form a specific predefined mesh based on grid_type
    auto success_mesh = ho_2d.setup_mesh();
	if (success_mesh != 1) {
		std::cout << std::endl << "Mesh file (" << input_file_name << ") not processed correctly!" << std::endl;
        return success_mesh;
    }
    auto success_allocate = ho_2d.allocate_arrays();
	if (success_allocate != 1) {
		std::cout << std::endl << "Arrays not allocated correctly!" << std::endl;
        return success_allocate;
    }

    //set the local coordinates of the solution points, geometry points
    ho_2d.setup_sps_gps();
    ho_2d.read_process_sample_points();

    //form the sps lagrangian shape function (&derivatives), gps shape functions (&derivatives), form Radau bases
    ho_2d.form_bases();
    ho_2d.form_metrics();
    ho_2d.setup_IC_BC_SRC();

    //form the LHS matrix in Laplace discretization in Eigen/AMGCL/HYPRE format. Done only once.
    // (The type of BC comes from BC_switch_psi, so the type of BC for psi can not change over time
    ho_2d.form_Laplace_operator_matrix();

    //solves the vorticity and stream functions equations for all the times
    ho_2d.solve_vorticity_streamfunction();

    // dealloc and quit
    ho_2d.release_memory();
}

