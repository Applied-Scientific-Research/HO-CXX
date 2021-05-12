/*
 * main.cpp - Driver code for HO-CXX, a high-order solver in vorticity variables
 *
 * (c)2018-21 Applied Scientific Research, Inc.
 *            Mohammad Haji <mhajit@gmail.com>
 *            Based on work by
 *            Adrin Gharakhani <adrin@applied-scientific.com>
 *            Christoper Stone <cpstone@gmail.com>
 */

#include "HO_2D.hpp"
#include "Mesh.hpp"

#include <iostream>


int main()
{
    const std::string input_file_name = "input.dat"; //the file listing all the settings
    HO_2D ho_2d;
    ho_2d.read_input_file(input_file_name); //true means the OCC format. zero means the buildin format of GMSH
    int success_mesh = ho_2d.setup_mesh(); //if problem type is 10 then read from file, otherwise form a specific predefined mesh based on grid_type
    char success_allocate = ho_2d.allocate_arrays();
    ho_2d.setup_sps_gps(); //set the local coordinates of the solution points, geometry points
    ho_2d.read_process_sample_points();
    ho_2d.form_bases(); //form the sps lagrangian shepe function (&derivatives), gps shape functions (&derivatives), form Radau bases
    ho_2d.form_metrics();
    ho_2d.setup_IC_BC_SRC();
    ho_2d.form_Laplace_operator_matrix(); //form the LHS matrix in Laplace discretization in Eigen/AMGCL/HYPRE format. Done only once. (The type of BC comes from BC_switch_psi, so the type of BC for psi can not change over time
    ho_2d.solve_vorticity_streamfunction(); //solves the vorticity and stream functions equations for all the times
    ho_2d.release_memory();
}

