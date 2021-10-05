/*
 * main.cpp - Driver code for HO-CXX, a high-order solver in vorticity variables
 *
 * (c)2018-21 Applied Scientific Research, Inc.
 *            Mohammad Haji <mhajit@gmail.com>
 *            Based on work by
 *            Adrin Gharakhani <adrin@applied-scientific.com>
 *            Christoper Stone <cpstone@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "HO_2D.hpp"

#include <iostream>


int main(int argc, char const *argv[])
{
	// file listing all the settings, default name
	std::string input_file_name = "input.dat";

	// grab the input file from the command line
	if (argc == 2) {
		input_file_name = argv[1];
	} else if (argc > 2) {
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "  " << argv[0] << " [input.dat]" << std::endl << std::endl;
		return -1;
	}

	// the kitchen sink class
	HO_2D ho_solver;

	// return true means gmsh file made with opencascade kernel, zero means the built-in format of GMSH
	ho_solver.read_input_file(input_file_name);

	auto error_mesh = ho_solver.setup_mesh();
	if (error_mesh) {
		std::cout << std::endl << "Mesh file (" << input_file_name << ") not processed correctly!" << std::endl;
		return error_mesh;
	}

	auto error_allocate = ho_solver.allocate_arrays();
	if (error_allocate) {
		std::cout << std::endl << "Arrays not allocated correctly!" << std::endl;
		return error_allocate;
	}

	// set the local coordinates of the solution points, geometry points
	ho_solver.setup_sps_gps();

	// read and process the sampling points
	ho_solver.read_process_sample_points();

	// form the sps lagrangian shape function (&derivatives), gps shape functions (&derivatives), form Radau bases
	ho_solver.form_bases();
	ho_solver.form_metrics();
	ho_solver.setup_IC_BC_SRC();

	// form the LHS matrix in Laplace discretization in Eigen/AMGCL/HYPRE format. Done only once.
	// (The type of BC comes from BC_switch_psi, so the type of BC for psi can not change over time
	ho_solver.form_Laplace_operator_matrix();

	// solves the vorticity and stream functions equations for all the times
	ho_solver.solve_vorticity_streamfunction();

	// return success
	return 0;
}
