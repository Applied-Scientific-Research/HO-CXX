#include <iostream>
#include "calculation.h"
#include "preprocess.h"

int main()
{
    const std::string input_file_name = "D:/ASR/HO_CPP_2D_Mohammad/input.dat"; //the file listing all the settings
    HO_2D ho_2d;
    ho_2d.read_input_file(input_file_name);
    int success_mesh = ho_2d.setup_mesh(); //if problem type is 10 then read from file, otherwise form a specific predefined mesh based on grid_type
    char success_allocate = ho_2d.allocate_arrays();
    ho_2d.setup_sps_gps(); //set the local coordinates of the solution points, geometry points
    ho_2d.form_bases(); //form the sps lagrangian shepe function (&derivatives), gps shape functions (&derivatives), form Radau bases
    ho_2d.setup_IC_BC_SRC();
    ho_2d.form_metrics();
    


}
