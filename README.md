# HO-CXX
Research version of C++ High-Order solver for vorticity transport in 2D on unstructured meshes

## Build and run
You should have AMGCL (which needs boost-devel) or Eigen installed. Use the following commands to compile:

    mkdir build
    cd build
    ccmake ..
    make

To run a sample case, perform the following from the `build` directory:

    cp ../input/* ./
    ./HO-CXX.bin

The executable will always read the `input.dat` file in the current working directory. Make sure that the appropriate `*.msh` file is also in the cwd. Output will be to vtk-legacy-format files every 1000 steps by default.

This version makes 3 distinct cases for Cavity, BFS and flow around cylinder. See the `input.dat` file for details, or use the `input` directory for meshes and input files.

## Credits
This project is funded by the [National Institutes of Health (NIH)](https://www.nih.gov/) under grant number 1 R01 EB022180-01A1 ("A Fast High-Order CFD for Turbulent Flow Simulation in Cardio-Devices").

It makes use of [Eigen](http://eigen.tuxfamily.org/) and [amgcl](https://github.com/ddemidov/amgcl).

This code is mostly by Mohammad Hajit, based on earlier work of Adrin Gharakhani, both
of Applied Scientific Research, Inc. Refactoring and hybridization work by Mark Stock, also at ASR.

