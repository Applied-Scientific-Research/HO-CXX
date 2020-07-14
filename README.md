# HO-Fortran
Research version of High-Order solver for vorticity transport in 2D

## Build and run
Use the following commands to compile `Test_Package`:

    cd Test_Package
    make

To run a test case, edit `input.dat` and run the binary executable:

    ./2DVortTrans.bin

Output will be to legacy-format `Geometry.vtk` and `Vorticity*.vtk` files.

## Notes
The code ran successfully for Re = 100, Knod = 2, dt = 1.e-5, and NX = NY = 128. However, all other simulations using higher Re and/or Knod
eventually blow up even at low mesh resolution and smaller dt. The culprit seems to be convection, which goes unstable. More work is 
required to establish the proper stability criteria, as well as double check for potential bugs

Everything is in VTM.f90 for now, as the program is small enough to work with just one file.

Only first-order time integration is in place for now!

## Credits
The majority of this code is by Adrin Gharakhani of Applied Scientific Research, Inc. 
It makes use of the APLLES solver by Christopher Stone, and XML and Lapack libraries.
