# HO-Fortran
Current research version of the HO solver

First "successful" attempt using flow in a cavity

Use the folowing to compile Test_Package

g++  -O2 -fopenmp -I./include -lgfortran aplles_interface_fortran.F90 aplles_interface.cxx VTM.f90 -L./libs -laplles_mp -lxml2 -llapack

The code ran successfully for Re = 100, Knod = 2, dt = 1.e-5, and NX = NY = 128. However, all other simulations using higher Re and/or Knod
eventually blow up even at low mesh resolution and smaller dt. The culprit seems to be convection, which goes unstable. More work is 
required to establish the proper stability criteria, as well as double check for potential bugs

I've put everything in VTM.f90 for now, cuz the file is small enough to work with just one file

NOTE: Only first-order time integration is in place for now!
