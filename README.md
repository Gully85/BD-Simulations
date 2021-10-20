BD-Simulation

This is a C++ code that performs a many-particle simulation in a periodic simulation box.

To run: You need all .cpp and .h files and the Makefile. Use make sim. The gcc compiler must be installed.
The analysis part (make obs, make snapshots, and combinations of these. See Makefile) also need the .py and .plt files, Python 3 and Gnuplot must be installed for it.

The most important file for the user is parameter.h, all parameters are set there. See comments in that file. The code will execute timesteps in a Langevin time integration scheme. dynamik_methoden.cpp defines the calculation of forces between all particle types.

Output of the sim are files pos1_<x>.txt and pos2_<x>.txt with x from 0 to (runs-1), runs being defined in parameter.h. This means the positions of particles, species 1/2, in run x. In the file, there time and x/y coordinate for each particle.