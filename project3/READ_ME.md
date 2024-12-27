The MD_code is a Fortran program that simulates the molecular dynamics (MD) of a system of atoms. 
The program uses the classical Lennard-Jones potential to model interactions between atoms and integrates the system's motion over time 
using the Verlet integration algorithm. The simulation tracks the motion of atoms in 3D space, 
computes potential and kinetic energies, and writes the atomic coordinates to an output file at each time step.

The code expects an input file that contains the number of atoms, as well as their initial coordinates and masses,
and simulates the dynamics for 1000 steps.

Strucure of the code:
You have to introduce the name of the input file having this structure:
 1) First line: show the number of atoms
 2) The coordinates of the atoms in 3D space (x, y, z) in amstrons  and their masses

Then the program will automatically compute the molecular dynamics and give a trajectories file with the coordinates of the atoms at each time step.

The program includes several functions and subroutines designed to calculate the various parameters required for the molecular dynamics simulation.

First, it computes the distances between each pair of atoms, followed by calculating the system's potential energy using the Lennard-Jones potential. It then determines the accelerations and updates the velocities and the next coordiantes. These calculations are integrated into the Verlet algorithm to simulate the motion of atoms over time.

Subrutines and functions:

read_Natoms :  read the number of atoms of the system and their coordinates and masses 
compute_distances :  Calculates the pairwise distances between atoms.
V: Computes the Lennard-Jones potential energy for the system.
T: Computes the total kinetic energy based on atom velocities.
total_energy: Computes the total energy (sum of potential and kinetic energy).
compute_acc :  Computes the accelerations based on the interatomic forces using the Lennard-Jones potential.
verlet_algorithm : performs the Verlet integration to iteratively update atomic positions and velocities. At each time step, computes and records the system's potential energy, kinetic energy, and total energy. Additionally, generates an output file containing the molecular dynamics trajectory, including atomic coordinates and energy values for visualization and analysis.
itoa :  A helper function to convert an integer to a string (used in the trajectory file).
 
  

Key Features:

Molecular Dynamics Simulation: Uses the Lennard-Jones potential to calculate interatomic forces and accelerations.

Verlet Integration: Utilizes the Verlet algorithm to update atom positions and velocities over time.

Energy Calculation: Computes and displays the potential energy, kinetic energy, and total energy at each step.

Trajectory Output: Saves the atom positions in an XYZ file format for visualization 

Dynamic Memory Allocation: Dynamically allocates memory for various data structures based on the number of atoms.