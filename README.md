# Dendrimers
Contain files related to Dendrimer project 

# Analysis
Done using VMD after obtaining trajectories by LAMMPS 

# Displacement.tcl
This tcl script can find average displacement of selected molecules (here phenol) in Angstom per seconds within user defined radius (20 here). 
Data will be written to the displacements_within20ofDen_COM.dat file.

# MSD.tcl
This script can be use to find Diffusion coefficient of one type of selected molecule within user defined radius. Here we calculated diffusion coefficient of phenol beads within 20 Angstrom of radius from center of mass of reference molecule. Reference molecule is "Segname D1" and phenols are "resname PHNL". Replace those two selections as your need.
Data will be written to the MSD.dat file.
