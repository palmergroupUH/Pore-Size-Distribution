# Pore-Size-Distribution
Computes pore size distribution using Torquato and Gubbin's pore size method.
This code assumes porous medium made of uniform particles, and uses particle trajectories in dcd format.
It takes particle and probe size as inputs
It uses Linked-list method to compute minimum distance of a point in the void space from nearby porous medium particles
It uses solvopt non-linear optimization code to compute Gubbin's pore size
