# Hydroflex
Temporary repository to relativistic SPH code

## Author: José Hugo Elsas
## Colaborators: Rafael Derradi, Takeshi Kodama, Tomoi Koide

Related publications: 
   - Gaspar Elsas, J. H. , Koide, T. and Kodama, T. . Noether’s Theorem of Relativistic-Electromagnetic
Ideal Hydrodynamics. Brazilian Journal of Physics (Printed), v. 45, p. 334-339, 2015.

# Contents: 

The objective of this code is to provide a Smoothed Particle Hydrodynamics code for Ideal Relativistic Hydrodynamics, based in the paper by Aguiar et. al. [1]. The original objective was to integrate an out-of-equilibrium dust-like component to a main SPH hydro-solver, to be able to solve thermalization+hydrodynamics together in the time evolution. 

The main part of the code is the SPH solver, that is split along its components.

- hydro-mcv: main executable file
- freezeout_extern : freeze-out post processing executable file
- SPH-ODEsolver: main time evolution solver of the code, contains most of the critical time evolution and physics components
- SPH-lista: fast linked-list pair search algorithm implementation, which enables the efficient execution of the hydrodynamics code. It changes the two-point force computation from a always $O(N^2)$ algorithm to an average $O(N \log N)$, worst case $O(N^2)$ algorithm, which can run much faster and enable reasonable cases to be studied.
- SPH-state: Equation of State (EoS) and Freeze-Out implementations. These include EoS for pure gas (Landau-Type), simple QGP and dust EoS, besides a constant time Freeze-Out computation. 
- SPHtypes : Implements the main datatypes used in the code
- kernel: Smoothing Kernels definition file
- utilities: Some extra pratical utilities for the code

[1]: C E Aguiar, T Kodama, T Osada & Y Hama (2001). Smoothed particle hydrodynamics for relativistic heavy-ion collisions. Journal of Physics G: Nuclear and Particle Physics, 27, 75.

Obs:

- In case makefile don't work, try creating a obj folder on the home directory. It's necessary to store the .o files temporarely
