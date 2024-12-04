# TuringPatterns1DHybridGrowth
Code to support the paper ``Turing patterns under hybrid apical-uniform domain growth"<br>
by William A. Waters and Robert A. Van Gorder<br>
Department of Mathematics and Statistics, University of Otago, P.O. Box 56, Dunedin 9054, New Zealand<br>
robert.vangorder@otago.ac.nz

We briefly discuss specific MATLAB code used to run the simulations shown in the paper. In carrying out simulations shown in the paper, we run simulations of the fully nonlinear reaction-diffusion systems using the finite difference method. We simulate these under relevant the static-domain formulations and map the solutions back to the time-varying domains to generate plots. We also simulate the system of ODEs obtained through the Turing instability analysis for each of the various Theorems and Corollaries in the paper, where relevant. All of these are implemented via MATLAB. 

There are several specific m files, and their brief descriptions are:<br>
Uniform_Solver.m - Code used to solve the uniform growth problem<br>
uniform_condition_calc.m - Code to generate plots of instability regions for the uniform growth problem<br>
Apical_Solver.m - Code used to solve the apical growth problem and associated linear instability problem<br>
Hybrid_Solver.m - Code used to solve the one--compartment hybrid growth model and associated linear instability problem<br>
Hybrid_script.m - Code used to input parameters for multi--compartment hybrid growth<br>
Hybrid_Solver_Function.m - Code used to solve the multi--compartment hybrid growth model and associated linear instability problem<br>
extract_gen_fourier_coff.m - Code called on to compute Fourier coefficients

More specifically...

Uniform_Solver.m solves the reaction-diffusion equation under Schnackenberg reaction kinetics and diffusion parameters and on a uniformly evolving domain. The problem is mapped to a stationary domain and is then solved using the method of lines and a finite difference discretisation on the stationary space domain using the pdepe routine. Coefficients of the spatial modes corresponding to this numerical simulation are calculated and the dominant mode at each timestep is recorded. The linear instability result from Corollary 3.1 is obtained (using uniform_condition_calc.m) and compared with the spatial frequency of the dominant mode from the full PDE simulations. 

Apical_Solver.m solves the reaction-diffusion equation under Schnackenberg reaction kinetics and diffusion parameters and on an apically evolving domain. The problem is mapped to a stationary domain and is then solved using the method of lines and a finite difference discretisation on the stationary space domain using the pdepe routine. Coefficients of the spatial modes corresponding to this numerical simulation are calculated and the dominant mode at each timestep is recorded. The linear instability result from Corollary 3.2 is obtained by solving the relevant ODEs for the linear mode coefficients (using ode45) and compared these with the spatial frequency of the dominant mode from the full PDE simulations. 

Hybrid_Solver.m solves the reaction-diffusion equation under Schnackenberg reaction kinetics and diffusion parameters and on a space domain with the simplest hybrid growth. Here there is one compartment, with apical nodes at each boundary enclosing a region of uniform evolution. The problem is mapped to a stationary domain and is then solved using the method of lines and a finite difference discretisation on the stationary space domain using the pdepe routine. Coefficients of the spatial modes corresponding to this numerical simulation are calculated and the dominant mode at each timestep is recorded. The linear instability result from Corollary 4.1 is obtained by solving the relevant ODEs for the linear mode coefficients (using ode45) and compared these with the spatial frequency of the dominant mode from the full PDE simulations. 

Hybrid_Solver_Function.m solves the reaction-diffusion equation under Schnackenberg reaction kinetics and diffusion parameters and on a space domain with multi-compartment hybrid growth. There are multiple compartments, with apical nodes at each external and interior boundary enclosing a region of uniform evolution. The user specifies the positions of apical nodes, and hence the number of compartments. The multi-compartment hybrid evolvution problem is mapped to a stationary domain and is then solved using the method of lines and a finite difference discretisation on the stationary space domain using the pdepe routine. Coefficients of the spatial modes corresponding to this numerical simulation are calculated and the dominant mode at each timestep is recorded. The linear instability result from Corollary 5.1 is obtained by solving the relevant ODEs for the linear mode coefficients (using ode15s) and compared these with the spatial frequency of the dominant mode from the full PDE simulations. We use ode15s rather than ode45 for the linear instability ODEs as we find the multi-compartment case particularly involved to simulate, and ode15s speeds up the computation of the stiff ODE solutions nicely for this more involved problem. 

All parameters (physical and numerical) for the multi-compartment hybrid growth model are entered into the script Hybrid_script.m which then calls Hybrid_Solver_Function.m to perform the calculations. Parameters to be entered include:<br>
node_initals: Row vector of apical growth node locations as t=0 (include 0 and 1)<br>
s_dots: Cell array of apical growth functions at each node<br>
S_ints: Cell array of uniform growth functions in each compartment<br>
final_time: Run the simulations until this time<br>
time_res: How many timesteps used in PDE simulation<br>
space_res: How many spacesteps used in PDE simulation<br>
max_k: Where to truncate the series in the ODE system <br>

We can vary the space- and timesteps to test the convergence of the solutions for each problem, and we have done this in several cases for various domain evolutions. The solutions appear to converge to what one suspects is a true solution for the reaction-diffusion system in each case. 

PDE solutions for the relevant reaction-diffusion systems are plotted for the activator, u, in all cases. The heat map runs from blue (smallest or zero value) and yellow (maximal value). In Figures 1, 2, and 4 we plot mode amplitudes of the PDE solutions calculated from the numerical simulations, as these are useful for detecting dominant spatial modes at each timestep. (We obtain this information in the simulations of Section 5 and Figures 6, 7, and 8, yet do not use the space to plot it there, rather just extracting the dominant mode.) In Figures 2, 4, 6, 7, and 8 we plot the amplitudes of the linear modes present in the perturbation of a base state near the Turing instability at each point in time, which allows us to understand which spatial modes contribute to the spatial instability leading to pattern formation. As expected from the theory, the envelope of unstable modes changes over time as the domain evolves. We plot the dominant spatial mode from the simulation of the PDE systems in red against these unstable modes, for sake of comparison.

All MATLAB code principally written by W. A. Waters and modified in minor ways by R. A. Van Gorder.
