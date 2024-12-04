# TuringPatterns1DHybridGrowth
Code to support the paper ``Turing patterns under hybrid apical-uniform domain growth"<br>
by William A. Waters and Robert A. Van Gorder<br>
Department of Mathematics and Statistics, University of Otago, P.O. Box 56, Dunedin 9054, New Zealand<br>
robert.vangorder@otago.ac.nz

We briefly discuss codes used to run simulations shown in the paper. In carrying out simulations shown in the paper, we run simulations of the fully nonlinear reaction-diffusion systems using the finite difference method. We simulate these under relevant the static-domain formulations and map the solutions back to the time-varying domains to generate plots. We also simulate the system of ODEs obtained through the Turing instability analysis for each of the various Theorems and Corollaries in the paper, where relevant. All of these are implemented via MATLAB. 

There are several specific m files, and their brief descriptions are:<br>
Uniform_Solver.m - Code used to solve the uniform growth problem<br>
uniform_condition_calc.m - Code to generate plots of instability regions for the uniform growth problem<br>
Apical_Solver.m - Code used to solve the apical growth problem and associated linear instability problem<br>
Hybrid_Solver.m - Code used to solve the one--compartment hybrid growth model and associated linear instability problem<br>
Hybrid_script.m - Code used to input parameters for multi--compartment hybrid growth<br>
Hybrid_Solver_Function.m - Code used to solve the multi--compartment hybrid growth model and associated linear instability problem<br>
extract_gen_fourier_coff.m - Code called on to compute Fourier coefficients

The most involved code is that used to solve the multi-compartment form of hybrid growth described in Section 5 of the paper. To generate the simulations shown in Figures 6, 7, and 8 we run Hybrid_script.m, with the parameters to be entered including:<br>
node_initals: Row vector of apical growth node locations as t=0 (include 0 and 1)<br>
s_dots: Cell array of apical growth functions at each node<br>
S_ints: Cell array of uniform growth functions in each compartment<br>
final_time: Run the simulations until this time<br>
time_res: How many timesteps used in PDE simulation<br>
space_res: How many spacesteps used in PDE simulation<br>
max_k: Where to truncate the series in the ODE system 

As an example, to generate the plots in Figure 6(a), we enter the parameters:<br>
 time_res=1001;<br>
 space_res=201;<br>
 node_initals = [0,1/2,1];<br>
 s_dots = {@(t) 0*t+0.04, @(t) 0*t+0.04, @(t) 0*t+0.04};<br>
 S_ints = {@(t) 0*t, @(t) 0*t};<br>
 max_k=60;<br>
 final_time=100;<br>
We choose 1001 timesteps and 201 space steps at each timestep, with three apical nodes initially at $\xi$ = 0, 0.5, and 1. s_dots and S_ints describe the evolution of the apical nodes and uniform regions, as discussed in the paper. The maximal number of modes displayed in the instability analysis is 60, and the simulation is run out to a final time of $t = 100$.

The parameters are then used in a function Hybrid_Solver_Function.m, taking in the parameters above like<br>
Hybrid_Solver_Function(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)<br>
and this function which actually carries out the simulations. Lines 7 to 27 read in parameters. Lines 30-56 solve the reaction-diffusion system (using Schnakenberg kinetics) on a static domain (in $\xi$) subject to random initial conditions and no-flux boundary conditions, as indicated in the paper. The solution method is an implementation of the method of lines and finite difference method in space using pdepe. Code used to plot the PDE solution is given on lines 58 to 83 of the code. We then extract the dominant mode in the PDE solution, for later comparison with the linear instability analysis, using code on lines 86 through 133. The rest of the code solves the ODE system obtained from the linear instability analysis of the spatially heterogeneous linearised problem, and this is shown on lines 136 through 250. The ODE system has terms which undergo exponential growth and we find that the ODE solver ode15s is effective. The rest of the function, lines 253 through 338, contained useful sub-functions and quantities. 

Similar yet simpler m files were used to obtain solutions shown in other figures in the paper. Uniform_Solver.m was used for purely uniform growth of Section 3(a) as shown in Figure 1 while Apical_Solver.m was used for the purely apical growth of Section 3(b) as shown in Figure 2. The single-compartment form of hybrid growth, with two apical nodes bounding a single uniform growth region, was simulated using Hybrid_Solver.m, with the corresponding results shown in Figure 4. These files still use pdepe to solve the relevant reaction diffusion systems, yet ode45 is used for any ODE systems in the tests for the Turing instability, as the relevant ODEs are simpler than in the multi-compartment hybrid case (where ode15s results in faster run time). The instability plots are generated in a similar manner for all cases except for the uniform case as there the instability results are purely analytic under Corollary 3.1, and utilise the implementation in uniform_condition_calc.m rather than the simulation of a system of ODEs.

PDE solutions for the relevant reaction-diffusion systems are plotted for the activator, u, in all cases. The heat map runs from blue (smallest or zero value) and yellow (maximal value). In Figures 1, 2, and 4 we plot mode amplitudes of the PDE solutions calculated from the numerical simulations, as these are useful for detecting dominant spatial modes at each timestep. (We obtain this information in the simulations of Section 5, but do not use the space to plot it there, rather just extracting the dominant mode.) In Figures 2, 4, 6, 7, and 8 we plot the amplitudes of the linear modes present in the perturbation of a base state near the Turing instability at each point in time, which allows us to understand which spatial modes contribute to the spatial instability leading to pattern formation. As expected from the theory, the envelope of unstable modes changes over time as the domain evolves. We plot the dominant spatial mode from the FDM simulation of the PDE in red against these unstable modes, for sake of comparison.

All code principally written by W. A. Waters and modified in minor ways by R. A. Van Gorder.
