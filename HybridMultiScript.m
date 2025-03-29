%% Script to run using Hybrid_Solver_Function(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k):
%% Parameters are specified for multi-compartment hybrid growth model

% Attatch appropriate suffix to function to specify which Kinetics to use
% _S: Schnakenberg
% _GM: Gierer–Meinhardt
% _FHN: FitzHugh-Nagumo

% Parameter Definitions
% node_initals: Row vector of apical growth node locations at time t=0 (include 0 and 1 as the initial boundaries for the scaled domain)
% s_dots: Cell array of apical growth functions at each node
% S_ints: Cell array of uniform growth functions in each compartment
% final_time: Run the simulations until this time
% time_res: How many time-steps used in PDE simulation (ODE automatic)
% space_res: How many space_steps used in PDE simulation (ODE automatic)
% max_k: Where to trunate the series in the ODE system 


%% Also used to run Uniform_Solver(strain_rate, final_time, time_res, space_res, max_k):

% Attatch appropriate suffix to function to specify which Kinetics to use
% _S: Schnakenberg
% _GM: Gierer–Meinhardt
% _FHN: FitzHugh-Nagumo

% Parameter Definitions
% strain_rate: Uniform growth function
% final_time: Run the simulations until this time
% time_res: How many time-steps used in PDE simulation (ODE automatic)
% space_res: How many space_steps used in PDE simulation (ODE automatic)
% max_k: Where to trunate the series in the ODE system 

%% Code examples for Hybrid Solver

% We provide all examples shown in the paper in PARAMETERS.txt
% Other examples are provided here
% Uncomment parameters for a specific example to run it

%% Hybrid Testing cases (low res for speed)

% Pure uniform - switch solvers to ode115 before running
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0*t, @(t) 0*t};
% S_ints = {@(t) 0*t+0.022};
% max_k=40;
% final_time=100;

% % Left apical - use ode15s
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0.1+0*t, @(t) 0*t};
% S_ints = {@(t) 0*t};
% max_k=40;
% final_time=100;

% % Right apical - use ode15s
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0*t, @(t) 0.1+0*t};
% S_ints = {@(t) 0*t};
% max_k=40;
% final_time=100;

% % Both apical - use ode15s
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0.05+0*t, @(t) 0.05+0*t};
% S_ints = {@(t) 0*t};
% max_k=40;
% final_time=100;

% % Apical plus uniform growth - use ode15s
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0.01+0*t, @(t) 0.07+0*t};
% S_ints = {@(t) 0.01+0*t};
% max_k=50;
% final_time=100;

% % Apical plus uniform decay - use ode15s
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0.1+0*t, @(t) 0.1+0*t};
% S_ints = {@(t) -0.02+0*t};
% max_k=50;
% final_time=100;

% % 7(a) from paper - use ode15s
% time_res=201;
% space_res=101;
% node_initals = [0,1];
% s_dots = {@(t) 0.05+0*t, @(t) 0.05+0*t};
% S_ints = {@(t) -0.02*(t>50)};
% max_k=45;
% final_time=100;

%% Choose kinetics and run simulations for set parameters

% tic
% Hybrid_Solver_Function_FHN(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)
% toc
% tic
% Hybrid_Solver_Function_GM(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)
% toc
% tic
% Hybrid_Solver_Function_S(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)
% toc

%% Testing cases for Uniform Solver

% Constant strain
% time_res=201;
% space_res=101;
% strain = @(t) 0.02+0*t;
% max_k=45;
% final_time=100;

% Positive but decreasing strain
% time_res=201;
% space_res=101;
% strain = @(t) 0.04-0.0004*t;
% max_k=45;
% final_time=100;

% Positive then negitive strain
% time_res=201;
% space_res=101;
% strain = @(t) 0.02*(t<=100) - 0.02*(t>100) ;
% max_k=40;
% final_time=200;

% Oscilating strain
% time_res=201;
% space_res=101;
% strain = @(t) 0.02+0.02*sin(t/5);
% max_k=50;
% final_time=100;

%% Choose kinetics and run uniform simulations for set parameters

% Uniform_Solver_S(strain, final_time, time_res, space_res, max_k)
% Uniform_Solver_GM(strain, final_time, time_res, space_res, max_k)
% Uniform_Solver_FHN(strain, final_time, time_res, space_res, max_k)

















