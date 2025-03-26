%% Script to run using Hybrid_Solver_Function(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k):
%% Parameters are specified for multi-compartment hybrid growth model

% Parameter Definitions
% node_initals: Row vector of apical growth node locations at time t=0 (include 0 and 1 as the initial boundaries for the scaled domain)
% s_dots: Cell array of apical growth functions at each node
% S_ints: Cell array of uniform growth functions in each compartment
% final_time: Run the simulations until this time
% time_res: How many time-steps used in PDE simulation (ODE automatic)
% space_res: How many space_steps used in PDE simulation (ODE automatic)
% max_k: Where to trunate the series in the ODE system 

% We provide all examples shown in the paper
% Uncomment parameters for a specific example to run it

% Apical and mixed growth cases shown in Figure 6 of the paper

% Fig 6a
% time_res=1001;
% space_res=201;
% node_initals = [0,1/2,1];
% s_dots = {@(t) 0*t+0.04, @(t) 0*t+0.04, @(t) 0*t+0.04};
% S_ints = {@(t) 0*t, @(t) 0*t};
% max_k=60;
% final_time=100;

% Fig 6b
% time_res=1001;
% space_res=201;
% node_initals = [0,1/3,2/3,1];
% s_dots = {@(t) 0*t+0.01, @(t) 0*t+0.03, @(t) 0*t+0.05,@(t) 0*t+0.07};
% S_ints = {@(t) 0*t, @(t) 0*t, @(t) 0*t};
% max_k=80;
% final_time=70;

% Fig 6c
% time_res=1001;
% space_res=201;
% node_initals = [0,1/5,2/5,3/5,4/5,1];
% s_dots = {@(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01};
% S_ints = {@(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t};
% max_k=80;
% final_time=200;

% Fig 6d
% time_res=1001;
% space_res=201;
% node_initals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
% s_dots = {@(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01};
% S_ints = {@(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t};
% max_k=100;
% final_time=200;

% Fig 6e
% time_res=1001;
% space_res=201;
% node_initals = [0,1/2,1];
% s_dots = {@(t) 0*t, @(t) 0*t, @(t) 0*t};
% S_ints = {@(t) 0*t+0.02, @(t) 0*t+0.03};
% max_k=60;
% final_time=100;

% Fig 6f
% time_res=1001;
% space_res=201;
% node_initals = [0,1/25,24/25,1];
% s_dots = {@(t) 0*t, @(t) 0*t, @(t) 0*t, @(t) 0*t};
% S_ints = {@(t) 0*t+0.04, @(t) 0*t+0.01, @(t) 0*t+0.04};
% max_k=80;
% final_time=120;

% Fig 6g
% time_res=1001;
% space_res=201;
% node_initals = [0,1/5,2/5,3/5,4/5,1];
% s_dots = {@(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01};
% S_ints = {@(t) 0*t+0.02, @(t) 0*t+0.01, @(t) 0*t+0.02, @(t) 0*t+0.01, @(t) 0*t+0.02};
% max_k=80;
% final_time=100;

% Fig 6h
% time_res=1001;
% space_res=201;
% node_initals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
% s_dots = {@(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01, @(t) 0*t+0.01};
% S_ints = {@(t) 0*t+0.02, @(t) 0*t+0.01, @(t) 0*t+0.02, @(t) 0*t+0.01, @(t) 0*t+0.02, @(t) 0*t+0.01, @(t) 0*t+0.02, @(t) 0*t+0.01, @(t) 0*t+0.02, @(t) 0*t+0.01};
% max_k=60;
% final_time=60;


% Mixed growth and decay cases shown in Figure 7 of the paper

% Fig 7a
% time_res=1001;
% space_res=201;
% node_initals = [0,1/2,1];
% apfun = @(t) (0*t + 0.05);
% unifun = @(t) (0*t - 0.2*(t>50));
% s_dots = {apfun, apfun, apfun};
% S_ints = {unifun, unifun};
% max_k=40;
% final_time=70;

% Fig 7b
% time_res=1001;
% space_res=201;
% node_initals = [0,1/3,2/3,1];
% apfun = @(t) (0*t + 0.05);
% unifun = @(t) (0*t - 0.2*(t>50));
% s_dots = {apfun, apfun, apfun, apfun};
% S_ints = {unifun, unifun, unifun};
% max_k=40;
% final_time=70;

% Fig 7c
% time_res=1001;
% space_res=201;
% node_initals = [0,1/2,1];
% apfun = @(t) (0*t + 0.05*(t<50) - 0.2*(t>50));
% unifun = @(t) (0*t + 0.05*(t>50));
% s_dots = {apfun, apfun, apfun};
% S_ints = {unifun, unifun};
% max_k=40;
% final_time=80;

% Fig 7d
% time_res=1001;
% space_res=201;
% node_initals = [0,1/3,2/3,1];
% apfun = @(t) (0*t + 0.05*(t<50) - 0.2*(t>50));
% unifun = @(t) (0*t + 0.05*(t>50));
% s_dots = {apfun, apfun, apfun, apfun};
% S_ints = {unifun, unifun, unifun};
% max_k=40;
% final_time=80;


% Transport cases shown in Figure 8 of the paper

% Fig 8a
time_res=2001;
space_res=201;
node_initals = [0,1/2,1];
s_dots = {@(t) 0*t+0.1, @(t) 0*t-0.1, @(t) 0*t+0.1};
S_ints = {@(t) 0*t+0.02, @(t) 0*t+0.02};
max_k=100;
final_time=160;

% Fig 8b
% time_res=2001;
% space_res=201;
% node_initals = [0,1/2,1];
% s_dots = {@(t) 0*t-0.1, @(t) 0*t+0.1, @(t) 0*t-0.1};
% S_ints = {@(t) 0*t+0.02, @(t) 0*t+0.02};
% max_k=100;
% final_time=160;

% Fig 8c
% time_res=2001;
% space_res=201;
% node_initals = [0,1/2,1];
% s_dots = {@(t) 0*t+0.1, @(t) 0*t+0.1, @(t) 0*t+0.1};
% S_ints = {@(t) 0*t-0.02, @(t) 0*t-0.02};
% max_k=100;
% final_time=200;

% Fig 8d
% time_res=1001;
% space_res=201;
% node_initals = [0,1/2,1];
% s_dots = {@(t) 0*t-0.01, @(t) 0*t-0.01, @(t) 0*t-0.01};
% S_ints = {@(t) 0*t+0.040, @(t) 0*t+0.040};
% max_k=100;
% final_time=200;


% time_res=501;
% space_res=101;
% node_initals = [0,0.5,1];
% 
% s_dots = {@(t) 0.05*(t<50)-0.15*(t>=50), @(t) 0.05*(t<50)-0.15*(t>=50), @(t) 0.05*(t<50)-0.15*(t>=50)};
% S_ints = {@(t) 0.05*(t>=50),@(t) 0.05*(t>=50)};
% max_k=35;
% final_time=72;


Hybrid_Solver_Function_S(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)

Hybrid_Solver_Function_GM(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)

Hybrid_Solver_Function_FHN(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)
















