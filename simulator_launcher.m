clc; clear all; close all;
%% Configuration

addpath('mesh_processing');
addpath('propagators');
addpath('stopping_criteria');
addpath('timestep_control');
addpath('problems/energy_2D_fvm');
addpath('results_manager');
addpath('utils');
addpath('constants')

nodes_file = FilePaths.NODES;
cells_file = FilePaths.CELLS;
bc_files = FilePaths.BOUNDARY_CONDITIONS;

%% Mesh processing

[mesh] = mesh_processor(nodes_file, cells_file, bc_files);

%% Problem configuration

w0 = 300 .* ones(length(mesh.volume), 1);
t0 = 0.;

cp = 1004.5; % J / (kg K)
cv = cp / 1.4; % J / (kg K)
rho = 1.225; % kg / (m ^ 3)
k = 0.024; % W / (m K)

speed = @(x, y, t) [0., 0.01];
convective_method = @(w, t, mesh, cp, cv, rho, speed) ...
    dummy_fluxes(w, t, mesh, cp, cv, rho, k, speed);
conductive_method = @(w, t, mesh, cv, rho, k) ...
    dummy_fluxes(w, t, mesh, cp, cv, rho, k, speed);
boundary_condition_functions = {};

problem = @(w, t) energy_2D_fvm(w, t, mesh, cp, cv, rho, k, speed, ...
    convective_method, conductive_method, boundary_condition_functions);

% Stop the computation at 1.5 s. The stopping condition should be a
% function of the state vector and the current time, so configure it
% accordingly using an anonymous function that calls stop_at_time with 1.5
% s as stopping time.
condition = @(w, t) stop_at_time(t, 1.5);

% Use a sampling period of 0.01 s. Configure the results manager using an
% anonymous function, so the solver only sees a function that depends on
% the state vector, the time and the previous results structure.
manager = @(w, t, old_results) sample_results(w, t, old_results, 0.01);

% Use an implicit time integrators
impl_propagator = @bw_euler;

dt_fun = @(w, t) constant_dt(w, t, 0.01);

%% Solve with an implicit method.
impl_results = solver(w0, t0, problem, impl_propagator, ...
    dt_fun, condition, manager);
