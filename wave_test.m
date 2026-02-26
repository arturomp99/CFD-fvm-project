%% 1D wave equation tester

%% Paths

clear all;
close all;

addpath('problems/wave');
addpath('results_manager');
addpath('propagators');
addpath('stopping_criteria');
addpath('timestep_control');

%% Configuration

% Domain length and number of cells
L = 1;
n = 100;

% Compute geometry and mesh
% Cell interface positions
x_interfaces = linspace(0, L, n + 1)';
% Cell centroid coordinates
x_centroids = (x_interfaces(1:end-1) + x_interfaces(2:end)) / 2;
% Cell lengths
dx_cells = (x_interfaces(2:end) - x_interfaces(1:end - 1));

% Areas
S = ones(n + 1, 1);
% Cell volumes
V = dx_cells .* (S(1:end-1) + S(2:end)) / 2;

% Propagation speed
c = 1;

% Initial conditions
t0 = 0;
w0 = zeros(n, 1);
w0((x_centroids < (L / 3)) & (x_centroids > (L / 5))) = 2;

% Boundary condition fluxes function
bcs = @(V, S, c, t) open_ends(V, S, c);

% Internal fluxes method
internal_fluxes = @upwind;

% Problem function. Is has to be a function of the state vector and the
% current time. We can configure it using an anonymous (or lambda) function
% which calls the wave equation function passing it the mesh geometry, the
% function that computes the internal fluxes and the function that computes
% the boundary condition fluxes.
problem = @(w, t) wave(w, t, V, S, c, internal_fluxes, bcs);

% Configure the time step using a CFL of 1.
CFL = 1;
dt_fun = @(w, t) CFL_dt(CFL, c, dx_cells);

% Use two time integrators
expl_propagator = @fw_euler;
impl_propagator = @bw_euler;

% Stop the computation at 1.5 s. The stopping condition should be a
% function of the state vector and the current time, so configure it
% accordingly using an anonymous function that calls stop_at_time with 1.5
% s as stopping time.
condition = @(w, t) stop_at_time(t, 1.5);

% Use a sampling period of 0.01 s. Configure the results manager using an
% anonymous function, so the solver only sees a function that depends on
% the state vector, the time and the previous results structure.
manager = @(w, t, old_results) sample_results(w, t, old_results, 0.01);

%% Solve with an explicit method.
expl_results = solver(w0, t0, problem, expl_propagator, ...
    dt_fun, condition, manager);

%% Solve with an implicit method.
impl_results = solver(w0, t0, problem, impl_propagator, ...
    dt_fun, condition, manager);

%% Call the postprocessing functions.
[X, T] = meshgrid(x_centroids, expl_results(:, 1));
figure(1);
subplot(1, 2, 1);
contourf(X, T, expl_results(:, 2:end));
colorbar();

[X, T] = meshgrid(x_centroids, expl_results(:, 1));
subplot(1, 2, 2);
contourf(X, T, expl_results(:, 2:end));
colorbar();