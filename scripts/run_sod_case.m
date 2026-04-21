function out = run_sod_case(nnodes, dt_value, propagator_in)
%RUN_SOD_CASE Run the 1D Sod problem for a given mesh and fixed dt.
%
%   out = run_sod_case(nnodes, dt_value)
%   out = run_sod_case(nnodes, dt_value, propagator_in)
%
%   Inputs:
%   -------
%   nnodes : int
%       Mesh tag / number of nodes file identifier.
%   dt_value : double
%       Fixed time step used during the whole simulation.
%   propagator_in : function handle, optional
%       Temporal propagator to use. If omitted, Config.PROPAGATOR is used.
%
%   Output:
%   ------
%   out : struct
%       Structure with mesh info, primitive variables, runtime and results.

    % ---------------------------------------------------------
    % PATHS
    % ---------------------------------------------------------
    addpath('mesh_processing');
    addpath(genpath('utils'));
    addpath('constants');
    addpath('initial_conditions');
    addpath('problems/fvm_1D_euler');
    addpath(genpath('convective_flux'));
    addpath('stopping_criteria');
    addpath('timestep_control');
    addpath('propagators');
    addpath('results_manager');

    % ---------------------------------------------------------
    % OPTIONAL PROPAGATOR
    % ---------------------------------------------------------
    if nargin < 3 || isempty(propagator_in)
        propagator = Config.PROPAGATOR;
    else
        propagator = propagator_in;
    end

    % ---------------------------------------------------------
    % MESH PATHS
    % ---------------------------------------------------------
    [nodes_file, cells_file, bc_files] = mesh_case_paths(nnodes);

    % ---------------------------------------------------------
    % MESH
    % ---------------------------------------------------------
    cells = mesh_processor(nodes_file, cells_file, bc_files);
    num_cells = length(cells);

    centroids_x = get_cell_centroid_x(cells);
    centroids_x = centroids_x(:);

    % ---------------------------------------------------------
    % INITIAL CONDITION
    % ---------------------------------------------------------
    w0 = Config.INITIAL_CONDITIONS(centroids_x);

    % ---------------------------------------------------------
    % PROBLEM DEFINITION
    % ---------------------------------------------------------
    problem = @(state, time) fvm_1D_euler(state, cells, time);

    % Fixed time-step calculator
    timestep_calculator = @(w, t) constant_dt(w, t, dt_value);

    % Stopping condition and results manager
    stopping_condition = Config.STOPPING_CONDITION;
    manager = Config.RESULTS_MANAGER;

    % ---------------------------------------------------------
    % SOLVE + TIMING
    % ---------------------------------------------------------
    tic;
    results = solver( ...
        w0, ...
        Config.T0, ...
        problem, ...
        propagator, ...
        timestep_calculator, ...
        stopping_condition, ...
        manager ...
    );
    runtime = toc;

    % ---------------------------------------------------------
    % EXTRACT FINAL STATE
    % ---------------------------------------------------------
    if isempty(results)
        error('run_sod_case:EmptyResults', ...
            'The solver returned an empty results array.');
    end

    t_vec = results(:, 1);
    [~, idx_t] = min(abs(t_vec - Config.T_END));
    t_plot = t_vec(idx_t);

    rho  = results(idx_t, 2:num_cells + 1).';
    rhou = results(idx_t, num_cells + 2:2 * num_cells + 1).';
    E    = results(idx_t, 2 * num_cells + 2:3 * num_cells + 1).';

    % Primitive variables
    u = rhou ./ rho;
    p = (Air.GAMMA - 1) .* (E - 0.5 .* (rhou.^2) ./ rho);

    % Sort by x just in case
    [x_sorted, idx_sort] = sort(centroids_x);

    % ---------------------------------------------------------
    % OUTPUT
    % ---------------------------------------------------------
    out = struct();
    out.nnodes = nnodes;
    out.ncells = num_cells;
    out.x = x_sorted;
    out.rho = rho(idx_sort);
    out.u = u(idx_sort);
    out.p = p(idx_sort);
    out.t = t_plot;
    out.runtime = runtime;
    out.dt = dt_value;
    out.nsteps_est = round((Config.T_END - Config.T0) / dt_value);
    out.results = results;
    out.propagator = propagator;
end