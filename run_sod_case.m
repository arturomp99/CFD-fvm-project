function out = run_sod_case(nnodes, dt_value)
    %RUN_SOD_CASE Run the 1D Sod problem for a given mesh and fixed dt.

    % Paths
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

    [nodes_file, cells_file, bc_files] = mesh_case_paths(nnodes);

    % Mesh
    [cells, boundary_info] = mesh_processor(nodes_file, cells_file, bc_files);
    
    num_cells = length(cells);

    centroids = reshape([cells.centroid], 2, [])';
    centroids_x = centroids(:,1);

    % Initial condition
    w0 = Config.INITIAL_CONDITIONS(centroids_x);

    % Problem
    problem = @(state, time) fvm_1D_euler(state, cells, boundary_info);

    propagator = Config.PROPAGATOR;
    timestep_calculator = @(w, t) constant_dt(w, t, dt_value);
    stopping_condition = Config.STOPPING_CONDITION;
    manager = Config.RESULTS_MANAGER;

    % Solve + timing
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

    % Extract primitive variables at final sampled time
    t_vec = results(:,1);
    [~, idx_t] = min(abs(t_vec - Config.T_END));
    t_plot = t_vec(idx_t);

    rho  = results(idx_t, 2 : num_cells + 1).';
    rhou = results(idx_t, num_cells + 2 : 2 * num_cells + 1).';
    E    = results(idx_t, 2 * num_cells + 2 : 3 * num_cells + 1).';

    u = rhou ./ rho;
    p = (Config.GAMMA - 1) .* (E - 0.5 .* (rhou.^2) ./ rho);

    [x_sorted, idx] = sort(centroids_x);

    out.nnodes = nnodes;
    out.ncells = num_cells;
    out.x = x_sorted;
    out.rho = rho(idx);
    out.u = u(idx);
    out.p = p(idx);
    out.t = t_plot;
    out.runtime = runtime;
    out.dt = dt_value;
    out.nsteps_est = round((Config.T_END - Config.T0) / dt_value);
    out.results = results;
end