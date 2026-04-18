classdef Config

    properties (Constant)
        T0 = 0.0;
        T_END = 0.2;
        SAMPLE_DT = 1e-3;

        % Boundary condition types
        % 'open'     - transmissive/open boundary (waves pass through)
        % 'wall'     - solid wall boundary (reflective, zero normal velocity)
        % 'velocity' - constant velocity boundary (requires BOUNDARY_VELOCITIES)
        %
        % Configure per boundary surface. Order matches FilePaths.BOUNDARY_CONDITIONS.
        BOUNDARY_TYPES = { ...
                              'open', ...     % 1: bc_bottom
                              'velocity', ... % 2: bc_left
                              'open', ...     % 3: bc_right
                              'open', ...     % 4: bc_top
                              'open' ...      % 5: bc_whole_contour
                          };

        % Velocity values for 'velocity' boundary conditions [m/s]
        % Only used when the corresponding BOUNDARY_TYPES entry is 'velocity'.
        % Set to NaN for surfaces that don't use velocity BC.
        BOUNDARY_VELOCITIES = [NaN, 100.0, NaN, NaN, NaN];

        % Problem definition
        INITIAL_CONDITIONS = @(pos) ... % funcion de la posición
            uniform( ...
            Air.SEA_LEVEL_PRESSURE, ...
            Air.SEA_LEVEL_DENSITY, ...
            0.0, ...
            pos ...
        );
        %     sod( ...
        %     struct('left', 1.0, 'right', 0.1), ...
        %     struct('left', 1.0, 'right', 0.125), ...
        %     struct('left', 0.0, 'right', 0.0), ...
        %     0.5, ...
        %     pos ...
        % );

        SOURCE_TERMS = @(state, pos, t) ... % funcion de la posición
            point_source();

        BOUNDARY_CONDITIONS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            0;

        % Solver
        IS_SOURCE_IMPLICIT = false;

        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) ... % los interpoladores estan definidos en convective_flux\interpolators
            hllc_interpolator(state, cells);

        PROPAGATOR = @(state, time, d_time, problem) ...
            bw_euler(state, time, d_time, problem);

        TIMESTEP_CALCULATOR = @(w, t) ...
            constant_dt(w, t, 5e-5);

        STOPPING_CONDITION = @(w, t) ...
            stop_at_time(t, Config.T_END);

        RESULTS_MANAGER = @(w, t, old_results) ...
            sample_results(w, t, old_results, Config.SAMPLE_DT);

        % Solution
        DIFFUMINATED_VISUALIZATION = false;
    end

end
