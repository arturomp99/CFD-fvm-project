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
        % Configure the type for each boundary surface defined in FilePaths.
        % The order must match FilePaths.BOUNDARY_CONDITIONS.
        % Example: {'wall', 'open', 'open', 'wall', 'wall'}
        BOUNDARY_TYPES = { ...
            'open', ...   % bc_bottom
            'wall', ...   % bc_left (solid wall)
            'wall', ...   % bc_right (solid wall)
            'open', ...   % bc_top
            'open' ...    % bc_whole_contour
        };
        
        % Velocity values for 'velocity' boundary conditions [m/s]
        % Only used when the corresponding BOUNDARY_TYPES entry is 'velocity'.
        % Set to NaN for surfaces that don't use velocity BC.
        % Example: [NaN, 100.0, -50.0, NaN, NaN] for inlet at 100 m/s on surface 2
        BOUNDARY_VELOCITIES = [NaN, NaN, NaN, NaN, NaN];

        % Problem definition
        INITIAL_CONDITIONS = @(pos) ... % funcion de la posición
            sod( ...
            struct('left', 1.0, 'right', 0.1), ...
            struct('left', 1.0, 'right', 0.125), ...
            struct('left', 0.0, 'right', 0.0), ...
            0.5, ...
            pos ...
        );

        SOURCE_TERMS = @(state, pos, t) ... % funcion de la posición
            point_source();

        BOUNDARY_CONDITIONS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            0;

        % Solver
        IS_SOURCE_IMPLICIT = false;

        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells, boundary_info) ... % los interpoladores estan definidos en convective_flux\interpolators
            rusanov_interpolator(state, cells, boundary_info);

        PROPAGATOR = @(state, time, d_time, problem) ...
            fw_euler(state, time, d_time, problem);

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
