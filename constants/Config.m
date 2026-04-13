classdef Config

    properties (Constant)
        % Problem definition
        INITIAL_CONDITIONS = @(pos) ... % funcion de la posición
            uniform( ...
            Air.SEA_LEVEL_PRESSURE, ...
            Air.SEA_LEVEL_DENSITY, ...
            0., ...
            pos ...
        );
        %     sod( ...
        %     struct('left', 1., 'right', 0.1), ... % pressure
        %     struct('left', 10., 'right', 0.125), ... % density
        %     struct('left', 0.1, 'right', 0.1), ... % velocity
        %     0.5, ...
        %     pos ...
        % );
        SOURCE_TERMS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            point_source();
        BOUNDARY_CONDITIONS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            0;

        % Solver
        IS_SOURCE_IMPLICIT = false;
        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) ... % los interpoladores estan definidos en convective_flux\interpolators
            upwind_interpolator(state, cells);
        PROPAGATOR = @(state, time, d_time, problem) ...
            bw_euler(state, time, d_time, problem);
        TIMESTEP_CALCULATOR = @(w, t) ...
            constant_dt(w, t, 0.01);
        STOPPING_CONDITION = @(w, t) ...
            stop_at_time(t, 1.5);
        RESULTS_MANAGER = @(w, t, old_results) ...
            sample_results(w, t, old_results, 0.01);

        % Solution
        DIFFUMINATED_VISUALIZATION = false;
    end

end
