classdef Config

    properties (Constant)
        GAMMA = 1.4;

        T0 = 0.0;
        T_END = 0.2;
        SAMPLE_DT = 1e-3;

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

        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) ... % los interpoladores estan definidos en convective_flux\interpolators
            rusanov_interpolator(state, cells);

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
