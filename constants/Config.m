classdef Config

    properties (Constant)
        % Problem definition
        INITIAL_CONDITIONS = @(pos) ... % funcion de la posición
            sod( ...
            struct('left', 1., 'right', 0.1), ... % pressure
            struct('left', 1., 'right', 0.125), ... % density
            struct('left', 0., 'right', 0.0), ... % velocity
            0.5, ...
            pos ...
        );
        SOURCE_TERMS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            point_source();
        BOUNDARY_CONDITIONS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            0;

        % Solver
        IS_SOURCE_IMPLICIT = false;
        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) ... % los interpoladores estan definidos en convective_flux\interpolators
            upwind_interpolator(state, cells);
    end

end
