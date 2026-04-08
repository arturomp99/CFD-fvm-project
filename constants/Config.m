classdef Config

    properties (Constant)
        INITIAL_CONDITIONS = @(pos) uniform_velocity(50); % es una función que depende de la posición
        SOURCE_TERMS = @(state, pos, t) 0;
        BOUNDARY_CONDITIONS = @(state, pos, t) 0;
    end

end
