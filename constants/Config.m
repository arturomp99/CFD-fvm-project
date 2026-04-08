classdef Config

    properties (Constant)
        INITIAL_CONDITIONS = @(pos) ... % funcion de la posición
            uniform_velocity(0);
        SOURCE_TERMS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            ();
        BOUNDARY_CONDITIONS = @(state, pos, t) ... % funcion del vector de estado, la posición y del tiempo
            0;
    end

end
