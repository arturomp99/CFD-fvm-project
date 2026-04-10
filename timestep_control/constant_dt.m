function [dt] = constant_dt(w, t, dt0)
    %CONSTANT_DT Returns a fixed time step regardless of the current state.
    %   dt = CONSTANT_DT(w, t, dt0) is a trivial timestep calculator that
    %   always returns the pre-configured constant dt0. The state vector w and
    %   time t are accepted for interface compatibility with the solver but
    %   are not used.
    %
    %   Inputs:
    %   -------
    %   w   : column vector - Current state vector (unused).
    %   t   : double        - Current time (unused). [s]
    %   dt0 : double        - Constant time step to return. [s]
    %
    %   Outputs:
    %   --------
    %   dt : double - Time step. [s]

    dt = dt0;

end
