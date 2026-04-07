function [results] = solver( ...
        w0, ...
        t0, ...
        problem, ...
        propagator, ...
        dt_fun, ...
        stop_condition, ...
        manager)
    %   SOLVER: Solves a generic problem until a stopping condition is met.
    %
    %   Inputs
    %   ------
    %   w0 :
    %     Initial state vector. It is a column vector.
    %   t0 :
    %     Initial time. [s]
    %   problem :
    %     Function that computes the jacobian matrix and the independent terms
    %     vector of the problem. It accepts the state vector and the current
    %     time.
    %   propagator :
    %     Time integrator function. It accepts the state vector, the current
    %     time, problem function and the time step, and computes the state
    %     vector after integrating one time step.
    %   dt_fun :
    %     Timestep calculator. It accepts the state vector and the current
    %     time, and returns the time step in s.
    %   stop_condition :
    %     Stopping condition calculator. It accepts the state vector and the
    %     current time. Returns true if the simulation should stop.
    %   manager :
    %     Results manager function. Accepts the state vector, the time and the
    %     previous results structure, and returns an updated results structure.
    %     Can be used for online plotting, for data sampling, for adaptive mesh
    %     refinement...
    %
    %   Outputs
    %   -------
    %   results :
    %     Results structure, as returned by the results manager.
    results = [];
    w = w0;
    t = t0;

    while (stop_condition(w, t) == false)
        dt = dt_fun(w, t);
        w = propagator(w, t, dt, problem);
        t = t + dt;
        results = manager(w, t, results);
    end

end
