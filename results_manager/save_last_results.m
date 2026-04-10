function [new_results] = save_last_results(w, t, old_results)
    %SAVE_LAST_RESULTS Results manager that keeps only the most recent time step.
    %   new_results = SAVE_LAST_RESULTS(w, t, old_results) replaces the stored
    %   results with a column vector [t; w] at every call, so only the final
    %   state is retained at the end of the simulation.
    %
    %   Inputs:
    %   -------
    %   w           : column vector (N x 1) - Current state vector.
    %   t           : double                - Current time. [s]
    %   old_results : array                 - Previous results (discarded).
    %
    %   Outputs:
    %   --------
    %   new_results : column vector ((N+1) x 1) - [t; w] for the current step.

    new_results = [t; w];

end
