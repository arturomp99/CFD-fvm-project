function [new_results] = discard_results(w, t, old_results)
    %DISCARD_RESULTS Results manager that discards all data (no-op manager).
    %   new_results = DISCARD_RESULTS(w, t, old_results) always returns an empty
    %   array. Use this when you want the solver to run without storing any results.
    %
    %   Inputs:
    %   -------
    %   w           : column vector - Current state vector (unused).
    %   t           : double        - Current time (unused). [s]
    %   old_results : array         - Previous results (unused).
    %
    %   Outputs:
    %   --------
    %   new_results : [] - Always empty.

    new_results = [];

end
