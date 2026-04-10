function [dt] = CFL_dt(CFL, c, dx)
    %CFL_DT Computes the time step from the CFL stability condition.
    %   dt = CFL_DT(CFL, c, dx) calculates the maximum stable time step
    %   dt = CFL * min(dx / c) across all cells.
    %
    %   Inputs:
    %   -------
    %   CFL : double
    %     Courant-Friedrichs-Lewy number (typically <= 1 for explicit schemes).
    %   c : double
    %     Wave propagation speed. [m/s]
    %   dx : vector
    %     Cell sizes. [m]
    %
    %   Outputs:
    %   --------
    %   dt : double
    %     Time step satisfying the CFL condition. [s]

    dt = CFL * min(dx ./ c);

end
