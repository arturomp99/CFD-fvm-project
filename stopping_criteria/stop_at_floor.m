function [stop] = stop_at_floor(w, t)
    %STOP_AT_FLOOR Stops the simulation when the projectile reaches the ground (y <= 0).
    %   stop = STOP_AT_FLOOR(w, t) returns true when the y-position component of the
    %   state vector becomes negative, indicating the projectile has hit the floor.
    %   Intended for use with the cannon problem where w = [x; y; vx; vy].
    %
    %   Inputs:
    %   -------
    %   w : column vector - Current state vector; w(2) is the y-position. [m, m, m/s, m/s]
    %   t : double        - Current time (unused). [s]
    %
    %   Outputs:
    %   --------
    %   stop : logical - True when the projectile is at or below y = 0.

    stop = true;

    if (w(2) >= 0.)
        stop = false;
    end

end
