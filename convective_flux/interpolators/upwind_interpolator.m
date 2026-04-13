function [A, b] = upwind_interpolator(state, cells)
    %UPWIND_INTERPOLATOR
    %   First-order finite-volume discretisation for 1D Euler
    %   using Rusanov (local Lax-Friedrichs) numerical flux.
    %
    %   Returns:
    %       d(state)/dt = A*state + b
    %   Here A = 0 and b contains the nonlinear spatial residual.

    num_cells = length(cells);
    gamma = Config.GAMMA;

    A = sparse(3 * num_cells, 3 * num_cells);
    b = zeros(3 * num_cells, 1);

    rho  = state(1:num_cells);
    rhou = state(num_cells + 1 : 2 * num_cells);
    E    = state(2 * num_cells + 1 : 3 * num_cells);

    x = zeros(num_cells, 1);
    for i = 1:num_cells
        x(i) = cells(i).centroid(1);
    end

    rhs_rho  = zeros(num_cells, 1);
    rhs_rhou = zeros(num_cells, 1);
    rhs_E    = zeros(num_cells, 1);

    for i = 1:num_cells
        neighbours = get_neighour_cells(cells(i), cells);

        Ui = [rho(i); rhou(i); E(i)];

        % Left face
        if isempty(neighbours.left)
            UL = Ui;  % transmissive BC
            xL = x(i) - 0.5 * local_dx(i, x, cells);
        else
            jL = neighbours.left;
            UL = [rho(jL); rhou(jL); E(jL)];
            xL = 0.5 * (x(jL) + x(i));
        end

        UR_left = Ui;
        F_left = rusanov_flux(UL, UR_left, gamma);

        % Right face
        if isempty(neighbours.right)
            UR = Ui;  % transmissive BC
            xR = x(i) + 0.5 * local_dx(i, x, cells);
        else
            jR = neighbours.right;
            UR = [rho(jR); rhou(jR); E(jR)];
            xR = 0.5 * (x(i) + x(jR));
        end

        UL_right = Ui;
        F_right = rusanov_flux(UL_right, UR, gamma);

        dx_i = xR - xL;

        rhs_i = -(F_right - F_left) / dx_i;

        rhs_rho(i)  = rhs_i(1);
        rhs_rhou(i) = rhs_i(2);
        rhs_E(i)    = rhs_i(3);
    end

    b = [rhs_rho; rhs_rhou; rhs_E];
end


function dx = local_dx(i, x, cells)
    neighbours = get_neighour_cells(cells(i), cells);

    if ~isempty(neighbours.left) && ~isempty(neighbours.right)
        dx = 0.5 * (x(neighbours.right) - x(neighbours.left));
    elseif ~isempty(neighbours.right)
        dx = x(neighbours.right) - x(i);
    elseif ~isempty(neighbours.left)
        dx = x(i) - x(neighbours.left);
    else
        dx = 1.0;
    end

    dx = max(dx, 1e-12);
end


function F = rusanov_flux(UL, UR, gamma)
    FL = euler_flux(UL, gamma);
    FR = euler_flux(UR, gamma);

    alpha = max(max_wave_speed(UL, gamma), max_wave_speed(UR, gamma));

    F = 0.5 * (FL + FR) - 0.5 * alpha * (UR - UL);
end


function F = euler_flux(U, gamma)
    rho  = max(U(1), 1e-12);
    rhou = U(2);
    E    = U(3);

    u = rhou / rho;
    p = (gamma - 1) * (E - 0.5 * rhou^2 / rho);
    p = max(p, 1e-12);

    F = [ ...
        rhou; ...
        rhou * u + p; ...
        u * (E + p) ...
    ];
end


function a = max_wave_speed(U, gamma)
    rho  = max(U(1), 1e-12);
    rhou = U(2);
    E    = U(3);

    u = rhou / rho;
    p = (gamma - 1) * (E - 0.5 * rhou^2 / rho);
    p = max(p, 1e-12);

    c = sqrt(gamma * p / rho);
    a = abs(u) + c;
end