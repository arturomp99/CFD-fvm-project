function [A, b] = rusanov_interpolator(state, cells, boundary_info)
    %RUSANOV_INTERPOLATOR
    %   First-order finite-volume residual for 1D Euler equations
    %   using Rusanov (local Lax-Friedrichs) numerical flux.
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     State vector [density; momentum; total energy].
    %   cells : struct array (1 x N)
    %     Mesh cells with geometry and boundary information.
    %   boundary_info : struct
    %     Structure with .boundary_types cell array specifying 'open' or 'wall'
    %     for each boundary surface.
    %
    %   Output format:
    %       d(state)/dt = A*state + b
    %   We return A = 0 and b = RHS(state).

    if nargin < 3
        boundary_info = struct('boundary_types', {{}});
    end

    num_cells = length(cells);
    gamma = Air.GAMMA;

    A = sparse(3 * num_cells, 3 * num_cells);
    b = zeros(3 * num_cells, 1);

    rho = state(1:num_cells);
    rhou = state(num_cells + 1:2 * num_cells);
    E = state(2 * num_cells + 1:3 * num_cells);

    x = zeros(num_cells, 1);

    for i = 1:num_cells
        x(i) = cells(i).centroid(1);
    end

    % Work in sorted x-order
    [x_sorted, perm] = sort(x);

    rho_sorted = rho(perm);
    rhou_sorted = rhou(perm);
    E_sorted = E(perm);

    % Map boundary info to sorted order
    cells_sorted = cells(perm);

    rhs_rho = zeros(num_cells, 1);
    rhs_rhou = zeros(num_cells, 1);
    rhs_E = zeros(num_cells, 1);

    for i = 1:num_cells
        Ui = [rho_sorted(i); rhou_sorted(i); E_sorted(i)];

        % Left state at left face
        if i == 1
            % Apply boundary condition
            [bc_type, bc_params] = get_boundary_type( ...
                cells_sorted(i), ...
                'left', ...
                boundary_info ...
            );
            UL = apply_bc_state(Ui, bc_type, -1, bc_params);
            xL = x_sorted(i) - 0.5 * (x_sorted(i + 1) - x_sorted(i));
        else
            UL = [rho_sorted(i - 1); rhou_sorted(i - 1); E_sorted(i - 1)];
            xL = 0.5 * (x_sorted(i - 1) + x_sorted(i));
        end

        % Right state at right face
        if i == num_cells
            % Apply boundary condition
            [bc_type, bc_params] = get_boundary_type(cells_sorted(i), 'right', boundary_info);
            UR = apply_bc_state(Ui, bc_type, +1, bc_params);
            xR = x_sorted(i) + 0.5 * (x_sorted(i) - x_sorted(i - 1));
        else
            UR = [rho_sorted(i + 1); rhou_sorted(i + 1); E_sorted(i + 1)];
            xR = 0.5 * (x_sorted(i) + x_sorted(i + 1));
        end

        F_left = rusanov_flux(UL, Ui, gamma);
        F_right = rusanov_flux(Ui, UR, gamma);

        dx_i = xR - xL;

        rhs_i =- (F_right - F_left) / dx_i;

        rhs_rho(i) = rhs_i(1);
        rhs_rhou(i) = rhs_i(2);
        rhs_E(i) = rhs_i(3);
    end

    % Map back to original cell ordering
    invperm = zeros(num_cells, 1);
    invperm(perm) = 1:num_cells;

    rhs_rho_orig = rhs_rho(invperm);
    rhs_rhou_orig = rhs_rhou(invperm);
    rhs_E_orig = rhs_E(invperm);

    b = [rhs_rho_orig; rhs_rhou_orig; rhs_E_orig];
end

function F = rusanov_flux(UL, UR, gamma)
    FL = euler_flux(UL, gamma);
    FR = euler_flux(UR, gamma);

    alpha = max(max_wave_speed(UL, gamma), max_wave_speed(UR, gamma));

    F = 0.5 * (FL + FR) - 0.5 * alpha * (UR - UL);
end

function F = euler_flux(U, gamma)
    rho = U(1);
    rhou = U(2);
    E = U(3);

    rho = max(rho, 1e-10);

    u = rhou / rho;
    p = (gamma - 1) * (E - 0.5 * rhou ^ 2 / rho);
    p = max(p, 1e-10);

    F = [
         rhou;
         rhou * u + p;
         u * (E + p)
         ];
end

function a = max_wave_speed(U, gamma)
    rho = U(1);
    rhou = U(2);
    E = U(3);

    rho = max(rho, 1e-10);

    u = rhou / rho;
    p = (gamma - 1) * (E - 0.5 * rhou ^ 2 / rho);
    p = max(p, 1e-10);

    c = sqrt(gamma * p / rho);
    a = abs(u) + c;
end
