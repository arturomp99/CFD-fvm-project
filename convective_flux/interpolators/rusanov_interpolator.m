function [A, b] = rusanov_interpolator(state, cells)
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
    %
    %   Output format:
    %       d(state)/dt = A*state + b
    %   We return A = 0 and b = RHS(state).

    num_cells = length(cells);
    gamma = Air.GAMMA;

    A = sparse(3 * num_cells, 3 * num_cells);
    b = zeros(3 * num_cells, 1);

    [ ...
         x_sorted, ...
         perm, ...
         rho_sorted, ...
         rhou_sorted, ...
         E_sorted ...
     ] = sort_state_by_cell_centroid_x(state, cells);

    rhs_rho = zeros(num_cells, 1);
    rhs_rhou = zeros(num_cells, 1);
    rhs_E = zeros(num_cells, 1);

    for i = 1:num_cells
        Ui = [rho_sorted(i); rhou_sorted(i); E_sorted(i)];

        % First cell
        if i == 1
            % First cell
            UL = Config.LEFT_BOUNDARY_CONDITION(Ui);
            xL = x_sorted(i) - 0.5 * (x_sorted(i + 1) - x_sorted(i));
        else
            UL = [rho_sorted(i - 1); rhou_sorted(i - 1); E_sorted(i - 1)];
            xL = 0.5 * (x_sorted(i - 1) + x_sorted(i));
        end

        % Last cell
        if i == num_cells
            % Last cell
            UR = Config.RIGHT_BOUNDARY_CONDITION(Ui);
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

    rhs_rho_orig = zeros(num_cells, 1);
    rhs_rhou_orig = zeros(num_cells, 1);
    rhs_E_orig = zeros(num_cells, 1);

    rhs_rho_orig(perm) = rhs_rho;
    rhs_rhou_orig(perm) = rhs_rhou;
    rhs_E_orig(perm) = rhs_E;

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
