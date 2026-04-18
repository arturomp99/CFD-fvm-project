function [A, b] = hllc_interpolator(state, cells)
    %HLLC_INTERPOLATOR
    %   First-order finite-volume residual for 1D Euler equations
    %   using HLLC approximate Riemann solver.
    %
    %   Output format:
    %       d(state)/dt = A*state + b
    %   Here A = 0 and b = RHS(state).

    num_cells = length(cells);

    A = sparse(3 * num_cells, 3 * num_cells);
    b = zeros(3 * num_cells, 1);

    rho = state(1:num_cells);
    rhou = state(num_cells + 1:2 * num_cells);
    E = state(2 * num_cells + 1:3 * num_cells);

    x = zeros(num_cells, 1);

    for i = 1:num_cells
        x(i) = cells(i).centroid(1);
    end

    [x_sorted, perm] = sort(x);
    rho_sorted = rho(perm);
    rhou_sorted = rhou(perm);
    E_sorted = E(perm);

    rhs_rho = zeros(num_cells, 1);
    rhs_rhou = zeros(num_cells, 1);
    rhs_E = zeros(num_cells, 1);

    for i = 1:num_cells
        Ui = [rho_sorted(i); rhou_sorted(i); E_sorted(i)];

        if i == 1
            UL = Ui;
            xL = x_sorted(i) - 0.5 * (x_sorted(i + 1) - x_sorted(i));
        else
            UL = [rho_sorted(i - 1); rhou_sorted(i - 1); E_sorted(i - 1)];
            xL = 0.5 * (x_sorted(i - 1) + x_sorted(i));
        end

        if i == num_cells
            UR = Ui;
            xR = x_sorted(i) + 0.5 * (x_sorted(i) - x_sorted(i - 1));
        else
            UR = [rho_sorted(i + 1); rhou_sorted(i + 1); E_sorted(i + 1)];
            xR = 0.5 * (x_sorted(i) + x_sorted(i + 1));
        end

        F_left = hllc_flux(UL, Ui, Air.GAMMA);
        F_right = hllc_flux(Ui, UR, Air.GAMMA);

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
