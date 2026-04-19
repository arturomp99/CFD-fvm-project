function [A, b] = parallel_hllc_interpolator(state, cells)
    %HLLC_INTERPOLATOR
    %
    %   Implementa el solver de Riemann aproximado HLLC para las ecuaciones de
    %   Euler 1D.
    %
    %   ALGORITMO HLLC:
    %   ===============
    %   1. Calcula velocidades de onda S_L, S_R, S_* usando estimaciones de Roe/Davis
    %   2. Identifica región del problema de Riemann según velocidades de onda
    %   3. Retorna flujo apropiado según la región:
    %      - F_L si S_L > 0 (supersónico izquierdo)
    %      - F_*L si S_L ≤ 0 < S_* (subsónico izquierdo)
    %      - F_*R si S_* ≤ 0 < S_R (subsónico derecho)
    %      - F_R si S_R ≤ 0 (supersónico derecho)
    %
    %   Input
    %   ---------------------
    %   state : double (3*N×1)
    %       Vector de estado [ρ₁...ρₙ; (ρu)₁...(ρu)ₙ; E₁...Eₙ]
    %
    %   cells : struct array (1×N)
    %       Estructura de células con campos: .centroid, .volume, etc.
    %
    %   Output
    %   --------
    %   A : sparse (3*N×3*N)
    %       Matriz jacobiana (=0 para esquemas no lineales)
    %
    %   b : double (3*N×1)
    %       Vector RHS con divergencia de flujos HLLC

    num_cells = length(cells);

    [A, b] = initialize_A_b(num_cells);

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

    U_all = [rho_sorted'; rhou_sorted'; E_sorted'];

    % Left faces
    UL_all = zeros(3, num_cells);
    xL_all = zeros(num_cells, 1);
    UL_all(:, 2:end) = U_all(:, 1:end - 1);
    UL_all(:, 1) = Config.LEFT_BOUNDARY_CONDITION(U_all(:, 1)); % Boundary
    xL_all(2:end) = 0.5 * (x_sorted(1:end - 1) + x_sorted(2:end));
    xL_all(1) = x_sorted(1) - 0.5 * (x_sorted(2) - x_sorted(1)); % Boundary

    % Right faces
    UR_all = zeros(3, num_cells);
    xR_all = zeros(num_cells, 1);
    UR_all(:, 1:end - 1) = U_all(:, 2:end);
    UR_all(:, end) = Config.RIGHT_BOUNDARY_CONDITION(U_all(:, end)); % Boundary
    xR_all(1:end - 1) = 0.5 * (x_sorted(1:end - 1) + x_sorted(2:end));
    xR_all(end) = x_sorted(end) + 0.5 * (x_sorted(end) - x_sorted(end - 1)); % Boundary

    for i = 1:num_cells
        U_i = U_all(:, i);
        UL_i = UL_all(:, i);
        UR_i = UR_all(:, i);
        xR_i = xR_all(i);
        xL_i = xL_all(i);

        F_left = hllc_flux(UL_i, U_i, Air.GAMMA);
        F_right = hllc_flux(U_i, UR_i, Air.GAMMA);

        dx_i = xR_i - xL_i;
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
