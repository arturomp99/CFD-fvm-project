function [A, b] = hllc_interpolator(state, cells)
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

    for i = 1:num_cells
        Ui = [rho_sorted(i); rhou_sorted(i); E_sorted(i)];

        % Left face
        if i == 1
            % First cell
            UL = Config.LEFT_BOUNDARY_CONDITION(Ui);
            xL = x_sorted(i) - 0.5 * (x_sorted(i + 1) - x_sorted(i));
        else
            UL = [rho_sorted(i - 1); rhou_sorted(i - 1); E_sorted(i - 1)];
            xL = 0.5 * (x_sorted(i - 1) + x_sorted(i));
        end

        % Right face
        if i == num_cells
            % Last cell
            UR = Config.RIGHT_BOUNDARY_CONDITION(Ui);
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
