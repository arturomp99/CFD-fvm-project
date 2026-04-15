function [A, b] = upwind_simple_interpolator(state, cells, boundary_info)
    %UPWIND_SIMPLE_INTERPOLATOR Assembles the convective flux matrix using a
    %   cell-local upwind scheme based on the sign of the bulk flow velocity.
    %
    %   Valid for SUPERSONIC flow only. When the flow velocity u exceeds the
    %   local speed of sound c, all three characteristic speeds (u-c, u, u+c)
    %   share the same sign, so a single upwind direction applies to all
    %   conserved variables. The scheme breaks down for subsonic flow where
    %   u-c < 0 < u+c and information propagates in both directions.
    %
    %   At each cell the physical flux Jacobian dF/dw is evaluated using the
    %   local primitive variables. The upwind flux is then linearised as:
    %     u > 0: dw_i/dt = -(J_i*w_i  - J_{i-1}*w_{i-1}) / dx
    %     u < 0: dw_i/dt = -(J_{i+1}*w_{i+1} - J_i*w_i)  / dx
    %
    %   State vector layout (stacked): [density(N); momentum(N); energy(N)].
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     Concatenated state vector [density; momentum; total energy].
    %   cells : struct array (1 x N)
    %     Mesh cells with 'connectivity', 'centroid', and 'volume' fields.
    %   boundary_info : struct (optional)
    %     Structure with .boundary_types cell array specifying 'open' or 'wall'
    %     for each boundary surface. (Currently not used in this interpolator)
    %
    %   Outputs:
    %   --------
    %   A : sparse matrix (3N x 3N)
    %     Upwind convective operator matrix (cell-local Jacobian based).
    %   b : column vector (3N x 1)
    %     Zero vector (no source terms added here).

    % Note: boundary_info is accepted for interface compatibility.

    gamma = Air.GAMMA;
    num_cells = length(cells);

    densities  = state(1:num_cells);
    momentums  = state(num_cells + 1 : 2*num_cells);
    energies   = state(2*num_cells + 1 : 3*num_cells);
    velocities = momentums ./ densities;
    pressures  = (gamma - 1) * (energies - 0.5 .* densities .* velocities .^ 2);
    enthalpies = (energies + pressures) ./ densities;

    A = sparse(3*num_cells, 3*num_cells);
    b = zeros(3*num_cells, 1);

    for i = 1:num_cells
        neighbours = get_neighour_cells(cells(i), cells);
        dx = cells(i).volume;  % cell length in 1D [m]

        cell_velocity  = velocities(i);
        cell_enthalpy  = enthalpies(i);
        cell_jacobian  = euler_flux_jacobian(cell_velocity, cell_enthalpy, gamma);

        % Row indices for cell i in the stacked state layout
        rows = [i, num_cells + i, 2*num_cells + i];

        if cell_velocity > 0
            % All waves travel right → upwind is the left cell.
            % dw_i/dt = -(J_i*w_i - J_{left}*w_{left}) / dx
            A(rows, rows) = A(rows, rows) - cell_jacobian / dx;

            if ~isempty(neighbours.left)
                left = neighbours.left;
                left_jacobian = euler_flux_jacobian(velocities(left), enthalpies(left), gamma);
                cols_left = [left, num_cells + left, 2*num_cells + left];
                A(rows, cols_left) = A(rows, cols_left) + left_jacobian / dx;
            end
        else
            % All waves travel left → upwind is the right cell.
            % dw_i/dt = -(J_{right}*w_{right} - J_i*w_i) / dx
            A(rows, rows) = A(rows, rows) + cell_jacobian / dx;

            if ~isempty(neighbours.right)
                right = neighbours.right;
                right_jacobian = euler_flux_jacobian(velocities(right), enthalpies(right), gamma);
                cols_right = [right, num_cells + right, 2*num_cells + right];
                A(rows, cols_right) = A(rows, cols_right) - right_jacobian / dx;
            end
        end

    end

end

% -------------------------------------------------------------------------

function J = euler_flux_jacobian(u, H, gamma)
    %EULER_FLUX_JACOBIAN Returns the 1D Euler flux Jacobian dF/dw at (u, H).
    %   Rows/columns correspond to [density, momentum, total energy].
    J = [0,                               1,                 0;
         (gamma - 3)/2 * u^2,             (3 - gamma)*u,     gamma - 1;
         (gamma - 1)/2 * u^3 - u * H,     H - (gamma-1)*u^2, gamma * u];
end

