function [A, b] = upwind_interpolator(state, cells)
    %UPWIND_INTERPOLATOR Assembles the convective flux matrix using first-order upwind scheme.
    %   [A, b] = UPWIND_INTERPOLATOR(state, cells) builds the spatial operator
    %   matrices for convective transport by applying the upwind (donor-cell)
    %   scheme: the flux at each face is taken from the upwind cell determined
    %   by the local flow velocity.
    %
    %   The state vector is ordered as [density (N); momentum (N); energy (N)].
    %   The resulting A matrix has a block-diagonal structure with three identical
    %   N x N blocks, one for each conserved variable.
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     Concatenated state vector [density; momentum; total energy].
    %   cells : struct array (1 x N)
    %     Mesh cells with 'connectivity' and 'centroid' fields.
    %
    %   Outputs:
    %   --------
    %   A : matrix (3N x 3N)
    %     Upwind convective operator matrix.
    %   b : column vector (3N x 1)
    %     Zero vector (no source terms added here).

    num_cells = length(cells);

    [A, b] = initialize_A_b(num_cells);

    for cell_index = 1:num_cells
        cell = cells(cell_index);
        neighbours = get_neighour_cells(cell, cells);

        cell_velocity = get_cell_velocity(cell_index, num_cells, state);
        left_neighbour_velocity = get_cell_velocity(neighbours.left, num_cells, state);
        right_neighbour_velocity = get_cell_velocity(neighbours.right, num_cells, state);

        [left_face_velocity, right_face_velocity] = ...
            get_velocity_at_faces(cell_velocity, left_neighbour_velocity, right_neighbour_velocity);

        if (cell_velocity > 0)
            A(cell_index, cell_index) = right_face_velocity;
            A(num_cells + cell_index, num_cells + cell_index) = right_face_velocity;
            A(2 * num_cells + cell_index, 2 * num_cells + cell_index) = right_face_velocity;
            A(cell_index, neighbours.left) = -left_face_velocity;
            A(num_cells + cell_index, num_cells + neighbours.left) = -left_face_velocity;
            A(2 * num_cells + cell_index, 2 * num_cells + neighbours.left) = -left_face_velocity;
        else
            A(cell_index, neighbours.right) = right_face_velocity;
            A(num_cells + cell_index, num_cells + neighbours.right) = right_face_velocity;
            A(2 * num_cells + cell_index, 2 * num_cells + neighbours.right) = right_face_velocity;
            A(cell_index, cell_index) = -left_face_velocity;
            A(num_cells + cell_index, num_cells + cell_index) = -left_face_velocity;
            A(2 * num_cells + cell_index, 2 * num_cells + cell_index) = -left_face_velocity;
        end

    end

end
