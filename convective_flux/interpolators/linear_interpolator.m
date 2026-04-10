function [A, b] = linear_interpolator(state, cells)
    %LINEAR_INTERPOLATOR Assembles the convective flux matrix using linear (central) interpolation.
    %   [A, b] = LINEAR_INTERPOLATOR(state, cells) builds the spatial operator
    %   matrices for convective transport using arithmetic averaging of cell
    %   velocities to compute face velocities (second-order central scheme).
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
    %     Central-difference convective operator matrix.
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

        % diagonal
        A(cell_index, cell_index) = (right_face_velocity - left_face_velocity) / 2;
        A(num_cells + cell_index, num_cells + cell_index) = (right_face_velocity - left_face_velocity) / 2;
        A(2 * num_cells + cell_index, 2 * num_cells + cell_index) = (right_face_velocity - left_face_velocity) / 2;

        % subdiagonal inferior
        A(cell_index, neighbours.left) = -left_face_velocity / 2;
        A(num_cells + cell_index, num_cells + neighbours.left) = -left_face_velocity / 2;
        A(2 * num_cells + cell_index, 2 * num_cells + neighbours.left) = -left_face_velocity / 2;

        % subdiagonal superior
        A(cell_index, neighbours.right) = right_face_velocity / 2;
        A(num_cells + cell_index, num_cells + neighbours.right) = right_face_velocity / 2;
        A(2 * num_cells + cell_index, 2 * num_cells + neighbours.right) = right_face_velocity / 2;
    end

end
