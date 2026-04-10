function [A, b] = linear_interpolator(state, cells)
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
