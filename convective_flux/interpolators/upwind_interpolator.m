function [A, b] = upwind_interpolator(state, cells)
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
            A(cell_index) = right_face_velocity;
            A(num_cells + cell_index) = right_face_velocity;
            A(2 * num_cells + cell_index) = right_face_velocity;
            A(neighbours.left) = -left_face_velocity;
            A(num_cells + neighbours.left) = -left_face_velocity;
            A(2 * num_cells + neighbours.left) = -left_face_velocity;
        else
            A(neighbours.right) = right_face_velocity;
            A(num_cells + neighbours.right) = right_face_velocity;
            A(2 * num_cells + neighbours.right) = right_face_velocity;
            A(cell_index) = -left_face_velocity;
            A(num_cells + cell_index) = -left_face_velocity;
            A(2 * num_cells + cell_index) = -left_face_velocity;
        end

    end

end
