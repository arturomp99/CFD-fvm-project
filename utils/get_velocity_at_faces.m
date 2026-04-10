function [left_face_velocity, right_face_velocity] = get_velocity_at_faces( ...
        cell_velocity, ...
        left_neighbour_velocity, ...
        right_neighbour_velocity ...
    )
    %GET_VELOCITY_AT_FACES Interpolates velocities at the left and right faces of a cell.
    %   [vL, vR] = GET_VELOCITY_AT_FACES(v_c, v_left, v_right) computes the
    %   face-centred velocities using arithmetic averaging between the cell and
    %   each of its neighbours.
    %
    %   Inputs:
    %   -------
    %   cell_velocity           : double  - Velocity at the current cell centre. [m/s]
    %   left_neighbour_velocity : double  - Velocity at the left neighbour cell centre. [m/s]
    %   right_neighbour_velocity: double  - Velocity at the right neighbour cell centre. [m/s]
    %
    %   Outputs:
    %   --------
    %   left_face_velocity  : double - Interpolated velocity at the left face. [m/s]
    %   right_face_velocity : double - Interpolated velocity at the right face. [m/s]

    left_face_velocity = (left_neighbour_velocity + cell_velocity) / 2;
    right_face_velocity = (right_neighbour_velocity + cell_velocity) / 2;
end
