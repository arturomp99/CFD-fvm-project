function face_avg_x = get_face_average_x(cell_nodes, face_idx)
    %GET_FACE_AVERAGE_X Compute the average x-coordinate of a face.
    %   face_avg_x = GET_FACE_AVERAGE_X(cell_nodes, face_idx) returns the
    %   average x-coordinate of the two vertices that form the face.
    %
    %   Inputs:
    %   -------
    %   cell_nodes : matrix (Vx2)
    %       Array of [x, y] coordinates for each vertex.
    %   face_idx : integer
    %       Index of the face (1 to num_vertices).
    %
    %   Outputs:
    %   --------
    %   face_avg_x : double
    %       Average x-coordinate of the face.

    num_vertices = size(cell_nodes, 1);
    
    v1_idx = face_idx;
    v2_idx = get_next_vertex_index(face_idx, num_vertices);
    
    x1 = cell_nodes(v1_idx, 1);
    x2 = cell_nodes(v2_idx, 1);
    
    face_avg_x = (x1 + x2) / 2;
end
