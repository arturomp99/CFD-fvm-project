function vector_to_next_vertex = get_vector_to_next_vertex( ...
        polygon_nodes, ...
        current_vertex_index, ...
        vertices_total_number ...
    )
    %GET_VECTOR_TO_NEXT_VERTEX Computes the edge vector from a vertex to the next one.
    %   v = GET_VECTOR_TO_NEXT_VERTEX(polygon_nodes, current_vertex_index, total)
    %   returns the vector from the current vertex to the next vertex in the
    %   polygon (with cyclic wrap-around).
    %
    %   Inputs:
    %   -------
    %   polygon_nodes         : N x 2 array - Coordinates of all polygon vertices [x, y]. [m]
    %   current_vertex_index  : integer      - 1-based index of the current vertex.
    %   vertices_total_number : integer      - Total number of vertices N.
    %
    %   Outputs:
    %   --------
    %   vector_to_next_vertex : 1x2 vector - Edge vector [dx, dy] to the next vertex. [m]

    next_vertex_index = get_next_vertex_index( ...
        current_vertex_index, ...
        vertices_total_number ...
    );
    vector_to_next_vertex = ...
        polygon_nodes(next_vertex_index, :) - polygon_nodes(current_vertex_index, :);
end
