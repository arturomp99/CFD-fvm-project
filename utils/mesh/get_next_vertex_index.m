function next_vertex_index = get_next_vertex_index( ...
        current_vertex_index, ...
        vertices_total_number ...
    )
    %GET_NEXT_VERTEX_INDEX Returns the index of the next vertex in a polygon, wrapping around.
    %   next = GET_NEXT_VERTEX_INDEX(current, total) returns current+1, or 1 when
    %   current equals total (cyclic traversal of polygon vertices).
    %
    %   Inputs:
    %   -------
    %   current_vertex_index  : integer - 1-based index of the current vertex.
    %   vertices_total_number : integer - Total number of vertices in the polygon.
    %
    %   Outputs:
    %   --------
    %   next_vertex_index : integer - 1-based index of the next vertex.

    next_vertex_index = mod(current_vertex_index, vertices_total_number) + 1;
end
