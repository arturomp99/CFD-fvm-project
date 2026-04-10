function normal = get_2D_vector_normal(edge_vector)
    %GET_2D_VECTOR_NORMAL Computes the unit normal of a 2D edge vector.
    %   normal = GET_2D_VECTOR_NORMAL(edge_vector) rotates the edge vector 90 degrees
    %   counter-clockwise and normalises the result to unit length.
    %
    %   Inputs:
    %   -------
    %   edge_vector : 1x2 vector - Vector along the edge [dx, dy].
    %
    %   Outputs:
    %   --------
    %   normal : 1x2 vector - Unit normal perpendicular to the edge.
    %     Returns the zero vector if the edge has zero length.

    % Perpendicular vector (rotate 90 degrees counterclockwise)
    normal = [-edge_vector(2), edge_vector(1)];

    % Normalize
    normal_length = norm(normal);

    if normal_length > 0
        normal = normal / normal_length;
    end

end
