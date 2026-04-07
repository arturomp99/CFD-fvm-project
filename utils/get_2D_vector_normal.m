function normal = get_2D_vector_normal(vector)
    % Perpendicular vector (rotate 90 degrees counterclockwise)
    normal = [-edge_vector(2), edge_vector(1)];

    % Normalize
    normal_length = norm(normal);

    if normal_length > 0
        normal = normal / normal_length;
    end

end
