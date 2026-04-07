function corrected_normal = correct_normal_direction( ...
        normal, ...
        edge_midpoint, ...
        centroid ...
    )
    % Vector from midpoint to centroid
    to_centroid = centroid - edge_midpoint;

    % If dot product is negative, flip the normal to point outward
    if dot(normal, to_centroid) < 0
        corrected_normal = -normal;
    else
        corrected_normal = normal;
    end

end
