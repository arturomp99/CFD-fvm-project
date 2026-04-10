function corrected_normal = correct_normal_direction( ...
        normal, ...
        edge_midpoint, ...
        centroid ...
    )
    %CORRECT_NORMAL_DIRECTION Flips a face normal so that it points outward from the cell.
    %   corrected_normal = CORRECT_NORMAL_DIRECTION(normal, edge_midpoint, centroid)
    %   checks whether the provided normal points away from the cell centroid by
    %   computing the dot product with the centroid-to-midpoint vector. If it points
    %   inward, the normal is negated.
    %
    %   Inputs:
    %   -------
    %   normal       : 1x2 vector - Unit normal to be corrected.
    %   edge_midpoint: 1x2 vector - Midpoint of the face/edge. [m]
    %   centroid     : 1x2 vector - Centroid of the owning cell. [m]
    %
    %   Outputs:
    %   --------
    %   corrected_normal : 1x2 vector - Outward-pointing unit normal.

    % Vector from midpoint to centroid
    to_centroid = centroid - edge_midpoint;

    % If dot product is negative, flip the normal to point outward
    if dot(normal, to_centroid) < 0
        corrected_normal = -normal;
    else
        corrected_normal = normal;
    end

end
