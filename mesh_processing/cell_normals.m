function [normals] = cell_normals(nodes, centroid)
    %CELL_NORMALS Compute outward unit normals for each face/edge of a polygonal cell.
    %   normals = CELL_NORMALS(nodes, centroid) calculates the outward-pointing
    %   unit normal vectors for each edge of the polygon.
    %
    %   Inputs:
    %       nodes    - N-by-2 array of node coordinates [x, y].
    %       centroid - 1-by-2 vector [cx, cy] of the cell centroid.
    %
    %   Output:
    %       normals - N-by-2 array where each row is the unit normal vector
    %                 for the corresponding edge, pointing outward.
    %
    %   For each edge from node i to node i+1, the normal is computed as the
    %   perpendicular to the edge vector, normalized, and oriented outward
    %   by checking the direction relative to the centroid.

    n = size(nodes, 1);
    normals = zeros(n, 2);

    for i = 1:n
        edge_vector = get_vector_to_next_vertex(nodes, i, n);

        normal = get_2D_vector_normal(edge_vector);

        % Midpoint of the edge
        midpoint = nodes(i, :) + edge_vector / 2;

        normal = correct_normal_direction(normal, midpoint, centroid);

        normals(i, :) = normal;
    end

end
