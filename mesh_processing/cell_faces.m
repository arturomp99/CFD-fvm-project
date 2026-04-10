function faces = cell_faces(cell_nodes)
    %CELL_FACES Compute the face connectivity for a polygonal cell.
    %   faces = CELL_FACES(cell_nodes) creates a matrix of face connectivity
    %   where each row contains the indices of two consecutive nodes.
    %
    %   Input:
    %       cell_nodes - N-by-2 array of node coordinates [x, y].
    %
    %   Output:
    %       faces - N-by-2 array where each row [i, j] represents a face
    %               connecting nodes i and j (with wrap-around).

    n = size(cell_nodes, 1);
    faces = zeros(n, 2, 2);

    for i = 1:n
        next_i = get_next_vertex_index(i, n);
        faces(i, 1, :) = cell_nodes(i, :);
        faces(i, 2, :) = cell_nodes(next_i, :);
    end

end
