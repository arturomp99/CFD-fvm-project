function [A] = cell_areas(cell_nodes)
    %CELL_AREAS Compute the lengths of each edge (face area in 2D) for a polygonal cell.
    %   A = CELL_AREAS(cell_nodes) calculates the length of each edge in the polygon.
    %
    %   Input:
    %       cell_nodes - N-by-2 array of node coordinates [x, y].
    %
    %   Output:
    %       A - 1-by-N vector of edge lengths, where A(i) is the length of the edge
    %           from node i to node i+1 (with wrap-around).

    n = size(cell_nodes, 1);
    A = zeros(1, n);

    for i = 1:n
        edge_vector = get_vector_to_next_vertex(cell_nodes, i, n);
        A(i) = norm(edge_vector);
    end

end
