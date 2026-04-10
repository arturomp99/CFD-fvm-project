function c = cell_centroid(cell_nodes)
    %CELL_CENTROID Compute the centroid of a polygonal cell given its node coordinates.
    %   centroid = CELL_CENTROID(cell_nodes) calculates the geometric centroid
    %   of a polygon defined by its vertices using MATLAB's polyshape.
    %
    %   Input:
    %       cell_nodes - N-by-2 array of node coordinates [x, y].
    %
    %   Output:
    %       centroid - 1-by-2 vector [cx, cy] representing the centroid coordinates.
    %
    %   This function uses MATLAB's polyshape to compute the centroid, which
    %   handles various polygon shapes and edge cases robustly.

    % TODO: Document that I've used this as a reference
    % https://es.mathworks.com/help/matlab/ref/polyshape.centroid.html
    n = size(cell_nodes, 1);

    if n < 3
        error('A polygon must have at least 3 vertices.');
    end

    % Create a polyshape from the node coordinates
    polygon = polyshape(cell_nodes(:, 1), cell_nodes(:, 2));

    % Get the centroid
    [c_x, c_y] = centroid(polygon);
    c = [c_x, c_y];
