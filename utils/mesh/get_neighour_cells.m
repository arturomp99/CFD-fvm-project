function neighbours = get_neighour_cells( ...
        cell, ...
        cells ...
    )
    %GET_NEIGHOUR_CELLS Classifies the neighbours of a cell as left or right.
    %   neighbours = GET_NEIGHOUR_CELLS(cell, cells) inspects the connectivity
    %   field of `cell` and, for each connected neighbour, decides whether it
    %   lies to the left (smaller x centroid) or to the right (larger x centroid)
    %   of the current cell.
    %
    %   Inputs:
    %   -------
    %   cell  : struct  - Current cell with 'centroid' and 'connectivity' fields.
    %   cells : struct array (1 x N) - Full array of mesh cells.
    %
    %   Outputs:
    %   --------
    %   neighbours : struct
    %     .left  - Index of the left neighbour cell.
    %     .right - Index of the right neighbour cell.

    neighbours_indices = cell.connectivity;

    for i = 1:length(neighbours_indices)
        neighbour_index = neighbours_indices(i);
        neighbour_cell = cells(neighbour_index);
        centroid_link_vec = cell.centroid - neighbour_cell.centroid;

        if (centroid_link_vec(1) < 0)
            neighbours.left = neighbour_index;
        else
            neighbours.right = neighbour_index;
        end

    end
