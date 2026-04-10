function neighbours = get_neighour_cells( ...
        cell, ...
        cells ...
    )

    neighbours_indices = cell.connectivity;

    for i = 1:neighbours_indices
        neighbour_index = neighbours_indices(i);
        neighbour_cell = cells(neighbour_index);
        centroid_link_vec = cell.centroid - neighbour_cell.centroid;

        if (centroid_link_vec < 0)
            neighbours.left = neighbour_index;
        else
            neighbours.right = neighbour_index;
        end

    end
