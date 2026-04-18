function nodes = cell_nodes(cell_nodes_indices, all_nodes_data)
    %CELL_NODES Return node coordinates [x,y] for one cell.
    %
    %   nodes = CELL_NODES(cell_nodes_indices, all_nodes_data)
    %   returns an N-by-2 array with node coordinates.

    num_nodes = length(cell_nodes_indices);
    nodes = zeros(num_nodes, 2);

    for k = 1:num_nodes
        node_index = cell_nodes_indices(k);

        if node_index < 1 || node_index > length(all_nodes_data)
            error('Node index %d exceeds loaded nodes (%d).', ...
                node_index, length(all_nodes_data));
        end

        node_coordinates = all_nodes_data{node_index};

        % Keep only finite numeric values
        node_coordinates = node_coordinates(isfinite(node_coordinates));

        if numel(node_coordinates) < 2
            error('Node %d does not contain at least two coordinates.', node_index);
        end

        nodes(k, :) = node_coordinates(1:2);
    end

end
