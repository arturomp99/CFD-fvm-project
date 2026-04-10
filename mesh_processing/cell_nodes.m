function nodes = cell_nodes(cell_nodes_indices, all_nodes_data)
    nodes = [];

    for node_index = cell_nodes_indices
        node_coordinates = all_nodes_data{node_index};
        nodes = [
                 nodes;
                 node_coordinates
                 ];
    end

end
