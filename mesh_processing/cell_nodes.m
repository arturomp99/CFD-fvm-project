function nodes = cell_nodes(cell_nodes_indices, all_nodes_data)
    %CELL_NODES Gathers the coordinate array for a cell from global node data.
    %   nodes = CELL_NODES(cell_nodes_indices, all_nodes_data) collects the
    %   coordinates of each node listed in cell_nodes_indices and returns them
    %   as an N x 2 matrix.
    %
    %   Inputs:
    %   -------
    %   cell_nodes_indices : row vector of integers
    %     1-based indices of the nodes that form this cell.
    %   all_nodes_data : cell array
    %     Global node coordinate data loaded from the mesh nodes file.
    %     all_nodes_data{i} is a 1x2 vector [x, y] for node i. [m]
    %
    %   Outputs:
    %   --------
    %   nodes : N x 2 array
    %     Coordinates [x, y] of the cell nodes in traversal order. [m]

    nodes = [];

    for node_index = cell_nodes_indices
        node_coordinates = all_nodes_data{node_index};
        nodes = [
                 nodes;
                 node_coordinates
                 ];
    end

end
