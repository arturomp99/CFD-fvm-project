function connectivity = compute_connectivity(input_mesh)
    %COMPUTE_CONNECTIVITY Compute the connectivity of a mesh and save it.
    %   The connectivity is saved as a property of the mesh called
    %   connectivity. It is a cell array (one element per cell), in which each
    %   element has a vector of cell or boundary condition indices.
    %
    %   Input
    %   -----
    %   input_mesh
    %     Mesh structure to be processed.
    %
    %   Output
    %   ------
    %   output_mesh
    %     Mesh structure processed. Includes the connectivity information.

    num_cells = length(input_mesh.volume);
    connectivity = cell(num_cells, 1);

    % Build a map from edge (sorted node pair) to list of cells that have it
    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for cell_index = 1:num_cells
        cell_nodes = input_mesh.cell_nodes{cell_index};
        num_nodes = length(cell_nodes);
        connectivity{cell_index} = zeros(1, num_nodes); % One per face/edge

        for cell_node_index = 1:num_nodes
            next_node = get_next_vertex_index(cell_node_index, num_nodes);
            n1 = cell_nodes(cell_node_index);
            n2 = cell_nodes(next_node);
            % Sort nodes to create unique key
            key = sprintf('%d_%d', min(n1, n2), max(n1, n2));

            if ~isKey(edge_map, key)
                edge_map(key) = [];
            end

            edge_map(key) = [edge_map(key), i];
        end

    end

    % Now, assign connectivity
    for i = 1:num_cells
        cell_nodes = input_mesh.cell_nodes{i};
        num_nodes = length(cell_nodes);

        for node_index = 1:num_nodes
            next_node_index = get_next_vertex_index(node_index, num_nodes);
            node = cell_nodes(node_index);
            next_node = cell_nodes(next_node_index);
            key = sprintf('%d_%d', min(n1, n2), max(n1, n2));
            neighbors = edge_map(key);
            neighbors = neighbors(neighbors ~= i); % Remove self

            if isempty(neighbors)
                conn = 0; % Boundary
            else
                conn = neighbors(1); % Neighbor cell index
            end

            connectivity{i}(j) = conn;
        end

    end

end
