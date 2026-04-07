function cell_data = process_cell(cell_nodes_indices, all_nodes_data)
    %PROCESS_CELL Compute geometric data for a single mesh cell.
    %   [cell_nodes, volume, centroid, areas, normals] = PROCESS_CELL(cell,
    %   all_nodes_data) builds the geometry of one cell from raw node data.
    %
    %   Inputs:
    %       cell_nodes_indices  - vector of node indices defining the cell polygon.
    %       all_nodes_data      - cell array of node coordinates loaded from the
    %                             mesh node file.
    %
    %   Outputs:
    %       cell - with properties
    %           cell_nodes - N-by-D array of the coordinates for the cell nodes.
    %           volume     - polygon area computed from the cell node coordinates.
    %           centroid   - geometric centroid of the cell.
    %           areas      - per-face/edge areas for the cell.
    %           normals    - outward normals for each face/edge.
    %
    %   The function reconstructs the cell node coocrdinates from the provided
    %   node indices, then computes the cell area, centroid, face areas, and
    %   normals using helper functions.

    cell_data.nodes = [];

    for node_index = cell_nodes_indices
        node_coordinates = all_nodes_data{node_index};
        cell_data.nodes = [
                           cell_data.nodes;
                           node_coordinates
                           ];
    end

    cell_data.volume = polyarea(cell_data.nodes(:, 1), cell_data.nodes(:, 2));
    cell_data.centroid = cell_centroid(cell_data.nodes);
    cell_data.areas = cell_areas(cell_data.nodes);
    cell_data.normals = cell_normals(cell_data.nodes, cell_data.centroid);
end
