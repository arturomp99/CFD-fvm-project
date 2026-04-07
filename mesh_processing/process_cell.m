function cell = process_cell(cell, all_nodes_data)
    %PROCESS_CELL Compute geometric data for a single mesh cell.
    %   [cell_nodes, volume, centroid, areas, normals] = PROCESS_CELL(cell,
    %   all_nodes_data) builds the geometry of one cell from raw node data.
    %
    %   Inputs:
    %       cell           - vector of node indices defining the cell polygon.
    %       all_nodes_data - cell array of node coordinates loaded from the
    %                        mesh node file.
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
    cell.nodes = [];

    for node_index = cell
        node_coordinates = all_nodes_data{node_index};
        cell.nodes = [
                      cell.nodes;
                      node_coordinates
                      ];
    end

    cell.volume = polyarea(cell.nodes(:, 1), cell.nodes(:, 2));
    cell.centroid = cell_centroid(cell.nodes);
    cell.areas = cell_areas(cell.nodes);
    cell.normals = cell_normals(cell.nodes, cell.centroid);
end
