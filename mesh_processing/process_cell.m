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

    % initialize the cell_data struct
    cell_data = struct( ...
        'nodes', [], ...
        'faces', [], ...
        'volume', 0, ...
        'centroid', [0, 0], ...
        'area', [], ...
        'normals', [0, 0] ...
    );

    cell_data.nodes = cell_nodes(cell_nodes_indices, all_nodes_data);
    cell_data.faces = cell_faces(cell_data.nodes);
    cell_data.volume = polyarea(cell_data.nodes(:, 1), cell_data.nodes(:, 2));
    cell_data.centroid = cell_centroid(cell_data.nodes);
    cell_data.area = cell_areas(cell_data.nodes);
    cell_data.normals = cell_normals(cell_data.nodes, cell_data.centroid);
end
