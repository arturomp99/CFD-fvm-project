function mesh_data = process_mesh_geometry( ...
        nodes_data, ...
        cells_data, ...
        bcs_data ...
    )
    % PROCESS_MESH_GEOMETRY Build geometric and topological data for the mesh.
    %   mesh_data = PROCESS_MESH_GEOMETRY(nodes_data, cells_data, bcs_data)
    %   converts raw node, cell, and boundary condition datasets into a
    %   structured mesh representation used by the solver.
    %
    %   Inputs:
    %       nodes_data - raw node coordinates loaded from the mesh node file.
    %       cells_data - cell connectivity data loaded from the mesh cell file.
    %       bcs_data   - boundary condition definitions loaded from BC files.
    %
    %   Output:
    %       mesh_data  - structure containing per-cell geometry fields:
    %           cell_nodes - node indices for each cell
    %           volume     - cell volume for each cell
    %           centroid   - cell centroid coordinates
    %           area       - face/edge areas for each cell
    %           normal     - outward normals for each face/edge
    %
    %   The function iterates over each cell in cells_data, computes the
    %   required geometric quantities with process_cell, and packages the
    %   results in mesh_data. Any error in a cell is rethrown with a
    %   descriptive message that includes the offending cell index.

    mesh_data = struct();
    num_cells = length(cells_data);

    for i = 1:num_cells

        try
            % Process the data of one single cell and save it in mesh_data.
            cell_data = process_cell(cells_data{i}, nodes_data);
            mesh_data.cell_nodes{i} = cell_data.nodes;
            mesh_data.volume(i) = cell_data.volume;
            mesh_data.centroid(i, :) = cell_data.centroid;
            mesh_data.area{i} = cell_data.areas;
            mesh_data.normal{i} = cell_data.normals;
        catch exception
            msg = sprintf('In cell number %d:\n  %s', ...
                i, exception.message);
            throw(MException('Simulator:MeshProcessorError', msg));
        end

    end

    mesh_data.num_cells = num_cells;
end
