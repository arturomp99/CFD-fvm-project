function [nodes_data, cells_data, bcs_data] = mesh_data_loader( ...
        nodes_file, ...
        cells_file, ...
        bc_files ...
    )
    % MESH_DATA_LOADER Loads mesh data from mesh data files.
    %
    %   Inputs
    %   ------
    %   nodes_file : string
    %     Name of the file that contains the coordinates of the nodes.
    %   cells_data : string
    %     Name of the file that contains the nodes forming each cell.
    %   bcs_data : cell(string)
    %     List of names of files that contain the nodes forming each boundary
    %     condition surface.
    %
    %   Outputs
    %   -------
    %   nodes_data : cell(double)
    %     Cell array with the coordinates of each node. The coordinates are
    %     expressed in m.
    %   cells_data : cell(double)
    %     Cell array with the node indices forming each cell.
    %   bcs_data : cell(double)
    %     Cell array with the nodes forming each boundary condition surface.

    nodes_data = data_loader( ...
        nodes_file, ...
        "\t", ...
        FilePaths.SKIP_NODES_FILE_HEADER, ...
        FilePaths.SKIP_NODES_FILE_FOOTER ...
    );

    cells_data = data_loader( ...
        cells_file, ...
        "\t", ...
        FilePaths.SKIP_CELLS_FILE_HEADER, ...
        FilePaths.SKIP_CELLS_FILE_FOOTER ...
    );

    bcs_data = cell(length(bc_files), 1);

    for i = 1:length(bc_files)
        % For these BC files, do not skip the footer by default
        current_bc_data = data_loader(bc_files{i}, "\t", 0, 0);

        if isempty(current_bc_data)
            error('Boundary condition file %s is empty after loading.', bc_files{i});
        end

        % If a BC file has several lines, concatenate all node indices
        bc_nodes = [];

        for k = 1:length(current_bc_data)
            bc_nodes = [bc_nodes, current_bc_data{k}]; %#ok<AGROW>
        end

        bcs_data{i} = unique(bc_nodes, 'stable');
    end

end
