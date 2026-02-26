function mesh_data = process_mesh_geometry(...
    nodes_data, cells_data, bcs_data)
% DOCUMENTAR DOCUMENTAR
    mesh_data = struct();
    for i=1:length(cells_data)
        try
            % Process the data of one single cell and save it in mesh_data.
            [cell_nodes, volume, centroid, area, normal] = ...
                process_cell(cells_data{i}, nodes_data);
            mesh_data.cell_nodes{i} = cell_nodes;
            mesh_data.volume(i) = volume;
            mesh_data.centroid(i, :) = centroid;
            mesh_data.area{i} = area;
            mesh_data.normal{i} = normal;
        catch exception
            msg = sprintf('In cell number %d:\n  %s', ...
                i, exception.message);
            throw(MException('Simulator:MeshProcessorError', msg));
        end
    end
end