function [cell_nodes, volume, centroid, areas, normals] = process_cell(...
    cell, all_nodes_data)
% DOCUMENTAR DOCUMENTAR
    cell_nodes = [];
    for node_index=cell
        node_coordinates = all_nodes_data{node_index};
        cell_nodes = [...
            cell_nodes;
            node_coordinates];
    end
    volume = polyarea(cell_nodes(:, 1), cell_nodes(:, 2));
    centroid = cell_centroid(cell_nodes);
    areas = cell_areas(cell_nodes);
    normals = cell_normals(cell_nodes, centroid);
end
