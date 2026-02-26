function [output_mesh] = compute_connectivity(input_mesh)
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

    output_mesh = input_mesh;
    output_mesh.connectivity = {};
    for i=1:length(output_mesh.volume)
        % DUMMY! TO BE COMPLETED!
        output_mesh.connectivity{i} = [0];
    end
    % Also: COMPUTE THE CONNECTIVITY TO THE BOUNDARY CONDITIONS!!!
end
