function is_on_boundary = is_face_on_boundary(node1_idx, node2_idx, bc_nodes)
    %IS_FACE_ON_BOUNDARY Check if a face lies on a boundary surface.
    %   is_on_boundary = IS_FACE_ON_BOUNDARY(node1_idx, node2_idx, bc_nodes)
    %   returns true if both nodes of the face belong to the boundary surface.
    %
    %   Inputs:
    %   -------
    %   node1_idx : integer
    %       Index of the first node of the face.
    %   node2_idx : integer
    %       Index of the second node of the face.
    %   bc_nodes : vector
    %       Node indices belonging to the boundary surface.
    %
    %   Outputs:
    %   --------
    %   is_on_boundary : logical
    %       True if both node1_idx and node2_idx are in bc_nodes.

    is_on_boundary = ismember(node1_idx, bc_nodes) && ismember(node2_idx, bc_nodes);
end
