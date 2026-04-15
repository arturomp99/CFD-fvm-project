function [node1_idx, node2_idx] = get_face_node_indices(cell_node_indices, face_idx)
    %GET_FACE_NODE_INDICES Get the two node indices that form a face.
    %   [node1_idx, node2_idx] = GET_FACE_NODE_INDICES(cell_node_indices, face_idx)
    %   returns the indices of the two nodes that define the specified face
    %   of a cell polygon.
    %
    %   Inputs:
    %   -------
    %   cell_node_indices : vector
    %       Node indices defining the cell polygon (in order).
    %   face_idx : integer
    %       Index of the face (1 to num_faces). Faces are defined by 
    %       consecutive node pairs, with the last face connecting the 
    %       last node back to the first.
    %
    %   Outputs:
    %   --------
    %   node1_idx : integer
    %       Index of the first node of the face.
    %   node2_idx : integer
    %       Index of the second node of the face.
    %
    %   Example:
    %   --------
    %   For a quadrilateral cell with node indices [10, 20, 30, 40]:
    %     Face 1: nodes 10 and 20
    %     Face 2: nodes 20 and 30
    %     Face 3: nodes 30 and 40
    %     Face 4: nodes 40 and 10 (wraps around)

    num_nodes = length(cell_node_indices);
    
    node1_idx = cell_node_indices(face_idx);
    next_idx = get_next_vertex_index(face_idx, num_nodes);
    node2_idx = cell_node_indices(next_idx);
end
