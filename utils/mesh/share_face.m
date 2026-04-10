function boolean = share_face(cell1, cell2)
    %SHARE_FACE Check if two cells share a common face.
    %   boolean = SHARE_FACE(cell1, cell2) compares the faces of two cells
    %   and returns true if they share at least one face (same two nodes),
    %   false otherwise.
    %
    %   Inputs:
    %       cell1, cell2 - structures with 'faces' field (N-by-2-by-2 arrays)
    %                      where faces are defined by two nodes with (x,y) coordinates
    %
    %   Output:
    %       boolean - true if cells share a face, false otherwise

    tolerance = 1e-10;

    faces1 = cell1.faces; % N-by-2-by-2 array
    faces2 = cell2.faces; % M-by-2-by-2 array

    % Loop through all faces of cell1
    for i = 1:size(faces1, 1)
        face1_node1 = squeeze(faces1(i, 1, :)); % First node of face i
        face1_node2 = squeeze(faces1(i, 2, :)); % Second node of face i

        % Loop through all faces of cell2
        for j = 1:size(faces2, 1)
            face2_node1 = squeeze(faces2(j, 1, :)); % First node of face j
            face2_node2 = squeeze(faces2(j, 2, :)); % Second node of face j

            % Check if faces are the same (same nodes, considering both directions)
            same_direction = (norm(face1_node1 - face2_node1) < tolerance) && ...
                (norm(face1_node2 - face2_node2) < tolerance);

            opposite_direction = (norm(face1_node1 - face2_node2) < tolerance) && ...
                (norm(face1_node2 - face2_node1) < tolerance);

            if same_direction || opposite_direction
                boolean = true;
                return;
            end

        end

    end

    boolean = false;
end
