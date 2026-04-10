function output_cells = compute_connectivity(input_cells)
    %COMPUTE_CONNECTIVITY Build adjacency information for mesh cells.
    %   output_cells = COMPUTE_CONNECTIVITY(input_cells) computes which cells
    %   are neighbors (share a common face) for each cell in the mesh.
    %
    %   Input:
    %       input_cells - array of cell structures, each with a 'faces' field
    %                     containing face definitions (N-by-2-by-2 arrays).
    %
    %   Output:
    %       output_cells - array of cell structures with added 'connectivity' field.
    %                      For each cell, connectivity(i) contains the indices
    %                      of all neighboring cells that share a face with cell i.
    %
    %   The function uses SHARE_FACE to determine if two cells are neighbors,
    %   then stores the neighbor indices in each cell's connectivity field.
    %
    %   Example:
    %       For a mesh with 3 cells where cells 1 and 2 share a face,
    %       output_cells(1).connectivity = [2]
    %       output_cells(2).connectivity = [1]
    %       output_cells(3).connectivity = []  % isolated cell

    output_cells = input_cells;
    num_cells = length(input_cells);

    for cell_index = 1:num_cells
        output_cells(cell_index).connectivity = [];
        cell = input_cells(cell_index);

        % Loop through all OTHER cells (skip self-comparison)
        for other_cell_index = 1:num_cells

            if other_cell_index ~= cell_index
                other_cell = input_cells(other_cell_index);

                if share_face(cell, other_cell)
                    output_cells(cell_index).connectivity = ...
                        [output_cells(cell_index).connectivity, other_cell_index];
                end

            end

        end

    end

end
