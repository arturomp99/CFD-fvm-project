classdef FilePaths

    properties (Constant)
        NODES = 'Meshes/coords_1D_34nodes.dat';
        CELLS = 'Meshes/cells_1D_34nodes_allsquares_sorted.dat';
        BOUNDARY_CONDITIONS = { ...
                                   'Meshes/bc_bottom_1D_34nodes.dat', ...
                                   'Meshes/bc_left_1D_34nodes.dat', ...
                                   'Meshes/bc_right_1D_34nodes.dat', ...
                                   'Meshes/bc_top_1D_34nodes.dat', ...
                                   'Meshes/bc_whole_contour_1D_34nodes.dat' ...
                               };

        SKIP_NODES_FILE_HEADER = false;
        SKIP_NODES_FILE_FOOTER = false;
        SKIP_CELLS_FILE_HEADER = false;
        SKIP_CELLS_FILE_FOOTER = false;
    end

end
