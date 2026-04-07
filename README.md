# Description

# Installation

# How to run

- Populate the Meshes directory with `.dat` files for
  - nodes
  - cells
  - boundary conditions
- Modify the paths from `FilePaths.m` in the constants subfolder. Example

```
NODES = 'Meshes/coords_1D_10nodes.dat';
        CELLS = 'Meshes/cells_1D_10nodes_allsquares_shuffled.dat';
        BOUNDARY_CONDITIONS = { ...
                                   'Meshes/bc_bottom_1D_10nodes.dat', ...
                                   'Meshes/bc_left_1D_10nodes.dat', ...
                                   'Meshes/bc_right_1D_10nodes.dat', ...
                                   'Meshes/bc_top_1D_10nodes.dat', ...
                                   'Meshes/bc_whole_contour_1D_10nodes.dat', ...
                               };
```
