function cell_velocity = get_cell_velocity(cell_index, num_cells, state)
    %GET_CELL_VELOCITY Computes the flow velocity for a given cell from the state vector.
    %   velocity = GET_CELL_VELOCITY(cell_index, num_cells, state) extracts the
    %   density and momentum from the state vector for the specified cell and
    %   returns the velocity as momentum / density.
    %
    %   Inputs:
    %   -------
    %   cell_index : integer
    %     1-based index of the target cell.
    %   num_cells : integer
    %     Total number of cells N. Used to locate momentum in the state vector.
    %   state : column vector (3*N x 1)
    %     Concatenated state vector [density (N); momentum (N); total energy (N)].
    %
    %   Outputs:
    %   --------
    %   cell_velocity : double
    %     Flow velocity at the cell centre. [m/s]

    cell_density = state(cell_index);
    cell_momentum = state(num_cells + cell_index);
    cell_velocity = cell_momentum / cell_density;

end
