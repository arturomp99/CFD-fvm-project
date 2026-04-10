function cell_velocity = get_cell_velocity(cell_index, num_cells, state)
    cell_density = state(cell_index);
    cell_momentum = state(num_cells + cell_index);
    cell_velocity = cell_momentum / cell_density;

end
