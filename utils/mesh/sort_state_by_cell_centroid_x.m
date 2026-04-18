function [x_sorted, perm, rho_sorted, rhou_sorted, E_sorted, cells_sorted] = sort_state_by_cell_centroid_x(state, cells)
    %SORT_STATE_BY_CELL_CENTROID_X Sorts the state and cells by cell centroid x-coordinate.
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     Stacked conservative state vector [density(N); momentum(N); energy(N)].
    %   cells : struct array (1 x N)
    %     Mesh cells with a .centroid field [x, y].
    %
    %   Outputs:
    %   --------
    %   x_sorted : column vector (N x 1)
    %     Sorted x-coordinates of the cell centroids.
    %   perm : column vector (N x 1)
    %     Permutation indices that sort the cells by x.
    %   rho_sorted : column vector (N x 1)
    %     Density entries reordered to follow x_sorted.
    %   rhou_sorted : column vector (N x 1)
    %     Momentum entries reordered to follow x_sorted.
    %   E_sorted : column vector (N x 1)
    %     Energy entries reordered to follow x_sorted.
    %   cells_sorted : struct array (1 x N)
    %     Cells reordered to match x_sorted.

    x = get_cell_centroid_x(cells);
    [x_sorted, perm] = sort(x);

    [rho, rhou, E] = state_vec2states(state);
    rho_sorted = rho(perm);
    rhou_sorted = rhou(perm);
    E_sorted = E(perm);

    cells_sorted = cells(perm);
end
