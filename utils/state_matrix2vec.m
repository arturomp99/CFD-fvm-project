function state_vec = state_matrix2vec(state_matrix)
    %STATE_MATRIX2VEC Converts an N x 3 state matrix into a flat column vector.
    %   state_vec = STATE_MATRIX2VEC(state_matrix) stacks the three columns of
    %   the per-cell state matrix into the concatenated format
    %   [density (N); momentum (N); total energy (N)].
    %
    %   Inputs:
    %   -------
    %   state_matrix : matrix (N x 3)
    %     Per-cell state variables [density, momentum, total energy].
    %
    %   Outputs:
    %   --------
    %   state_vec : column vector (3*N x 1)
    %     Concatenated state vector.

    state_vec = [
                 state_matrix(:, 1);
                 state_matrix(:, 2);
                 state_matrix(:, 3);
                 ];
end
