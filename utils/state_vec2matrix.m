function state_matrix = state_vec2matrix(state_vec)
    %STATE_VEC2MATRIX Reshapes the flat state vector into an N x 3 matrix.
    %   state_matrix = STATE_VEC2MATRIX(state_vec) converts the concatenated
    %   state vector [density (N); momentum (N); total energy (N)] into an
    %   N x 3 matrix where each row is one cell and columns are
    %   [density, momentum, total energy].
    %
    %   Inputs:
    %   -------
    %   state_vec : column vector (3*N x 1)
    %     Concatenated state vector.
    %
    %   Outputs:
    %   --------
    %   state_matrix : matrix (N x 3)
    %     Per-cell state variables [density, momentum, total energy].

    state_matrix = reshape(state_vec, [], 3);
end
