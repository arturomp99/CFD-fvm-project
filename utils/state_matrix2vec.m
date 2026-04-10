function state_vec = state_matrix2vec(state_matrix)
    state_vec = [
                 state_matrix(:, 1);
                 state_matrix(:, 2);
                 state_matrix(:, 3);
                 ];
end
