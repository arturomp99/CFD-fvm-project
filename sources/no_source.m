function [A, b] = no_source(centroids_x)
    num_cells = length(centroids_x);
    [A, b] = initialize_A_b(num_cells);
end
