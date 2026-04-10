function [A, b] = initialize_A_b(num_cells)
    A = zeros(num_cells * 3, num_cells * 3);
    b = zeros(num_cells * 3, 1);
end
