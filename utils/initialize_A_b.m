function [A, b] = initialize_A_b(num_cells)
    %INITIALIZE_A_B Allocates zeroed system matrices for the FVM Euler problem.
    %   [A, b] = INITIALIZE_A_B(num_cells) creates a (3N x 3N) zero matrix A and
    %   a (3N x 1) zero vector b, sized for N cells with three conserved variables
    %   (density, momentum, total energy).
    %
    %   Inputs:
    %   -------
    %   num_cells : integer
    %     Number of mesh cells N.
    %
    %   Outputs:
    %   --------
    %   A : matrix (3N x 3N)
    %     Zero-initialised spatial operator matrix.
    %   b : column vector (3N x 1)
    %     Zero-initialised independent terms vector.

    A = zeros(num_cells * 3, num_cells * 3);
    b = zeros(num_cells * 3, 1);
end
