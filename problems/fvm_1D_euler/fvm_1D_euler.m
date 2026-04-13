function [A, b] = fvm_1D_euler( ...
        state, ...
        cells ...
    )
    %FVM_1D_EULER Computes the spatial discretisation matrices for the 1D Euler equations.
    %   [A, b] = FVM_1D_EULER(state, cells) assembles the system matrices for the
    %   finite volume discretisation of the inviscid 1D Euler equations by
    %   computing the convective fluxes through each cell face.
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     Concatenated state vector [density; momentum; total energy] for all N cells.
    %   cells : struct array (1 x N)
    %     Mesh cell structures as produced by mesh_processor, containing geometry
    %     and connectivity information.
    %
    %   Outputs:
    %   --------
    %   A : matrix (3N x 3N)
    %     Spatial operator matrix such that dw/dt = A*w + b.
    %   b : column vector (3N x 1)
    %     Independent terms vector (source / boundary contributions).

    [A, b] = convective_flux(state, cells);

end
