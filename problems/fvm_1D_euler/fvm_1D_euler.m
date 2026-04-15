function [A, b] = fvm_1D_euler( ...
        state, ...
        cells, ...
        boundary_info ...
    )
    %FVM_1D_EULER Computes the spatial discretisation matrices for the 1D Euler equations.
    %   [A, b] = FVM_1D_EULER(state, cells, boundary_info) assembles the system matrices for the
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
    %   boundary_info : struct (optional)
    %     Structure containing boundary condition configuration:
    %       .boundary_types - cell array of BC types ('open' or 'wall') for each surface
    %
    %   Outputs:
    %   --------
    %   A : matrix (3N x 3N)
    %     Spatial operator matrix such that dw/dt = A*w + b.
    %   b : column vector (3N x 1)
    %     Independent terms vector (source / boundary contributions).

    if nargin < 3
        boundary_info = struct('boundary_types', {{}});
    end

    [A, b] = convective_flux(state, cells, boundary_info);

end
