function [A, b] = convective_flux(state, cells, boundary_info)
    %CONVECTIVE_FLUX Dispatches convective flux computation to the configured interpolator.
    %   [A, b] = CONVECTIVE_FLUX(state, cells, boundary_info) delegates the assembly of the convective
    %   flux matrices to whichever interpolator is set in Config.CONVECTIVE_FLUX_INTERPOLATOR
    %   (e.g. upwind_interpolator or linear_interpolator).
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     Concatenated state vector [density; momentum; total energy] for all N cells.
    %   cells : struct array (1 x N)
    %     Mesh cell structures containing geometry and connectivity information.
    %   boundary_info : struct (optional)
    %     Structure containing boundary condition configuration.
    %
    %   Outputs:
    %   --------
    %   A : matrix (3N x 3N)
    %     Convective flux operator matrix.
    %   b : column vector (3N x 1)
    %     Independent terms vector.

    if nargin < 3
        boundary_info = struct('boundary_types', {{}});
    end

    state_matrix = state_vec2matrix(state);
    [A, b] = Config.CONVECTIVE_FLUX_INTERPOLATOR(state, cells, boundary_info);
end
