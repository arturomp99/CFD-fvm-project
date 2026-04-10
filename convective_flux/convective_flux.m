function [A, b] = convective_flux(state, cells)
    state_matrix = state_vec2matrix(state);
    [A, b] = Config.CONVECTIVE_FLUX_INTERPOLATOR(state, cells);
end
