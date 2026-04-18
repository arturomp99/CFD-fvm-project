function [densities, momentums, energies] = state_vec2states(state_vec)
    %STATE_VEC2STATES Extrae las variables por cell del vector de estado
    %   Inputs:
    %   -------
    %   state_vec : column vector (3*N x 1)
    %     Vector estado concatenado que contiene densidad, momento y energía
    %
    %   Outputs:
    %   --------
    %   densities : column vector (N x 1)
    %     Densidad por celda.
    %   momentums : column vector (N x 1)
    %     Momento por celda.
    %   energies : column vector (N x 1)
    %     Energía total por celda.

    reshaped = reshape(state_vec, [], 3);
    densities = reshaped(:, 1);
    momentums = reshaped(:, 2);
    energies = reshaped(:, 3);
end
