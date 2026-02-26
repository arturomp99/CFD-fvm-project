function [A, b] = dummy_fluxes(w, t, mesh, cp, cv, rho, k, ...
        speed)
%DUMMY_FLUXES Summary of this function goes here
%   Detailed explanation goes here
    A = zeros(length(mesh.volume));
    b = zeros(length(mesh.volume), 1);
end
