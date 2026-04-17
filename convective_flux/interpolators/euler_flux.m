function F = euler_flux_1d(U, gamma)
    %EULER_FLUX_1D Physical flux for the 1D Euler equations.

    [rho, u, p, E] = cons2prim(U, gamma);

    F = [
         rho * u;
         rho * u ^ 2 + p;
         u * (E + p)
         ];
end
