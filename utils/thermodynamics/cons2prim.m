function [rho, u, p, E] = cons2prim(U, gamma)
    %CONS2PRIM Convert conserved variables to primitive variables.

    rho = max(U(1), 1e-10);
    rhou = U(2);
    E = max(U(3), 1e-10);

    u = rhou / rho;
    p = (gamma - 1) * (E - 0.5 * rhou ^ 2 / rho);
    p = max(p, 1e-10);
end
