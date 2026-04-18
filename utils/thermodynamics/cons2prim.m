function [rho, u, p, E] = cons2prim(U, gamma)
    %CONS2PRIM Convierte variables conservadas del vector estado en primitivas.
    %   Input
    %   -------
    %   U : vector (3 x 1)
    %     Variables conservadas [rho; rho*u; E].
    %   gamma : escalar
    %     Relación de calores específicos.
    %
    %   Salidas:
    %   -------
    %   rho : escalar
    %     Densidad.
    %   u : escalar
    %     Velocidad.
    %   p : escalar
    %     Presión.
    %   E : escalar
    %     Energía total.

    rho = max(U(1), 1e-10);
    rhou = U(2);
    E = max(U(3), 1e-10);
    u = rhou / rho;
    p = (gamma - 1) * (E - 0.5 * rhou ^ 2 / rho);
    p = max(p, 1e-10);
end
