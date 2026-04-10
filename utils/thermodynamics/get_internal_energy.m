function energy = get_internal_energy(density, pressure, velocity)
    %GET_INTERNAL_ENERGY Computes the total specific energy (internal + kinetic) per unit volume.
    %   energy = GET_INTERNAL_ENERGY(density, pressure, velocity) calculates the
    %   total energy E = rho * (c_v * T + v^2/2), where the temperature is derived
    %   from the ideal gas law.
    %
    %   Inputs:
    %   -------
    %   density  : double - Mass density rho. [kg/m^3]
    %   pressure : double - Static pressure p. [Pa]
    %   velocity : double - Flow velocity magnitude v. [m/s]
    %
    %   Outputs:
    %   --------
    %   energy : double - Total energy per unit volume E = rho*(c_v*T + v^2/2). [J/m^3]

    energy = density * (Air.C_V * get_temperature(density, pressure) + velocity ^ 2/2);
end
