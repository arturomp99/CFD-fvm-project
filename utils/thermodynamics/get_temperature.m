function temperature = get_temperature(density, pressure)
    %GET_TEMPERATURE Computes temperature from the ideal gas law.
    %   T = GET_TEMPERATURE(density, pressure) applies the ideal gas relation
    %   p = rho * R * T to obtain T = p / (R * rho).
    %
    %   Inputs:
    %   -------
    %   density  : double - Mass density rho. [kg/m^3]
    %   pressure : double - Static pressure p. [Pa]
    %
    %   Outputs:
    %   --------
    %   temperature : double - Static temperature T. [K]

    temperature = pressure / Air.R / density;
end
