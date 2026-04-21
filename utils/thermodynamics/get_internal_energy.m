function energy = get_internal_energy(density, pressure, velocity)
    %GET_INTERNAL_ENERGY Calcula energía total por unidad de volumen para gas ideal
    %
    %   FÍSICA:
    %   =======
    %   E = ρe + ½ρu² = ρ(e + ½u²)
    %
    %   donde:
    %   - e = Cv*T: energía interna específica [J/kg]
    %   - ½u²: energía cinética específica [J/kg]
    %   - E: energía total por unidad de volumen [J/m³]
    %
    %   Input
    %   ---------------------
    %   density : double
    %       Densidad del fluido ρ [kg/m³]
    %
    %   pressure : double
    %       Presión estática p [Pa = N/m²]
    %
    %   velocity : double
    %       Magnitud de velocidad u [m/s]
    %
    %   Ouptut
    %   --------
    %   energy : double
    %       Energía total por unidad de volumen E [J/m³]
    %

    energy = density * (Air.C_V * get_temperature(density, pressure) + velocity ^ 2/2);
end
