function temperature = get_temperature(density, pressure)
    %GET_TEMPERATURE Calcula temperatura usando la ecuación de estado de gas ideal
    %
    %   Aplica la ecuación de estado de gas ideal para obtener temperatura
    %   a partir de densidad y presión.
    %   
    %   ECUACIÓN DE ESTADO:
    %   ==================
    %   T = p/(ρR)
    %   
    %   Input
    %   ---------------------  
    %   density : double
    %       Densidad del fluido ρ [kg/m³]
    %       
    %   pressure : double
    %       Presión estática p [Pa]
    %
    %   Output
    %   --------
    %   temperature : double
    %       Temperatura absoluta T [K]
    %       

    temperature = pressure / Air.R / density;
end
