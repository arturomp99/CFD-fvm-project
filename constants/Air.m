classdef Air

    properties (Constant)
        % Propiedades termodinámicas fundamentales
        GAMMA = 1.4; % Relación de calores específicos [-]
        C_P = 1004.5; % Calor específico a presión constante [J/(kg·K)]
        C_V = 1004.5/1.4; % Calor específico a volumen constante [J/(kg·K)]
        R = 1004.5 - (1004.5/1.4); % Constante específica del gas [J/(kg·K)]

        % Condiciones atmosféricas
        SEA_LEVEL_PRESSURE = 101325; % Presión estándar [Pa]
        SEA_LEVEL_DENSITY = 1.225; % Densidad estándar [kg/m³]
        SEA_LEVEL_TEMPERATURE = 288.15; % Temperatura estándar [K]
    end

end
