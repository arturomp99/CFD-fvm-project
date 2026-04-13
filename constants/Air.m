classdef Air

    properties (Constant)
        GAMMA = 1.4; % Relación de calores específicos
        C_P = 1004.5; % Calor específico presión constante [J/(kg·K)]
        C_V = 1004.5/1.4; % Calor específico volumen constante [J/(kg·K)]
        R = 1004.5 - (1004.5/1.4); % Constante de los gases ideales

        % Sea level
        SEA_LEVEL_PRESSURE = 101325; % Pa
        SEA_LEVEL_DENSITY = 1.225; % kg/m^3
        SEA_LEVEL_TEMPERATURE = 288.15; % K
    end

end
