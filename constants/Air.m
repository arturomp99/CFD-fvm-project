classdef Air

    properties (Constant)
        GAMMA = 1.4; % Relación de calores específicos
        C_P = 1004.5; % Calor específico presión constante [J/(kg·K)]
        C_V = 1004.5/1.4; % Calor específico volumen constante [J/(kg·K)]
        R = 1004.5 - (1004.5/1.4); % Constante de los gases ideales
    end

end
