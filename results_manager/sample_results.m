function [new_results] = sample_results(w, t, old_results, ...
        sampling_interval)
    %SAMPLE_RESULTS Gestor de resultados con muestreo temporal a intervalos fijos
    %
    %   Almacena el estado de la simulación a intervalos regulares de tiempo,
    %   creando una matriz de resultados completa para post-procesamiento
    %   y visualización.
    %   
    %   1. Verifica si ha pasado suficiente tiempo desde el último muestreo
    %   2. Si Δt >= sampling_interval, agrega [t, w'] a la matriz de resultados
    %   3. Si no, retorna la matriz sin cambios
    %
    %   Input
    %   ---------------------
    %   w : double (3*N×1)
    %       Vector de estado actual [ρ; ρu; E]
    %       
    %   t : double
    %       Tiempo actual de simulación [s]
    %       
    %   old_results : double (M×(3*N+1))
    %       Matriz de resultados previos
    %       
    %   sampling_interval : double
    %       Intervalo de muestreo deseado [s]
    %
    %   Output
    %   --------
    %   new_results : double ((M o M+1)×(3*N+1))
    %       Matriz actualizada de resultados

    new_results = old_results;

    if (size(old_results, 1) > 0)
        last_t = old_results(end, 1);
    else
        last_t = -1E100;
    end

    if ((t - last_t) >= sampling_interval)
        new_results = [old_results; [t, w']];
    end

end
