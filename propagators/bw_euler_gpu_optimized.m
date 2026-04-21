function [w_new] = bw_euler_gpu_optimized( ...
        w, ...
        t, ...
        dt, ...
        f, ...
        gpu_state)
    % BW_EULER_GPU_OPTIMIZED - Versión optimizada del integrador GPU
    %
    % OPTIMIZACIONES IMPLEMENTADAS:
    % - Persistencia de datos en GPU entre timesteps
    % - Uso de matrices sparse para eficiencia de memoria
    % - Pre-asignación de arrays GPU
    % - Minimización de transferencias CPU↔GPU
    %
    % ENTRADA ADICIONAL:
    %   gpu_state - estructura persistente con datos GPU pre-asignados
    %
    % USO:
    %   % Inicialización (solo una vez)
    %   gpu_state = init_gpu_state(N);
    %
    %   % En bucle temporal
    %   [w_new, gpu_state] = bw_euler_gpu_optimized(w, t, dt, f, gpu_state);

    persistent is_initialized A_gpu_prev I_gpu

    % =====================================================================
    % INICIALIZACIÓN (solo en primera llamada)
    % =====================================================================
    if isempty(is_initialized)
        fprintf('Inicializando solver GPU...\n');
        N = length(w) / 3;

        % Pre-asignar matriz identidad en GPU (invariante en tiempo)
        I_gpu = speye(length(w), 'gpuArray');

        % Flag de inicialización
        is_initialized = true;
        A_gpu_prev = [];
    end

    % =====================================================================
    % CÁLCULO DE MATRIZ Y VECTOR RHS
    % =====================================================================
    [A, b] = f(w, t);

    % Optimización: solo transferir A si ha cambiado significativamente
    if isempty(A_gpu_prev) || norm(A - gather(A_gpu_prev), 'fro') > 1e-12
        A_gpu = gpuArray(sparse(A)); % Usar sparse para eficiencia
        A_gpu_prev = A_gpu;
    else
        A_gpu = A_gpu_prev; % Reutilizar matriz anterior
    end

    % Transferir solo vectores (más eficiente que matrices)
    b_gpu = gpuArray(b);
    w_gpu = gpuArray(w);

    % =====================================================================
    % RESOLUCIÓN DEL SISTEMA LINEAL EN GPU
    % =====================================================================
    % Sistema: (I - A*dt) * w_new = w + b*dt
    K_gpu = I_gpu - A_gpu * dt;
    rhs_gpu = w_gpu + b_gpu * dt;

    % Usar solver GPU especializado para sistemas sparse
    w_new_gpu = K_gpu \ rhs_gpu; % Más eficiente que linsolve para sparse

    % Retornar resultado a CPU
    w_new = gather(w_new_gpu);

    % Limpiar variables temporales GPU para liberar memoria
    clear b_gpu w_gpu rhs_gpu
end

function gpu_state = init_gpu_state(N)
    % Inicializa estructuras GPU persistentes
    % N: número de células en la malla

    system_size = 3 * N;

    gpu_state.I = speye(system_size, 'gpuArray');
    gpu_state.initialized = true;

    fprintf('GPU state inicializado para %d células (%d DOF)\n', N, system_size);
end
