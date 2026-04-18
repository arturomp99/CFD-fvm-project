% Mensaje de finalización de la resolución numérica
% Detiene timer e informa tiempo total de simulación
function finishing_solver_msg()
    execution_time = toc;
    fprintf("Finalizada en %f.5 s", execution_time);
end
