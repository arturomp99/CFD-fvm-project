% Mensaje de finalización del procesado de malla
% Detiene timer e informa tiempo transcurrido
function finishing_mesh_processing_msg()
    execution_time = toc;
    fprintf("Finalizada en %f.5 s", execution_time);
end
