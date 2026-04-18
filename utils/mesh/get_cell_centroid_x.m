function x = get_cell_centroid_x(cells)
    %GET_CELL_CENTROID_X Extrae las coordenadas x de los centroides de las celdas
    %
    %   Inputs:
    %   -------
    %   cells : struct array (1 x N)
    %     Vector de celdas que tienen una propiedad .centroid definida como [x, y].
    %
    %   Outputs:
    %   --------
    %   x : column vector (N x 1)
    %     Coordenadas x de las celdas

    if isempty(cells)
        x = zeros(0, 1);
        return;
    end

    centroids = reshape([cells.centroid], 2, [])';
    x = centroids(:, 1);
end
