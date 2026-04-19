function cells = process_mesh_geometry( ...
        nodes_data, ...
        cells_data, ...
        bcs_data ...
    )
    %PROCESS_MESH_GEOMETRY Procesamiento completo de geometría de malla y identificación de fronteras
    %
    %   Esta función es el núcleo del procesamiento geométrico, transformando
    %   datos en estructuras de células
    %   completas con toda la información geométrica necesaria para FVM.
    %
    %   1. Inicializa estructuras de células vacías
    %   2. Para cada célula:
    %      a) Extrae coordenadas de nodos
    %      b) Calcula centroide geométrico
    %      c) Calcula área/volumen usando polyarea()
    %      d) Define caras y calcula normales salientes
    %      e) Identifica caras que están en fronteras
    %   3. Construye información de superficies de frontera
    %
    %   ESTRUCTURA DE SALIDA:
    %   ====================

    %
    %   Input
    %   ---------------------
    %   nodes_data : double (M×2)
    %       Coordenadas [x, y] de todos los nodos de la malla
    %
    %   cells_data : cell array (N×1)
    %       Cada elemento es vector de índices de nodos que forman una célula
    %
    %   bcs_data : cell array
    %       Cada elemento es vector de índices de nodos en una superficie de frontera
    %
    %   Output
    %   --------
    %   cells : struct array (1×N)
    %       Estructuras de células completamente procesadas
    %       cells(i) contiene:
    %           - .nodes: coordenadas [x,y] de vértices
    %           - .centroid: coordenadas [cx, cy] del centro
    %           - .volume: área de la célula [m²]
    %           - .faces: definiciones de caras como endpoints
    %           - .area: longitudes de aristas [m]
    %           - .normals: vectores normales unitarios salientes
    %

    % initialize the cells vector
    num_cells = length(cells_data);
    cell_data = struct( ...
        'nodes', [], ...
        'faces', [], ...
        'volume', 0, ...
        'centroid', [0, 0], ...
        'area', [], ...
        'normals', [0, 0] ...
    );
    cells = repmat(cell_data, 1, num_cells);

    if Config.IS_MESH_PROCESSING_PARALLEL

        parfor cell_index = 1:num_cells
            % Process the data of one single cell (including boundary faces)
            cells(cell_index) = process_cell(cells_data{cell_index}, nodes_data, bcs_data);
        end

    else

        for cell_index = 1:num_cells
            % Process the data of one single cell (including boundary faces)
            cells(cell_index) = process_cell(cells_data{cell_index}, nodes_data, bcs_data);
        end

    end

end
