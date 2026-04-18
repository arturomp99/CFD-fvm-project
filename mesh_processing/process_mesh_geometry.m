function [cells, boundary_info] = process_mesh_geometry( ...
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
    %           - .boundary_faces: índices de caras en frontera
    %           - .boundary_surface_ids: mapeo a superficies específicas
    %           
    %   boundary_info : struct
    %       Información de superficies de frontera:
    %       - .surface_names: nombres de superficies
    %       - .surface_nodes: nodos por superficie
    %       

    % initialize the cells vector
    num_cells = length(cells_data);
    cell_data = struct( ...
        'nodes', [], ...
        'faces', [], ...
        'volume', 0, ...
        'centroid', [0, 0], ...
        'area', [], ...
        'normals', [0, 0], ...
        'boundary_faces', [], ...  % indices of faces that are boundaries
        'boundary_surface_ids', [] ... % which boundary surface each boundary face belongs to
    );
    cells = repmat(cell_data, 1, num_cells);

    % Process boundary surface node sets
    num_bc_surfaces = length(bcs_data);
    boundary_info = struct();
    boundary_info.surface_names = cell(num_bc_surfaces, 1);
    boundary_info.surface_nodes = cell(num_bc_surfaces, 1);
    
    for bc_idx = 1:num_bc_surfaces
        boundary_info.surface_names{bc_idx} = sprintf('surface_%d', bc_idx);
        boundary_info.surface_nodes{bc_idx} = bcs_data{bc_idx};
    end

    for cell_index = 1:num_cells
        % Process the data of one single cell (including boundary faces)
        cells(cell_index) = process_cell(cells_data{cell_index}, nodes_data, bcs_data);
    end

end
