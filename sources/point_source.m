function [A, b] = point_source(source_x, mesh)
    % initialize
    [A, b] = initialize_A_b(length(cells_centroid_x));

    is_implicit = Config.IS_SOURCE_IMPLICIT;

end
