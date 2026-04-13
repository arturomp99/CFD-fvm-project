function dt = euler_cfl_dt(w, t, cells, CFL)
 
    num_cells = length(cells);

    rho  = w(1:num_cells);
    rhou = w(num_cells + 1 : 2 * num_cells);
    E    = w(2 * num_cells + 1 : 3 * num_cells);

    rho = max(rho, 1e-10);

    u = rhou ./ rho;
    p = (Air.GAMMA - 1) .* (E - 0.5 .* (rhou.^2) ./ rho);
    p = max(p, 1e-10);

    c = sqrt(Air.GAMMA .* p ./ rho);
    lambda = abs(u) + c;

    x = zeros(num_cells, 1);
    for i = 1:num_cells
        x(i) = cells(i).centroid(1);
    end

    % Muy importante: eliminar x repetidas
    x_unique = unique(sort(x));

    if numel(x_unique) < 2
        error('Not enough distinct x-coordinates to compute CFL timestep.');
    end

    dx_min = min(diff(x_unique));
    dx_min = max(dx_min, 1e-6);

    max_lambda = max(lambda);
    max_lambda = max(max_lambda, 1e-10);

    dt = CFL * dx_min / max_lambda;
end