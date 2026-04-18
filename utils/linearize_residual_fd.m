function [A, b] = linearize_residual_fd(residual_fun, w)
    %LINEARIZE_RESIDUAL_FD
    %   Linearise a nonlinear residual R(w) around the current state:
    %       R(w) ~= A*w + b
    %   using forward finite differences.

    n = length(w);
    R0 = residual_fun(w);

    A = sparse(n, n);

    for j = 1:n
        delta = 1e-7 * max(1.0, abs(w(j)));

        w_pert = w;
        w_pert(j) = w_pert(j) + delta;

        Rj = residual_fun(w_pert);

        A(:, j) = (Rj - R0) / delta;
    end

    b = R0 - A * w;
end
