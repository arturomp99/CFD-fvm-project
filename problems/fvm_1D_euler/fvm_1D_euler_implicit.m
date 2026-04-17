function [A, b] = fvm_1D_euler_implicit(state, cells)
    %FVM_1D_EULER_IMPLICIT
    %   Linearised implicit form of the 1D Euler residual:
    %       d(state)/dt = A*state + b

    residual_fun = @(w) residual_only(w, cells);

    [A, b] = linearize_residual_fd(residual_fun, state);
end

function R = residual_only(state, cells)
    [~, R] = rusanov_interpolator(state, cells);
end
