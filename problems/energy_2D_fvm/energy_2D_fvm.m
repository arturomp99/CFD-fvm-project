function [A, b] = energy_2D_fvm(w, t, mesh, ...
    cp, cv, rho, k, speed, convective_method, conductive_method, ...
    boundary_condition_functions);
%ENERGY_2D_FVM Summary of this function goes here
%   Detailed explanation goes here

    [A_conv, b_conv] = convective_method(w, t, mesh, cp, cv, rho, ...
        speed);
    [A_cond, b_cond] = conductive_method(w, t, mesh, cv, rho, k);
    A_BC = zeros(length(mesh.volume));
    b_BC = zeros(length(mesh.volume), 1);
    for BC=boundary_condition_functions
        [A_, b_] = BC(w, t, mesh, cp, cv, rho, k, speed);
        A_BC = A_BC + A_;
        b_BC = b_BC + b_;
    end
    A = A_conv + A_cond + A_BC;
    b = b_conv + b_cond + b_BC;
end
