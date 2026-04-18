function out = richardson_gci_three_grids(x3, phi3, x2, phi2, x1, phi1)
    %RICHARDSON_GCI_THREE_GRIDS
    % Grid 3 = coarse, grid 2 = medium, grid 1 = fine.
    % All comparisons are performed on the coarse grid x3.

    % Interpolate medium and fine solutions onto the coarse grid
    phi2_on_3 = interp1(x2, phi2, x3, 'linear');
    phi1_on_3 = interp1(x1, phi1, x3, 'linear');

    % Characteristic mesh sizes
    h3 = mean(diff(x3));
    h2 = mean(diff(x2));
    h1 = mean(diff(x1));

    r32 = h3 / h2;
    r21 = h2 / h1;

    if abs(r32 - r21) > 1e-6
        warning('Refinement ratios are not exactly equal: r32=%.6f, r21=%.6f', r32, r21);
    end

    r = 0.5 * (r32 + r21);

    % L1 norm errors
    E32 = mean(abs(phi3 - phi2_on_3));
    E21 = mean(abs(phi2_on_3 - phi1_on_3));

    % Observed order
    p = log(E32 / E21) / log(r);

    % Richardson extrapolation on coarse grid
    phi_ext = phi1_on_3 + (phi1_on_3 - phi2_on_3) / (r ^ p - 1);

    % Relative errors
    eps32 = E32 / max(mean(abs(phi2_on_3)), 1e-12);
    eps21 = E21 / max(mean(abs(phi1_on_3)), 1e-12);

    Fs = 1.25; % standard safety factor for 3+ grids

    GCI32 = 100 * Fs * eps32 / (r ^ p - 1);
    GCI21 = 100 * Fs * eps21 / (r ^ p - 1);

    % Asymptotic range check
    asymptotic_ratio = GCI32 / (r ^ p * GCI21);

    out.p = p;
    out.E32 = E32;
    out.E21 = E21;
    out.GCI32 = GCI32;
    out.GCI21 = GCI21;
    out.asymptotic_ratio = asymptotic_ratio;

    out.phi3 = phi3;
    out.phi2_on_3 = phi2_on_3;
    out.phi1_on_3 = phi1_on_3;
    out.phi_ext = phi_ext;
end
