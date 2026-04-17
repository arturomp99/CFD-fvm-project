function F = hllc_flux(UL, UR, gamma)
    %HLLC_FLUX Compute HLLC numerical flux for 1D Euler equations.

    [rhoL, uL, pL, EL] = cons2prim(UL, gamma);
    [rhoR, uR, pR, ER] = cons2prim(UR, gamma);

    FL = euler_flux_1d(UL, gamma);
    FR = euler_flux_1d(UR, gamma);

    aL = sqrt(gamma * pL / rhoL);
    aR = sqrt(gamma * pR / rhoR);

    SL = min(uL - aL, uR - aR);
    SR = max(uL + aL, uR + aR);

    denom = rhoL * (SL - uL) - rhoR * (SR - uR);

    if abs(denom) < 1e-12
        denom = 1e-12;
    end

    SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom;

    if 0 <= SL
        F = FL;
        return;
    end

    if SL <= 0 && 0 <= SM
        rhoStarL = rhoL * (SL - uL) / (SL - SM);
        EStarL = rhoStarL * (EL / rhoL + (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
        UStarL = [rhoStarL; rhoStarL * SM; EStarL];
        F = FL + SL * (UStarL - UL);
        return;
    end

    if SM <= 0 && 0 <= SR
        rhoStarR = rhoR * (SR - uR) / (SR - SM);
        EStarR = rhoStarR * (ER / rhoR + (SM - uR) * (SM + pR / (rhoR * (SR - uR))));
        UStarR = [rhoStarR; rhoStarR * SM; EStarR];
        F = FR + SR * (UStarR - UR);
        return;
    end

    F = FR;
end
