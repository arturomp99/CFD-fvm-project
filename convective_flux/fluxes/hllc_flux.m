function F = hllc_flux(UL, UR, gamma)
    % HLLC_FLUX - Flujo numérico HLLC para ecuaciones de Euler 1D
    %
    %   Implementa el solver de Riemann HLLC (Harten-Lax-van Leer-Contact)
    %   para resolver el problema de Riemann local entre dos estados.
    %   Captura choques, rarefacciones y discontinuidades de contacto.
    %
    % Input
    %   UL     - Vector de estado conservativo izquierdo [ρ; ρu; E] [3×1]
    %   UR     - Vector de estado conservativo derecho [ρ; ρu; E] [3×1]
    %   gamma  - Razón de calores específicos γ = cp/cv [-]
    %
    % Ouput
    %   F      - Flujo numérico HLLC en la interfaz [3×1]
    %
    % ALGORITMO HLLC:
    %   Divide el dominio en 4 regiones separadas por 3 ondas:
    %   - Onda izquierda SL (choque o rarefacción)
    %   - Discontinuidad de contacto SM (velocidad constante)
    %   - Onda derecha SR (choque o rarefacción)
    %

    % PASO 1: CONVERSIÓN A VARIABLES PRIMITIVAS
    [rhoL, uL, pL, EL] = cons2prim(UL, gamma); % Estado izquierdo
    [rhoR, uR, pR, ER] = cons2prim(UR, gamma); % Estado derecho

    % PASO 2: FLUJOS FÍSICOS EN LOS ESTADOS ORIGINALES
    FL = euler_flux_1d(UL, gamma); % Flujo físico izquierdo
    FR = euler_flux_1d(UR, gamma); % Flujo físico derecho

    % PASO 3: VELOCIDADES DEL SONIDO Y ESTIMACIÓN DE ONDAS
    aL = sqrt(gamma * pL / rhoL); % Velocidad del sonido izquierda
    aR = sqrt(gamma * pR / rhoR); % Velocidad del sonido derecha

    % Estimación de velocidades de ondas extremas (Davis, 1988)
    % SL: velocidad onda izquierda, SR: velocidad onda derecha
    SL = min(uL - aL, uR - aR); % Velocidad mínima (onda izquierda)
    SR = max(uL + aL, uR + aR); % Velocidad máxima (onda derecha)

    % PASO 4: VELOCIDAD DE DISCONTINUIDAD DE CONTACTO
    % Cálculo de SM basado en conservación de masa y momentum
    % Fórmula de Toro (1999) para velocidad de contacto
    denom = rhoL * (SL - uL) - rhoR * (SR - uR);

    % Protección contra división por cero en casos degenerados
    if abs(denom) < 1e-12
        denom = 1e-12;
    end

    % SM: velocidad de la discontinuidad de contacto
    SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom;

    % =====================================================================
    % PASO 5: DETERMINACIÓN DE REGIÓN Y FLUJO CORRESPONDIENTE
    % =====================================================================
    % El dominio se divide en 4 regiones según la posición relativa de las ondas

    % REGIÓN 1: Completamente a la izquierda (0 ≤ SL)
    % Solo influye el estado izquierdo original
    if 0 <= SL
        F = FL; % Flujo físico izquierdo
        return;
    end

    % REGIÓN 2: Entre onda izquierda y contacto (SL ≤ 0 ≤ SM)
    % Estado intermedio izquierdo U*L influenciado por onda SL
    if SL <= 0 && 0 <= SM
        % Densidad en región estrella izquierda (conservación de masa)
        rhoStarL = rhoL * (SL - uL) / (SL - SM);

        % Energía total en región estrella (conservación de energía)
        EStarL = rhoStarL * (EL / rhoL + (SM - uL) * (SM + pL / (rhoL * (SL - uL))));

        % Estado conservativo intermedio izquierdo
        UStarL = [rhoStarL; rhoStarL * SM; EStarL];

        % Flujo HLLC usando relación de Rankine-Hugoniot
        F = FL + SL * (UStarL - UL);
        return;
    end

    % REGIÓN 3: Entre contacto y onda derecha (SM ≤ 0 ≤ SR)
    % Estado intermedio derecho U*R influenciado por onda SR
    if SM <= 0 && 0 <= SR
        % Densidad en región estrella derecha (conservación de masa)
        rhoStarR = rhoR * (SR - uR) / (SR - SM);

        % Energía total en región estrella (conservación de energía)
        EStarR = rhoStarR * (ER / rhoR + (SM - uR) * (SM + pR / (rhoR * (SR - uR))));

        % Estado conservativo intermedio derecho
        UStarR = [rhoStarR; rhoStarR * SM; EStarR];

        % Flujo HLLC usando relación de Rankine-Hugoniot
        F = FR + SR * (UStarR - UR);
        return;
    end

    % REGIÓN 4: Completamente a la derecha (SR ≤ 0)
    % Solo influye el estado derecho original
    F = FR; % Flujo físico derecho
end
