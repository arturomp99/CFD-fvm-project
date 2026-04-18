function boundary_state = open_to_atmosphere(face_state, gamma)
    % OPEN_TO_ATMOSPHERE - Boundary condition for open to atmosphere
    %
    %   Sets boundary to atmospheric conditions: P = 101325 Pa, rho = 1.225 kg/m³, u = 0
    %   This represents a boundary open to standard atmospheric conditions.
    %
    % Input
    %   face_state  - State vector from interior [ρ; ρu; E] [3×1] (unused)
    %   gamma       - Specific heat ratio γ = cp/cv [-]
    %
    % Output
    %   boundary_state - Atmospheric state [ρ; ρu; E] [3×1]

    % Atmospheric conditions
    P_atm = Air.SEA_LEVEL_PRESSURE; % Atmospheric pressure [Pa]
    rho_atm = Air.SEA_LEVEL_DENSITY; % Atmospheric density [kg/m³]
    u_atm = 0.0; % Zero velocity at atmosphere

    % Compute total energy for atmospheric state
    internal_energy = P_atm / (rho_atm * (gamma - 1));
    E_atm = rho_atm * (internal_energy + 0.5 * u_atm ^ 2);

    % Return atmospheric state
    boundary_state = [rho_atm; rho_atm * u_atm; E_atm];
end
