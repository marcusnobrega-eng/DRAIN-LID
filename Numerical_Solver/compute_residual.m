% =========================================================================
% 📂 File Location : Numerical_Solver/compute_residual.m
%
% 🔁 FUNCTION: compute_residual
% Purpose    : Computes the mixed-form Richards Equation residual vector F(h)
%              and Darcy fluxes at all vertical interfaces (including BCs).
%
% Inputs:
%   h           – Current pressure head profile               [Nz x 1]
%   h_old       – Previous pressure head profile              [Nz x 1]
%   theta       – Current water content                       [Nz x 1]
%   theta_old   – Previous time step water content            [Nz x 1]
%   K           – Hydraulic conductivity profile              [Nz x 1]
%   params      – Struct with dt, dz, S_s, Nz, BC flags
%   top_bc      – Top BC value: flux or head [m/s or m]
%   bottom_bc   – Bottom BC value: flux or head [m/s or m]
%   source_term – Source/sink term [m³/m³/s]                  [Nz x 1]
%
% Outputs:
%   F      – Residual vector for Newton iteration             [Nz x 1]
%   q      – Vertical fluxes at all interfaces                [1 x (Nz+1)]
%   q_bot  – Bottom boundary flux (q(1))                      [scalar]
%   q_top  – Top boundary flux    (q(Nz+1))                   [scalar]
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
%
% =========================================================================
% 💧 1D Vertical Soil Column Schematic (Nz nodes, Nz+1 fluxes)
%
%   Coordinate system:
%       z =  0       → Top boundary (soil surface)
%       z = -L       → Bottom boundary (deepest node)
%       z increases upward 📈 (from -L (i = 1) to 0 (i = Nz))
%
%   Flux sign convention (from Darcy’s law):
%       • q < 0 → downward flow (with gravity)
%       • q > 0 → upward flow   (against gravity)
%
%     z =  0  (Top)
%     │             ← Top boundary flux (e.g., Neumann or computed)
%     │    ┌───↓ q(Nz+1)──┐
%     │    │   Node h(Nz) │ ← Top node (surface)
%     │    ├───↓ q(Nz)────┤
%  .  │    │   Node h(Nz-1)
% z.  │    ├───↓ q(Nz-1)──┤
%  .  │    │     ...
%     │    ├───↓ q(3)─────┤
%     │    │   Node h(2)
%     │    ├───↓ q(2)─────┤
%     │    │   Node h(1)  │ ← Bottom node
%     │    ├───↓ q(1)─────┤
%     │
%     Bottom boundary flux (e.g., Dirichlet, Neumann, Free)
%     z = -L  (Bottom)
%
%% ----------------------------- THEORY ------------------------------------
%
% 📘 Mixed-Form Richards Equation (1D Vertical)
%
% The function computes the residual vector F(h) for Newton-Raphson by evaluating
% the mass balance of water in each soil control volume (cell) in a vertically
% discretized soil column. The formulation follows the **mixed-form** of Richards:
%
%   ∂θ(h)/∂t + S_s ∂h/∂t = ∂/∂z [ K(h) (∂h/∂z + 1) ] + S(z,t)
%
% where:
%   - θ(h): volumetric water content [m³/m³]
%   - h   : pressure head [m]
%   - K(h): unsaturated hydraulic conductivity [m/s]
%   - S_s : specific storage [1/m]
%   - S   : source/sink term [1/s]
%
% ---------------------------------------------
% 🔁 Discrete Control Volume Residual (node i):
%
%   F_i = [θᵢⁿ⁺¹ - θᵢⁿ]/Δt + S_sᵢ [hᵢⁿ⁺¹ - hᵢⁿ]/Δt
%         + (qᵢ₊₁ - qᵢ)/Δz - Sᵢ
%
% where:
%   - qᵢ = -K_face * (∂h/∂z + 1)  → Darcy flux [m/s]
%   - qᵢ₊₁ and qᵢ are fluxes at the interfaces above and below node i
%   - Residual is defined so F_i = 0 at steady state or convergence
%
% ---------------------------------------------
% ⚠️ Boundary Nodes:
%
% - Node 1 (Bottom): applies "dirichlet", "neumann", "free", or "noflow"
% - Node Nz (Top): same logic, adjusted for surface boundary processes
%
% - Flux array q(i) is defined at interface *above* node i.
%   → Total of Nz+1 fluxes for Nz nodes.
%
% --------------------------------------------------------------------------
%% =========================================================================
function [F, q, q_bot, q_top] = compute_residual(h, h_old, theta, theta_old, K, ...
                                                 params, top_bc, bottom_bc, source_term, top_bc_type_used)

    %% === 1. Initialize Arrays ============================================
    Nz = params.Nz;
    dz = params.dz;
    dt = params.dt;

    F = zeros(1, Nz);
    q = zeros(1, Nz + 1);  % q(i) = flux *above* node i

    %% === 2. Apply Bottom Boundary Condition (Node 1) =====================
    switch params.bottom_bc_type
        case "dirichlet"
            % Enforced head condition
            F(1) = bottom_bc - h(1);  % Enforce h(1)
            dhdz_bot = (h(2) - h(1)) / dz(1);
            K_bot = 0.5 * (K(1) + K(2));
            q(1) = -K_bot * (dhdz_bot + 1);  % Used for mass balance

        case "neumann"
            % Enforced flux condition
            q(1) = bottom_bc;
            K_bot = 0.5 * (K(1) + K(2));
            dhdz_bot = (h(2) - h(1)) / dz(1);
            F(1) = bottom_bc + K_bot * (dhdz_bot + 1);

        % case "neumann"
        %     % Apply enforced Neumann flux at bottom interface
        %     q(1) = bottom_bc;
        % 
        %     % Compute internal flux just above bottom node
        %     K_face = 0.5 * (K(1) + K(2));
        %     dhdz   = (h(2) - h(1)) / dz(1);
        %     q(2)   = -K_face * (dhdz + 1);
        % 
        %     % Full residual at bottom node (mass balance)
        %     dqdz = -(q(2) - q(1)) / dz(1);
        % 
        %     F(1) = (theta(1) - theta_old(1)) / dt ...
        %          + params.S_s(1) * (h(1) - h_old(1)) / dt ...
        %          - dqdz ...
        %          - source_term(1);

        case "free"
            % Unit gradient condition: ∂h/∂z = -1 → q = -K
            F(1) = h(2) - h(1);
            q(1) = -K(1);

        case "noflow"
            % No-flow: q(1) = 0 → residual from internal gradient
            q(1) = 0;
            F(1) = (h(2) - h(1)) / dz(1) + 1;

        otherwise
            error("Unknown bottom_bc_type: %s", params.bottom_bc_type);
    end

    %% === 3. Apply Top Boundary Condition (Node Nz) =======================
        % switch params.top_bc_type
        %     case "dirichlet"
        %         % Enforced head condition
        %         F(Nz) = top_bc - h(Nz);
        %         dhdz_top = (h(Nz) - h(Nz-1)) / dz(Nz);
        %         K_top = 0.5 * (K(Nz) + K(Nz-1));
        %         q(Nz+1) = -K_top * (dhdz_top + 1);
        % 
        %     case "neumann"
        %         % Enforced flux
        %         q(Nz+1) = top_bc;
        %         K_top = 0.5 * (K(Nz) + K(Nz-1));
        %         dhdz_top = (h(Nz) - h(Nz-1)) / dz(Nz);
        %         F(Nz) = top_bc + K_top * (dhdz_top + 1);

        % case "neumann"
        %     % Apply enforced Neumann flux directly
        %     q(Nz+1) = top_bc;
        % 
        %     % Compute flux just below top node
        %     K_face = 0.5 * (K(Nz-1) + K(Nz));
        %     dhdz   = (h(Nz) - h(Nz-1)) / dz(Nz);
        %     q(Nz)  = -K_face * (dhdz + 1);
        % 
        %     % Full residual at top node (like any interior node)
        %     dqdz = -(q(Nz+1) - q(Nz)) / dz(Nz);
        % 
        %     F(Nz) = (theta(Nz) - theta_old(Nz)) / dt ...
        %           + params.S_s(Nz) * (h(Nz) - h_old(Nz)) / dt ...
        %           - dqdz ...
        %           - source_term(Nz);
        % 
        % 
        % 
        %     otherwise
        %         error("Unknown top_bc_type: %s", params.top_bc_type);
        % end

    switch top_bc_type_used
        case "dirichlet"
            % Enforce head directly
            F(Nz) = top_bc - h(Nz);
            dhdz_top = (h(Nz) - h(Nz-1)) / dz(Nz);
            K_top = 0.5 * (K(Nz) + K(Nz-1));
            q(Nz+1) = -K_top * (dhdz_top + 1);

        case "neumann"
            % Flux-specified boundary
            q(Nz+1) = top_bc;
            K_top = 0.5 * (K(Nz) + K(Nz-1));
            dhdz_top = (h(Nz) - h(Nz-1)) / dz(Nz);
            F(Nz) = top_bc + K_top * (dhdz_top + 1);

        otherwise
            error("Unknown top_bc_type_used: %s", top_bc_type_used);
    end

    %% === 4. Compute Internal Face Gradients & Fluxes =====================
    dh    = diff(h);                        % [Nz-1 x 1]
    dhdz  = dh ./ dz(1:end-1);              % Vertical gradients
    K_face = 0.5 * (K(1:end-1) + K(2:end)); % Interface K
    q(2:Nz) = -K_face .* (dhdz + 1);        % Darcy flux at internal faces

    %% === 5. Compute Flux Divergence and Residuals ========================
    dqdz = -diff(q) ./ dz;  % ∂q/∂z from q(i+1) - q(i)

    % Mass balance at internal nodes
    int_nodes = 2:(Nz-1);
    F(int_nodes) = (theta(int_nodes) - theta_old(int_nodes)) / dt ...
        + params.S_s(int_nodes) .* (h(int_nodes) - h_old(int_nodes)) / dt ...
        - dqdz(int_nodes) ...
        - source_term(int_nodes);

    %% === 6. Output Boundary Fluxes =======================================
    q_bot = q(1);
    q_top = q(end);
end 