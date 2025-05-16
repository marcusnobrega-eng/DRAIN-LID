%% =========================================================================
% 📂 File Location: Physical_Functions_And_Utilities/get_boundary_values.m
%
% 🔁 BOUNDARY CONDITION EVALUATION FUNCTION
% Description : Computes top and bottom boundary condition values (head or flux)
%               and handles ponding dynamics based on infiltration/exfiltration
%               capacity at the surface.
%
% Inputs:
%   h                - Current pressure head profile            [Nz x 1]
%   params           - Model parameter struct (must include soil & BC config)
%   t                - Current simulation time                  [scalar, s]
%   ponding_depth_prev - Ponded water depth from previous step [scalar, m]
%
% Outputs:
%   top_bc        - Top boundary condition value (flux or head) [m/s or m]
%   bottom_bc     - Bottom boundary condition value              [m/s or m]
%   ponding_depth - Updated surface ponded water depth           [m]
%
% Author  : Marcus Nóbrega, Ph.D.
% Updated : May 2025
%
% z = 0         ← Top of domain (soil surface)
% │
% │    ┌─────── q(Nz+1) ───────→  (Top boundary flux: rainfall, ET)
% │    │
% │    │   h(Nz) ← Top node          (Dirichlet or Neumann BC)
% │    ├─────── q(Nz)
% │    │
% │    │   h(Nz-1)
% │    ├─────── q(Nz-1)
% │    │
% │    │     ...
% │    │
% │    │   h(2)
% │    ├─────── q(2)
% │    │
% │    │   h(1) ← Bottom node       (Dirichlet, Neumann, Free, or NoFlow)
% │    └─────── q(1)
% │
%
% z = -L        ← Bottom of domain

function [top_bc, bottom_bc, ponding_depth] = get_boundary_values(h, params, t, ponding_depth_prev)

%% === ⚙️ TOP BOUNDARY CONDITION ======================================

if params.top_bc_type == "neumann"
    % ⛅ Time-varying surface flux (e.g., rainfall or ET)
    surface_flux = interp1(params.surface_flux_time, params.surface_flux_vals, ...
        t, 'linear', 'extrap');  % [m/s]

    % Extract surface conductivity and gradient
    K         = K_vgm(h, params.Ks, params.theta_r, params.theta_s, ...
        params.alpha, params.n, params.m);
    dhdz_top  = ((max(h(end), ponding_depth_prev)) - h(end-1)) / params.dz(end);
    K_face    = 0.5 * (K(end) + K(end-1));  % Upwinded face conductivity

    % Compute downward infiltration capacity (Richards flux)
    q_capacity = -K_face * (dhdz_top + 1);  % [m/s]

    if surface_flux < 0
        % 🌧️ Rainfall case (downward flux)
        if surface_flux < q_capacity
            % Rainfall exceeds capacity → ponding
            excess = q_capacity - surface_flux;
            ponding_depth = ponding_depth_prev + excess * params.dt;
            top_bc = q_capacity;
        else
            % Rainfall within capacity
            top_bc = surface_flux;
            ponding_depth = max(0, ponding_depth_prev + (q_capacity - surface_flux) * params.dt);
        end
    
    elseif surface_flux > 0
        alpha = feddes_alpha(h(end), params.h_lim_upper, params.h_lim_down); % Alpha factor in evaporation
        q_ETP = alpha * surface_flux;
        % 💨 Evaporation case (upward flux)
        if q_ETP > abs(q_capacity)
            % Evaporation exceeds capillary supply
            evap_cap = abs(q_capacity);
            top_bc = evap_cap;
            ponding_depth = max(0, ponding_depth_prev - (q_ETP - evap_cap) * params.dt);
        else
            top_bc = q_ETP;
            ponding_depth = max(0, ponding_depth_prev);  % No change
        end

    else
        % 🛑 Zero surface flux
        top_bc = 0;
        ponding_depth = ponding_depth_prev;
    end

else
    % 📏 Dirichlet case: fixed surface pressure head
    top_bc = params.top_bc_value;
    ponding_depth = max(top_bc, 0);  % Optional: ponding = max(0, h_surface)
end

%% === ⚙️ BOTTOM BOUNDARY CONDITION ===================================

if params.bottom_bc_type == "neumann"
    bottom_bc = interp1(params.bottom_flux_time, params.bottom_flux_vals, ...
        t, 'linear', 'extrap');  % [m/s]

elseif params.bottom_bc_type == "dirichlet"
    bottom_bc = params.bottom_bc_value;

elseif params.bottom_bc_type == "free"
    % Unit gradient condition (∂h/∂z + 1 = 0)
    % Implemented inside residual as q = -K
    bottom_bc = NaN;

elseif params.bottom_bc_type == "noflow"
    % Conditional seepage face (handled in residual)
    bottom_bc = NaN;

else
    error('❌ Unknown bottom boundary type: %s', params.bottom_bc_type);
end
end
