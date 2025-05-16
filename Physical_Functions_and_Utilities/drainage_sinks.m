% =========================================================================
% 📂 File Location : Numerical_Solver/drainage_sinks.m
%
% 🔁 FUNCTION: drainage_sinks
% Purpose    : Computes orifice and spillway drainage flows and their
%              analytical derivatives for use in the 1D Richards Equation.
%              Includes volume-limited caps to ensure physical realism.
%
% ----------------------------- THEORY -------------------------------------
%
% Drainage is modeled using:
%
%   - Orifice Flow    : q = Ko · max(h - h0, 0)^e
%   - Spillway Flow   : q = c  · max(h - h_spill, 0)^e
%
% Each flow is capped based on the maximum physically allowed volume loss:
%
%   - Orifice Limit   : q ≤ h·(θ_s - θ_r) / Δt, h >= 0
%   - Spillway Limit  : q ≤ h / Δt, h >= 0
%
% The function returns:
%   - Total drainage (negative sink)
%   - Analytical derivatives d(q)/dh including cap effects
%
% --------------------------------------------------------------------------
%
% Inputs:
%   h               – Pressure head vector [Nz x 1]
%   theta_diff      – θ_s - θ_r vector [Nz x 1]
%   dz              – Spatial step vector [Nz x 1]
%   dt              – Time step scalar
%   K_orifice       – Orifice coefficient scalar
%   exp_orifice     – Orifice exponent scalar
%   spillway_enabled– Logical flag to enable spillway
%   c_spillway      – Spillway coefficient scalar
%   h_spill         – Spillway threshold scalar
%   exp_spillway    – Spillway exponent scalar
%
% Outputs:
%   q_total         – Total drainage flow vector (sink) [Nz x 1]
%   dqdh_orifice      – Analytical derivative dq/dh for orifice [Nz x 1]
%   dqdh_spill        – Analytical derivative dq/dh for spillway [Nz x 1]
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
% =========================================================================
function [q_total, dqdh_orifice, dqdh_spill, q_orifice, q_spill] = drainage_sinks( ...
    h, theta_diff, dz, dt, ...
    K_orifice, exp_orifice, ...
    spillway_enabled, c_spillway, ...
    h_spill, exp_spillway)

    %% === 1. Initialization ===============================================
    Nz = length(h);
    q_spill    = zeros(1, Nz);
    dqdh_spill   = zeros(1, Nz);

    %% === 2. Orifice Drainage =============================================
    % Raw flow: q = K_o * max(h, 0)^e
    min_depth = 1e-4;               % Minimum depth for derivatives
    h_eff_orifice = max(h - min_depth, 0);
    q_orifice = K_orifice .* (h_eff_orifice).^exp_orifice; % [m^3/s / m^2]
    h_eff_orifice_linearized = max(h, 0);
    % Derivative of raw flow
    active_raw = h_eff_orifice > 0;
    dqdh_orifice = ...
        K_orifice .* exp_orifice .* (h_eff_orifice_linearized).^(exp_orifice - 1);
    dqdh_orifice(~active_raw) = 0;
    % Maximum allowable flux
    dqdh_max = max(h_eff_orifice_linearized .* theta_diff / dt);
    % Final derivative
    dqdh_orifice = min(dqdh_orifice, dqdh_max);
    % Final derivative (min rule)
    % is_capped = q_raw_orifice >= q_cap_orifice;
    % dqdh_orifice(~is_capped) = dq_raw_orifice(~is_capped);
    % dqdh_orifice(is_capped)  = dq_cap_orifice(is_capped);
    % dqdh_orifice(~active_raw) = 0;

    %% === 3. Spillway Drainage (Optional) ================================
    if spillway_enabled
        % Raw flow: q = c * max(h - h_spill, 0)^e
        min_depth = 1e-4;               % Minimum depth for derivatives
        h_eff_spill = max(h - h_spill - min_depth, 0); % [m]
        active_spill = h_eff_spill > 0;

        h_eff_spill_linearized = max(h, 0); 

        q_spill = c_spillway * (h_eff_spill).^exp_spillway; % [m3/s/m2] (with original h_eff)
        q_spill(~active_spill) = 0;
        dqdh_spill = c_spillway * exp_spillway * (h_eff_spill_linearized).^(exp_spillway - 1);  % with linearized depth      
        dqdh_spill(~active_spill) = 0;
    end

    %% === 4. Combine Total Drainage ======================================
    q_total     = -(q_orifice + q_spill);   % Negative sign: drainage = sink
    dqdh_orifice  = -dqdh_orifice;
    dqdh_spill    = -dqdh_spill;
end
