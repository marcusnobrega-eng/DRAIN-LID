%% =========================================================================
% 📂 File Location : Physical_Functions_And_Utilities/mass_balance_check.m
%
% 🔍 FUNCTION: mass_balance_check
% Purpose    : Evaluates mass conservation by comparing change in water storage
%              with net fluxes and source terms.
%
% ----------------------------- THEORY -------------------------------------
%
% For a given time step Δt, mass balance is:
%
%     ΔStorage = Net_Input
%
% Where:
%   ΔStorage     = [S(t+Δt) + Ponding(t+Δt)] - [S(t) + Ponding(t)]
%   Net_Input    = ∫(q_in - q_out) dt + ∫S(z,t) dz dt
%
% Discretized as:
%
%   ΔS = sum(θᵢⁿ⁺¹·Δzᵢ) - sum(θᵢⁿ·Δzᵢ) + ΔPonding
%
%   Net_Input = [-q_top_avg + q_bot_avg]·Δt + avg(Source)·Δt
%
% Mass balance error:
%   E_mass = ΔStorage - Net_Input
%
% --------------------------------------------------------------------------
%
% Inputs:
%   h_old, h_new          – Pressure head at t and t+Δt       [Nz x 1]
%   ponding_prev, ponding_depth – Surface ponding at t, t+Δt [scalar]
%   q_prev, q_now         – Interface fluxes at t, t+Δt       [Nz+1 x 1]
%   t, dt                 – Current time and time step        [s]
%   mb_error_cumulative   – Accumulated error                 [m]
%   params                – Model parameter struct
%   cumulative_net_flux_prev – Previously accumulated net inflow [m]
%
% Outputs:
%   mb_error              – Mass balance error (ΔS - Net Input) [m]
%   mb_error_cumulative   – Updated cumulative mass error       [m]
%   cumulative_net_flux   – Updated cumulative inflow            [m]
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
%% =========================================================================

function [mb_error, mb_error_cumulative, cumulative_net_flux] = ...
    mass_balance_check(h_old, h_new, ponding_prev, ponding_depth, ...
                       q_prev, q_now, t, dt, ...
                       mb_error_cumulative, params, cumulative_net_flux_prev, ...
                       Q_orifice_now, Q_spillway_now, ...
                       Q_orifice_prev, Q_spillway_prev, print_error)

    %% === 1. Compute θ at t and t+Δt =====================================
    theta_prev = theta_vgm(h_old, params.theta_r, params.theta_s, ...
                           params.alpha, params.n, params.m);
    theta_now  = theta_vgm(h_new, params.theta_r, params.theta_s, ...
                           params.alpha, params.n, params.m);

    %% === 2. Compute Soil Water Storage ==================================
    S_prev = sum(theta_prev(:) .* params.dz(:));  % [m]
    S_now  = sum(theta_now(:)  .* params.dz(:));  % [m]

    % Include ponding
    S_prev_total = S_prev + ponding_prev;
    S_now_total  = S_now  + ponding_depth;

    %% === 3. Compute Boundary Fluxes [m/s] ===============================
    q_top_now  = q_now(end);   q_bot_now  = q_now(1);
    q_top_prev = q_prev(end);  q_bot_prev = q_prev(1);

    q_top_avg = 0.5 * (q_top_now + q_top_prev);
    q_bot_avg = 0.5 * (q_bot_now + q_bot_prev);

    %% === 4. Compute Average Volumetric Source [m³/m³/s] =================
    source_now = interp1(params.source_times, params.source_profile', t, 'linear', 'extrap')';
    source_old = interp1(params.source_times, params.source_profile', t - dt, 'linear', 'extrap')';
    source_avg = 0.5 * (source_now + source_old);

    total_source = sum(source_avg .* params.dz(:)) * dt;  % [m]

    %% === 5. Include Drainage Sinks (Orifices and Spillway) ==============
    % Average Q_orifice and Q_spillway
    Q_orifice_avg  = 0.5 * (Q_orifice_now(:) + Q_orifice_prev(:));  % [m/s]
    Q_spillway_avg = 0.5 * (Q_spillway_now + Q_spillway_prev);      % [m/s]

    % Convert to total volume loss over time step
    orifice_volume  = sum(Q_orifice_avg .* dt);       % [m³/m²]
    spillway_volume = sum(Q_spillway_avg * dt);            % [m³/m²]
    total_sinks     = orifice_volume + spillway_volume;

    %% === 6. Compute Net Inflow [m] ======================================
    net_input = (-q_top_avg + q_bot_avg) * dt + total_source - total_sinks;
    cumulative_net_flux = cumulative_net_flux_prev + net_input;

    %% === 7. Compute Storage Change [m] ==================================
    delta_storage = S_now_total - S_prev_total;

    %% === 8. Compute Mass Balance Error [m] ==============================
    mb_error = delta_storage - net_input;
    mb_error_cumulative = mb_error_cumulative + mb_error;

    %% === 9. Report Mass Error (Optional Logging) ========================
    if abs(mb_error) > 1e-6 && print_error == 1
        fprintf('[t = %.2f min] ⚠️ Mass Balance Error: %.2e m³/m²\n', t / 60, mb_error);
    end
end
