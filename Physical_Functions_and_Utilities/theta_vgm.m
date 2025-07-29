%% =========================================================================
% 📂 File Location : Physical_Functions_And_Utilities/theta_vgm.m
%
% 🌱 FUNCTION: theta_vgm
% Purpose    : Computes volumetric water content θ(h) using the Van Genuchten
%              soil retention model for unsaturated media.
%
% Inputs:
%   h        – Pressure head [m]                      [vector or scalar]
%   theta_r  – Residual water content [m³/m³]
%   theta_s  – Saturated water content [m³/m³]
%   alpha    – Van Genuchten parameter [1/m]
%   n        – Van Genuchten parameter [-]
%   m        – Usually defined as m = 1 - 1/n
%
% Outputs:
%   theta    – Volumetric water content [m³/m³]
%
% Notes:
%   • Model applies for h ≤ 0 (unsaturated zone), but is extended to h > 0
%     by setting θ(h > 0) = θ_s
%   • Vectorized for fast computation
%
% Author     : Marcus Nóbrega, Ph.D.
% Updated    : May 2025
%% =========================================================================
function [theta, Se] = theta_vgm(h, theta_r, theta_s, alpha, n, m)

    % === Compute effective saturation Se = [(1 + (α|h|)^n)]^(-m)
    abs_h = abs(h);  % Ensure function is valid for all h
    Se = 1 ./ ((1 + (alpha .* abs_h).^n).^m);  % [0,1] (non-clipped)

    % === Compute volumetric water content
    theta = theta_r + (theta_s - theta_r) .* Se;

    % === Enforce saturation for h > 0
    theta(h > 0) = theta_s(h > 0);

    % === Optional: Clip θ to physical bounds (optional enhancement)
    % theta = min(max(theta, theta_r), theta_s);  % <- Uncomment if needed

end
