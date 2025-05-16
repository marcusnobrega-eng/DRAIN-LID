%% =========================================================================
% 📂 File Location : Physical_Functions_And_Utilities/K_vgm.m
%
% 💧 FUNCTION: K_vgm
% Purpose    : Computes unsaturated hydraulic conductivity K(h) [m/s]
%              using the Mualem–Van Genuchten formulation.
%
% Inputs:
%   h        – Pressure head [m]                         [vector or scalar]
%   Ks       – Saturated hydraulic conductivity [m/s]
%   theta_r  – Residual water content [m³/m³]
%   theta_s  – Saturated water content [m³/m³]
%   alpha    – Van Genuchten parameter [1/m]
%   n        – Van Genuchten parameter [-]
%   m        – Typically m = 1 - 1/n
%
% Output:
%   K        – Unsaturated hydraulic conductivity [m/s]
%
% Notes:
%   • Based on Van Genuchten–Mualem model:
%       K(h) = Ks · Se^0.5 · [1 - (1 - Se^{1/m})^m]^2
%   • Se = effective saturation = (θ - θ_r) / (θ_s - θ_r)
%
% Author     : Marcus Nóbrega, Ph.D.
% Updated    : May 2025
%% =========================================================================
function K = K_vgm(h, Ks, theta_r, theta_s, alpha, n, m)

    % === Step 1: Compute volumetric water content θ(h)
    theta = theta_vgm(h, theta_r, theta_s, alpha, n, m);  % [m³/m³]

    % === Step 2: Compute effective saturation Se = (θ - θ_r) / (θ_s - θ_r)
    Se = (theta - theta_r) ./ (theta_s - theta_r);  % [-]

    % === Step 3: Apply Mualem–Van Genuchten model
    % Handle dry-end precision to avoid imaginary/NaN values
    Se = max(Se, 1e-8);         % Avoid 0^negative
    Se = min(Se, 1.0);          % Cap at saturation

    term = (1 - Se.^(1./m)).^m; % Core nonlinearity
    K = Ks .* Se.^0.5 .* (1 - term).^2;  % [m/s]

end
