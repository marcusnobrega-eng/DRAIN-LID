%% =========================================================================
% 📂 File Location: Physical_Functions_And_Utilities/dtheta_vgm.m
%
% 🔁 FUNCTION: dtheta_vgm
% Purpose    : Computes the derivative dθ/dh of the volumetric water content
%              with respect to pressure head using the Van Genuchten model.
%
% --------------------------------------------------------------------------
% FORMULATION:
%
%   θ(h) = θ_r + (θ_s - θ_r) · S_e(h)           if h < 0
%   θ(h) = θ_s                                  if h ≥ 0
%
%   S_e(h) = [1 + (α·|h|)^n]^(-m)
%   dθ/dh = - (θ_s - θ_r) · α·n·m·(α·|h|)^(n-1) / [1 + (α·|h|)^n]^(m+1)
%
% --------------------------------------------------------------------------
% Inputs:
%   h        : Pressure head [m] (vector or scalar, negative for unsaturated)
%   theta_r  : Residual water content [–]
%   theta_s  : Saturated water content [–]
%   alpha    : Van Genuchten parameter [1/m]
%   n        : Van Genuchten parameter [–]
%   m        : Van Genuchten parameter [–] (usually m = 1 - 1/n)
%
% Output:
%   dtheta   : Derivative of θ with respect to h [1/m]
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : June 2025
%% =========================================================================

function dtheta = dtheta_vgm(h, theta_r, theta_s, alpha, n, m)

    % Initialize output with zeros
    dtheta = zeros(size(h));

    % Identify unsaturated domain (h < 0)
    idx = h < 0;
    h_unsat = h(idx);

    % Compute intermediate terms
    abs_ah = abs(alpha .* h_unsat);                             % |α·h|
    term   = (1 + abs_ah.^n);                                  % Denominator
    denom  = term.^(m + 1);                                    % [1 + (α·|h|)^n]^(m+1)
    num    = alpha .* n .* m .* abs_ah.^(n - 1);                  % Numerator

    % Compute dSe/dh
    dSe_dh = num ./ denom;                                     % [1/m]

    % Final expression for dθ/dh
    dtheta(idx) = -(theta_s - theta_r) .* dSe_dh;

end
