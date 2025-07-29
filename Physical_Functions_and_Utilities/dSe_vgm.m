function dSe = dSe_vgm(h, alpha, n, m)
% =========================================================================
% 🌱 FUNCTION: dSe_vgm
% Purpose   : Computes the derivative of the effective saturation (Se)
%             with respect to pressure head (h), based on the van Genuchten model.
%
% Equation  :
%   For h < 0 (unsaturated zone):
%     Se(h)   = [1 + (α·|h|)^n]^(-m)
%     dSe/dh  = -m·n·α^n·|h|^(n-1) · [1 + (α·|h|)^n]^(-m-1) · sign(h)
%
%   For h ≥ 0 (saturated zone): dSe/dh = 0
%
% Inputs:
%   h      – Pressure head [m] (can be scalar or vector)
%   alpha  – van Genuchten parameter, inverse air-entry suction [1/m]
%   n      – van Genuchten parameter (typically > 1)
%   m      – Shape parameter, m = 1 - 1/n
%
% Output:
%   dSe    – Derivative of effective saturation w.r.t. h [1/m]
%
% Author  : Marcus Nóbrega (2025)
% =========================================================================

    % Initialize output with zeros (saturated zone will remain 0)
    dSe = zeros(size(h));

    % Logical mask for unsaturated zone (where h < 0)
    mask = h < 0;

    % Compute absolute value of pressure head for unsaturated cells
    abs_h = abs(h(mask));

    % Compute Se(h) in unsaturated zone
    % Se = [1 + (α·|h|)^n]^(-m)
    Se = (1 + (alpha .* abs_h).^n).^(-m);

    % Compute derivative dSe/dh for h < 0:
    % dSe/dh = -m·n·α^n·|h|^(n-1) * [1 + (α·|h|)^n]^(-m-1)
    dSe(mask) = -m .* n .* alpha.^n .* abs_h.^(n - 1) ...
                .* (1 + (alpha .* abs_h).^n).^(-m - 1);

    % Correct the sign using sign(h) to maintain generality
    dSe(mask) = dSe(mask) .* sign(h(mask));  % Optional: sign(h) = -1 here

end
