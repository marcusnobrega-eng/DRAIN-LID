function dK = dK_vgm(h, Ks, theta_r, theta_s, alpha, n, m)
% =========================================================================
% ⛏️ FUNCTION: dK_vgm
% Purpose   : Compute dK/dh (derivative of hydraulic conductivity w.r.t. h)
%             using the Van Genuchten–Mualem model.
%
% Equation  : 
%   K(h) = Ks * Se^{1/2} * [1 - (1 - Se^{1/m})^m]^2
%   dK/dh = dK/dSe * dSe/dh
%
% Inputs:
%   h        – Pressure head [m] (negative in unsaturated zone)
%   Ks       – Saturated hydraulic conductivity [m/s]
%   theta_r  – Residual water content [–]
%   theta_s  – Saturated water content [–]
%   alpha    – Inverse of air-entry pressure [1/m]
%   n, m     – van Genuchten shape parameters (m = 1 - 1/n)
%
% Output:
%   dK       – Derivative dK/dh [m/s per m]
%
% Author    : Marcus Nóbrega (2025)
% =========================================================================

    % === Effective saturation and its derivative ==========================
    [~, Se] = theta_vgm(h, theta_r, theta_s, alpha, n, m);  % Effective saturation
    dSe = dSe_vgm(h, alpha, n, m) ./ (theta_s - theta_r);  % dSe/dh

    % === Component terms for dK/dh ========================================
    term1 = 0.5 * Se.^(-0.5) .* dSe;  % dA/dh = d/dh of Se^{1/2}
    
    term2 = (1 - (1 - Se.^(1/m)).^m).^2;  % B term
    
    term3 = Se.^(0.5);  % A term
    
    % dB/dh, expanded via chain rule
    dB_dSe = 2 .* (1 - (1 - Se.^(1/m)).^m) ...
               .* (1 - Se.^(1/m)).^(m - 1) ...
               .* Se.^((1/m) - 1) ...
               .* (-1/m);  % dB/dSe
    
    term4 = dB_dSe .* dSe;  % dB/dh = dB/dSe * dSe/dh
    
    % === Final expression for dK/dh =======================================
    dK = Ks * (term1 .* term2 + term3 .* term4);
    
    % Ensure zero derivative in saturated zone (h ≥ 0)
    dK(h >= 0) = 0;
end
