function J = compute_jacobian_analytic(h, h_old, params, top_val, bottom_val, source_term)
% =========================================================================
% 📦 FUNCTION: compute_jacobian_analytic
%
% Purpose:
%   Computes the analytical Jacobian J = ∂F/∂h for the 1D mixed-form
%   Richards Equation using Newton's method and variable dz (non-uniform mesh).
%
% Inputs:
%   h           – Current pressure head [Nz x 1]
%   h_old       – Previous head [Nz x 1]
%   params      – Struct with fields:
%                   .Nz, .dz (vector of size Nz-1), .dt
%                   .Ks, .theta_r, .theta_s, .alpha, .n, .m, .A
%                   .top_bc_type, .bottom_bc_type
%   top_val     – Top BC value
%   bottom_val  – Bottom BC value
%   source_term – Volumetric sources/sinks [Nz x 1]
%
% Output:
%   J – Jacobian matrix [Nz x Nz]
% =========================================================================

    %% === 1. Initialization =============================================
    Nz = params.Nz;
    dz_vec = params.dz;       % Vector of spacing: size Nz-1
    dt = params.dt;

    % === Van Genuchten functions ===
    theta = theta_vgm(h, params.theta_r, params.theta_s, ...
                         params.alpha, params.n, params.m);
    C = dtheta_dpsi(h, params);
    K = K_vgm(h, params.Ks, params.theta_r, params.theta_s, ...
                     params.alpha, params.n, params.m);
    dK_dh = dK_dpsi_vgm(h, params);

    %% === 2. Build Gradient and Divergence for non-uniform dz ==========
    % Gradient: (h_i+1 - h_i) / dz_i
    G_data = [-1 ./ dz_vec, 1 ./ dz_vec];
    G = spdiags(G_data, [0 1], Nz - 1, Nz);

    % Divergence: −G^T (so mass conservation holds)
    D = -G';

    %% === 3. Harmonic Mean of K and Average of dK ======================
    % K_face = 2 .* K(1:end-1) .* K(2:end) ./ (K(1:end-1) + K(2:end));
    K_face = 1/2 * (K(1:end-1) + K(2:end));
    dK_face = 0.5 * (dK_dh(1:end-1) + dK_dh(2:end));

    K_diag = spdiags(K_face, 0, Nz-1, Nz-1);
    dK_diag = spdiags(dK_face, 0, Nz-1, Nz-1);

    h_grad = (G * h')';

    %% === 4. Assemble Jacobian Terms ===================================
    J_time  = (1 / dt) * spdiags(C, 0, Nz, Nz);                             % Node-based
    J_flux  = -D * K_diag * G;                                             % Linear flow term
    J_Kgrad = -D * spdiags(h_grad .* dK_face, 0, Nz-1, Nz-1) * G;          % Nonlinear term
    J_Kz    = -D * spdiags(dK_face, 0, Nz-1, Nz-1) * G;                    % Gravity correction
    
    J = J_time + J_flux + J_Kgrad + J_Kz;


    %% === 5. Apply Boundary Conditions ================================

    % === Bottom Node (i = 1) ===========================================
    switch lower(params.bottom_bc_type)
        case 'dirichlet'
            J(1,:) = 0;
            J(1,1) = 1;

        case 'neumann'
            % F(1) = (h(2) - h(1))/dz(1) - q_bc
            J(1,:) = 0;
            J(1,1) = -1 / dz_vec(1);
            J(1,2) =  1 / dz_vec(1);

        case 'free'
            % Gravity flow → no change needed
    end

    % === Top Node (i = Nz) =============================================
    switch lower(params.top_bc_type)
        case 'dirichlet'
            J(Nz,:) = 0;
            J(Nz,Nz) = 1;

        case 'neumann'
            % F(Nz) = (h(Nz) - h(Nz-1)) / dz(end) - q_bc
            J(Nz,:) = 0;
            J(Nz,Nz-1) = -1 / dz_vec(end);
            J(Nz,Nz)   =  1 / dz_vec(end);

        case 'free'
            % No-flux or gravity → handled naturally
    end
end


function C = dtheta_dpsi(psi, p)
% =========================================================================
% Computes the derivative dθ/dψ (specific moisture capacity)
% using the Van Genuchten model.
%
% Inputs:
%   psi  – Pressure head [Nz x 1]
%   p    – Struct with fields:
%          .theta_r, .theta_s, .alpha, .n
%
% Output:
%   C    – Specific moisture capacity [Nz x 1]
% =========================================================================

    % Ensure all parameter fields can be vectors
    alpha    = p.alpha;     % [Nz x 1] or scalar
    n        = p.n;         % [Nz x 1] or scalar
    theta_s  = p.theta_s;
    theta_r  = p.theta_r;

    % Compute exponent
    beta = n;   % commonly used in Van Genuchten form

    % Ensure all are column vectors
    psi      = psi(:);
    abs_psi  = abs(psi);
    sgn_psi  = sign(psi);

    % Avoid zero divisions (add epsilon)
    epsilon = 1e-12;
    denom = (alpha + abs_psi.^beta).^2 + epsilon;

    % Derivative
    C = -alpha .* (theta_s - theta_r) .* ...
         (beta .* abs_psi.^(beta - 1)) ./ denom .* sgn_psi;
end

function dK = dK_dpsi_vgm(psi, p)
% =========================================================================
% Computes the derivative dK/dψ for the Van Genuchten–Mualem model.
%
% Inputs:
%   psi – Pressure head [Nz x 1]
%   p   – Struct with:
%         .alpha      - [Nz x 1] or scalar, Van Genuchten α [1/cm]
%         .n          - [Nz x 1] or scalar, Van Genuchten n
%         .m          - [Nz x 1] or scalar, usually 1 - 1/n
%         .Ks         - Saturated hydraulic conductivity [Nz x 1]
%
% Output:
%   dK  – Derivative of hydraulic conductivity w.r.t. ψ [Nz x 1]
% =========================================================================

    % Unpack parameters
    alpha = p.alpha;
    n     = p.n;
    m     = p.m;   % typically m = 1 - 1/n
    Ks    = p.Ks;
    L     = 0.5;   % Common pore-connectivity parameter

    % Ensure vectors
    psi     = psi;
    abs_psi = abs(psi);
    sgn_psi = sign(psi);
    epsilon = 1e-12;

    % === Effective saturation Se ===
    Se = (1 + (alpha .* abs_psi).^n) .^ (-m);
    Se = min(max(Se, epsilon), 1 - epsilon);  % Clamp Se away from 0/1 for stability

    % === Derivative dSe/dψ ===
    dSe = -m .* n .* alpha.^n .* abs_psi.^(n - 1) ./ ...
          (1 + (alpha .* abs_psi).^n).^(m + 1) .* sgn_psi;

    % === Intermediate terms for dK/dψ ===
    % K(ψ) = Ks * Se^L * [1 - (1 - Se^{1/m})^m]^2

    A = Se.^L;                                    % Se^L
    B = (1 - (1 - Se.^(1./m)).^m);                % [1 - (1 - Se^{1/m})^m]
    dA_dSe = L .* Se.^(L - 1);                    % d(Se^L)/dSe
    dB_dSe = (1 - Se.^(1./m)).^(m - 1) .* ...
             Se.^(1./m - 1);                      % dB/dSe (chain rule)

    % Apply full product + chain rule
    dK = Ks .* ( ...
        dSe .* dA_dSe .* B.^2 + ...
        A .* 2 .* B .* m .* dB_dSe .* dSe ...
    );

    % === Saturated zone: ψ ≥ 0 → dK/dψ = 0 ===
    dK(psi >= 0) = 0;
end



