% =========================================================================
% 📂 File Location : Numerical_Solver/compute_jacobian.m
%
% 🔁 FUNCTION: compute_jacobian
% Purpose    : Computes the Jacobian matrix J = dF/dh by forward finite 
%              differences, used in Newton-Raphson updates for the
%              mixed-form Richards Equation.
%
% ----------------------------- THEORY -------------------------------------
%
% The Richards equation is nonlinear due to the dependence of θ(h) and K(h)
% on the pressure head h. To solve it using Newton-Raphson, we need:
%
%   1. Residual function:      F(h)  ≈ 0
%   2. Jacobian matrix:        J = ∂F/∂h
%
% The Jacobian is defined as:
%
%        J(i,j) = ∂F_i / ∂h_j
%
% where F_i is the mass balance residual at node i. We approximate this
% derivative using **forward finite differences**:
%
%        J(:,j) ≈ [F(h + ε·e_j) - F(h)] / ε
%
% where ε is a small perturbation (e.g., 1e-6), and e_j is the unit vector
% with a 1 at the j-th position.
%
% Special treatment is given to Dirichlet and Neumann boundary nodes
% to enforce known derivatives.
%
% --------------------------------------------------------------------------
%
% Inputs:
%   h           – Current pressure head profile [Nz x 1]
%   h_old       – Previous time step pressure head [Nz x 1]
%   params      – Struct with model configuration
%   top_val     – Current top BC value (flux or head)
%   bottom_val  – Current bottom BC value (flux or head)
%   source_term – Volumetric source/sink term [Nz x 1]
%
% Output:
%   J – Sparse Jacobian matrix [Nz x Nz]
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
%% =========================================================================
function J = compute_jacobian(h, h_old, params, top_val, bottom_val, source_term_fixed)

    %% === 1. Initialization ===============================================
    Nz     = params.Nz;
    eps_fd = 1e-6;                      % Perturbation magnitude
    J      = spalloc(Nz, Nz, 3*Nz);     % Preallocate sparse matrix

    %% === 2. Unperturbed Residual (Base Case) =============================
    theta_old  = theta_vgm(h_old,  params.theta_r, params.theta_s, ...
                                     params.alpha, params.n, params.m);

    theta_base = theta_vgm(h,      params.theta_r, params.theta_s, ...
                                     params.alpha, params.n, params.m);

    K_base     = K_vgm(h,          params.Ks,      params.theta_r, ...
                                     params.theta_s, params.alpha, params.n, params.m);

    % === Apply structural drainage sinks before computing residual ==========
    [source_drainage, dq_orifice, dq_spill] = drainage_sinks( ...
        h, (params.theta_s - params.theta_r), params.dz, params.dt, ... % Calculated with original h
        params.K_orifice, params.exp_orifice, ...
        params.spillway_enabled, params.c_spillway, ...
        params.h_spill, params.exp_spillway);

    % === Store analytical derivatives of drainage ========================
    dq_drainage_total = dq_orifice + dq_spill;
    if sum(dq_drainage_total) < 0 
        ttt = 1;
    end

    % === Combine total source term =================================
    source_term_value = source_term_fixed + source_drainage;

    F0 = compute_residual(h, h_old, theta_base, theta_old, K_base, ...
                          params, top_val, bottom_val, source_term_value);

    %% === 3. Loop Over Nodes to Build Jacobian ============================
    for j = 1:Nz
        h_pert       = h;
        h_pert(j)    = h(j) + eps_fd;

        theta_pert = theta_vgm(h_pert, params.theta_r, params.theta_s, ...
                                          params.alpha, params.n, params.m);

        K_pert     = K_vgm(h_pert, params.Ks, params.theta_r, ...
                                      params.theta_s, params.alpha, params.n, params.m);

        % === Apply structural drainage sinks before computing residual ==========
        [source_drainage, ~, ~] = drainage_sinks( ...
            h_pert, (params.theta_s - params.theta_r), params.dz, params.dt, ... % Calculated with perturbated h
            params.K_orifice, params.exp_orifice, ...
            params.spillway_enabled, params.c_spillway, ...
            params.h_spill, params.exp_spillway);

        % === Combine total source term =================================
        source_term_value = source_term_fixed + source_drainage;
        
        % % === Compute Perturbed Residual without drainage ==============
        F_pert = compute_residual(h_pert, h_old, theta_pert, theta_old, ...
                                  K_pert, params, top_val, bottom_val, source_term_value); % Fixed source term

        % % Skip injection at BCs
        % if j ~= 1 && j ~= Nz
        %     F_pert = F_pert - dq_drainage_total;  % Subtract analytical drainage derivative
        % end

        J(:, j) = sparse((F_pert - F0) / eps_fd);  % Finite difference (without drainage)

        %% === 4. Enforce Analytical Derivatives at BC Nodes ===============
        % if j == 1
        %     switch params.bottom_bc_type
        %         case "dirichlet"
        %             J(:, j) = sparse([ -1; zeros(Nz-1, 1) ]);
        %         case "neumann"
        %             J(:, j) = sparse([ K_base(1)/params.dz(1); zeros(Nz-1, 1) ]);
        %         case "free"
        %             J(:, j) = sparse([ -1; zeros(Nz-1, 1) ]);
        %         case "noflow"
        %             J(:, j) = sparse([ -1 / params.dz(1); zeros(Nz-1, 1) ]);
        %     end
        % end
        % 
        % if j == Nz
        %     switch params.top_bc_type
        %         case "dirichlet"
        %             J(:, j) = sparse([ zeros(Nz-1, 1); -1 ]);
        %         case "neumann"
        %             J(:, j) = sparse([ zeros(Nz-1, 1); +K_base(end)/params.dz(end) ]);
        %         case "free"
        %             J(:, j) = sparse([ zeros(Nz-1, 1); +1 ]);
        %         case "noflow"
        %             J(:, j) = sparse([ zeros(Nz-1, 1); -1 / params.dz(end) ]);
        %     end
        % end
    end
end 


