%% =========================================================================
% 📂 File Location : Numerical_Solver/Main_Solver.m
%
% 🧠 MIXED-FORM RICHARDS SOLVER — Time-Stepping Engine
% -------------------------------------------------------------------------
% Description :
%   Solves the transient, one-dimensional mixed-form Richards Equation
%   using an implicit Newton-Raphson method with adaptive time stepping
%   and robust line search for convergence control.
%
% 💡 Governing Equation:
%   dθ/dt + S_s·dh/dt = d/dz [K(h)·(dh/dz + 1)] - S(z,t)
%
%   where:
%       θ(h)  : volumetric moisture content
%       h     : pressure head
%       K(h)  : unsaturated hydraulic conductivity
%       S_s   : specific storage [1/m]
%       S     : source/sink term [1/s]
%
% 🔁 Iterative Strategy:
%   • Newton-Raphson method
%   • Adaptive time stepping
%   • Line search for robust damping
%
% Author      : Marcus Nóbrega, Ph.D.
% Updated     : May 2025

%% =========================================================================

% === ⏱ Initialization ====================================================
t_end = params.Tmax;
t = params.dt;
tstep = 0;
save_count = 0;
save_index = 1;
q_prev = zeros(params.Nz + 1, 1);
cumulative_net_flux_prev = 0;
Q_orifice = 0;
Q_spillway = 0;
inflow_vol = 0;
outflow_vol = 0;
seepage_vol = 0;
evaporation_vol = 0;
infiltration_vol = 0;
time_step = zeros(size(head_out,2));

% --- Clogging capture setup (records only if clogging is active) ---------
clogging_active = isfield(params,'enable_clogging') && islogical(params.enable_clogging) ...
                  && params.enable_clogging && exist('clog','var') ...
                  && isfield(clog,'mask') && ~isempty(clog.mask);

if clogging_active
    mask_idx   = find(clog.mask);                          % all nodes that can clog
    nsave_est  = ceil(params.Tmax/params.save_interval) + 2;
    theta_s_out      = nan(numel(mask_idx), nsave_est);    % θs for cloggable nodes
    Ks_out           = nan(numel(mask_idx), nsave_est);    % Ks for cloggable nodes
    Ks_node1_series  = nan(1, nsave_est);                  % optional proxy at node 1
else
    mask_idx         = [];
    theta_s_out      = [];
    Ks_out           = [];
    Ks_node1_series  = [];
end
saved_this_step = false;    % toggled when we save a time slice


% === TIME LOOP ===========================================================
while t <= t_end

    %% === 💾 Store Previous Step ==========================================
    h_old = h;
    h_new = h_old;
    tstep = tstep + 1;
    index_failure = 0;

    if t >= 2765700*60
        ttt = 1;
    end


    ponding_prev   = ponding_depth;
    [top_val, bottom_val, ponding_depth, top_bc_type_used] = get_boundary_values(h_old, params, t, ponding_prev, Delta, Gamma);
    
    %% === 🔁 Newton-Raphson + Line Search ================================
    for adapt_attempt = 1:20
        converged = false;
        source_term = interp1(params.source_times, params.source_profile', t, 'linear', 'extrap') ./ params.dz; % Original Source Term [1/sec]
        for k = 1:params.max_iters
            theta     = theta_vgm(h_new, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
            theta_old = theta_vgm(h_old, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
            K         = K_vgm(h_new, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);

            Q_orifice_prev = Q_orifice; Q_spillway_prev = Q_spillway;
            % === Apply structural drainage sinks before computing residual ==========
            [source_drainage, ~, ~, Q_orifice, Q_spillway] = drainage_sinks( ...
                h_new, (theta - params.theta_r), params.dz, params.dt, ... % Calculated with h
                params.K_orifice, params.exp_orifice, ...
                params.spillway_enabled, params.c_spillway, ...
                params.h_spill, params.exp_spillway);

            % === Combine total source term =================================
            source_term_value = source_term + source_drainage;

            [F, q] = compute_residual(h_new, h_old, theta, theta_old, K, params, top_val, bottom_val, source_term_value, top_bc_type_used); % Original Residual with Drainage from h
            J      = compute_jacobian(h_new, h_old, params, top_val, bottom_val, source_term, top_bc_type_used); % Here we must perturbate drainage inside

            delta = -(J \ F')';

            %% === 🔍 Line Search: Robust Convergence =====================
            lambda     = 1.0;  % Full Newton step
            lambda_min = 1e-8;
            beta       = 0.25;  % Backtracking factor
            eta        = 1e-4; % Armijo criterion
            res_norm_0 = norm(F);

            success = false;
            for ls_iter = 1:100
                h_trial = h_new + lambda * delta;

                theta_trial = theta_vgm(h_trial, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
                K_trial     = K_vgm(h_trial, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);

                % === Apply structural drainage sinks before computing residual ==========
                % [source_drainage, Q_orifice, Q_spillway] = drainage_sinks( ...
                %     h_trial, (theta - params.theta_r), params.dz, params.dt, ... % Calculated with h_trial
                %     params.K_orifice, params.exp_orifice, ...
                %     params.spillway_enabled, params.c_spillway, ...
                %     params.h_spill, params.exp_spillway);
                %
                % % === Combine total source term =================================
                % source_term_value = source_term + source_drainage;

                % F_trial     = compute_residual(h_trial, h_old, theta_trial, theta_old, K_trial, params, top_val, bottom_val, source_term_value);
                F_trial     = compute_residual(h_trial, h_old, theta_trial, theta_old, K_trial, params, top_val, bottom_val, source_term_value, top_bc_type_used);

                
                % Regular Line-Search Approach
                threshold_it = eps;
                if norm(F_trial) < threshold_it || norm(F_trial) < (1 - eta * lambda) * res_norm_0 
                    h_new = h_trial;
                    success = true;
                    break;
                else
                    lambda = lambda * beta;
                    if lambda < lambda_min
                        break;
                    end
                end
            end

            %% === ✅ Check Convergence + Mass Balance =====================
            print_error = 0;
            [mb_error, ~, ~, current_storage] = mass_balance_check( ...
                h_old, h_new, ponding_prev, ponding_depth, q_prev, q, t, params.dt, mb_error_cumulative, params, cumulative_net_flux_prev, Q_orifice, Q_spillway, Q_orifice_prev, Q_spillway_prev, top_bc_type_used, print_error);
            
            mb_tol = 1e-4;  % Set mass balance tolerance [adjust as needed]
            
            if norm(delta) < params.tol && norm(F) < params.tol && abs(mb_error) < mb_tol
                converged = true;
                break;
            end
        end

        %% === 📊 Mass Balance ================================================
        if tstep > 1 & converged
            print_error = 1;
            [mb_error, mb_error_cumulative, cumulative_net_flux] = mass_balance_check( ...
                h_old, h_new, ponding_prev, ponding_depth, q_prev, q, t, params.dt, mb_error_cumulative, params, cumulative_net_flux_prev, Q_orifice, Q_spillway, Q_orifice_prev, Q_spillway_prev, top_bc_type_used, print_error);
        else
            mb_error = 0;
            mb_error_cumulative = 0;
        end

        %% === ⏱ Adaptive Timestep Control ================================
        if converged
            if k <= params.n_up
                params.dt = min(params.dt * params.adapt_up, params.dt_max);
            elseif k >= params.n_down
                params.dt = max(params.dt * params.adapt_down, params.dt_min);
            end
            break;
        else
            params.dt = max(params.dt * params.adapt_down, params.dt_min);
            if params.dt <= params.dt_min
                warning('❌ Simulation failed: minimum dt reached without convergence.');
                break
            end
        end
    end

    %% === 📈 Plotting ====================================================
    if t >= plot_times(plot_index) || tstep == 1
        if tstep == 1
            figure('Color','w', 'Units','inches', 'Position',[1 1 7.5 7]);
            tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
        end
        clf;
        
        % Add this block:
        sim_percentage = 100 * t / params.Tmax;
    
        plot_soil_profiles(h_new, h_old, params, t, tstep, @theta_vgm, @K_vgm, q);
    
        % Set title on first tile 
        sgtitle(sprintf('Time = %.2f min — Step = %d — Simulation = %.2f days (%.1f%% Complete)', ...
            t / 60, tstep, t / 86400, sim_percentage), ...
            'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'Interpreter', 'none');

    
        drawnow;
        plot_index = plot_index + 1;
    end

    % clf;
    % plot_soil_profiles(h_new, h_old, params, t, tstep, @theta_vgm, @K_vgm, @compute_residual);
    % drawnow;
    % pause(0.0001)
    %% === 💾 Save Outputs ===============================================

    % Computing Current Inflow/Outflow Fluxes (For Free-Drainage and Neumann Only)
    inflow_vol = -q(end) * params.dt + inflow_vol; % [m]
    outflow_vol = -q(1) * params.dt + outflow_vol; % [m]
    seepage_vol = -min(q(1), 0) * params.dt + seepage_vol; % [m]
    infiltration_vol = -min(q(end), 0) * params.dt + infiltration_vol; % [m]
    evaporation_vol = max(q(end), 0) * params.dt + evaporation_vol; % [m]
    final_storage = current_storage;
    
    current_save_time = params.save_interval * (save_count);
    if t >= current_save_time || tstep == 1
        save_count = save_count + 1;

        theta_now = theta_vgm(h_new, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
        K_now     = K_vgm(h_new, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
        [~, q_now] = compute_residual(h_new, h_old, theta_now, theta_now, K_now, ...
            params, top_val, bottom_val, source_term, top_bc_type_used);

        head_out(:, save_count)     = h_new(:);
        theta_out(:, save_count)    = theta_now(:);
        flux_out(:, save_count)     = q_now(:);
        ponding_series(save_count)  = ponding_depth;
        outlet_flux(save_count)     = q_now(1);
        time_series(save_count)     = t;        

        % overall_mb_error = (inflow_vol - final_storage - outflow_vol)/inflow_vol*100

        time_step(save_count)       = params.dt;
        if params.bottom_bc_type == "noflow"
            seepage_flux(save_count) = max(0, q_now(1));
        else
            seepage_flux(save_count) = q_now(1);
        end

        Q_orifice_total(save_count)  = sum(Q_orifice(:));   % Total orifice flow [m³/m²/s]
        Q_spillway_total(save_count) = sum(Q_spillway(:));  % Total spillway flow [m³/m²/s]


        save_index = save_index + 1;
        % flag: this time slice was saved; record θs/Ks after clogging update
        saved_this_step = true;
    end



    %% === ✅ Accept Time Step ===========================================
    h = h_new;

    % --- Only update clogging if the flag **exists and is true** ---
    if isfield(params,'enable_clogging') && islogical(params.enable_clogging) && params.enable_clogging
        [params, clog] = update_clogging(params, clog, q, params.dt, t);
    end

    % --- If clogging is active and we saved this step, record θs and Ks now ---
    if clogging_active && saved_this_step
        theta_s_out(:, save_count)     = params.theta_s(mask_idx);
        Ks_out(:, save_count)          = params.Ks(mask_idx);
        if ~isempty(Ks_node1_series)
            Ks_node1_series(1, save_count) = params.Ks(1);  % optional proxy
        end
        saved_this_step = false;
    end

    % === ⏳ Snap Time Step to Next Save Interval ============================
    next_save_time = params.save_interval * save_count;

    % Check if current step would overshoot the next save time
    if t < next_save_time && (t + params.dt) > next_save_time
        params.dt = next_save_time - t;  % Reduce dt to land exactly on save point
    end

    t = t + params.dt;
    q_prev = q;
end

%% ✅ Completion Message
disp('✅ Simulation completed.');