function [time_series, ponding_series, head_out, flux_out] = Main_Solver_Function(media_params,baseline_path)
% =========================================================================
% 📂 Mixed-Form Richards Solver — Function Version
% =========================================================================
% Description:
%   Modular solver function for 1D mixed-form Richards equation.
%   Adapted to support automatic calibration via function calls.
%
% Inputs:
%   - media_params: Simulation and soil parameter structure
%
% Outputs:
%   - time_series: Time vector
%   - ponding_series: Ponding depth over time
%   - head_out: Pressure head profiles over time
%   - flux_out: Fluxes over time
%
% =========================================================================

load(baseline_path); % Loading input data

% Adjust Media Params
% === 3. MULTILAYER SOIL PROPERTIES ======================================
media_thicknesses = [0.203, 0.102, 0.0254];  % Bottom to top

media_props = struct( ...
    'alpha',   [media_params(1:3)], ...       % van Genuchten alpha (1/m)
    'n',       [media_params(4:6)], ...       % van Genuchten n
    'theta_r', [media_params(7:9)], ...       % Residual water content
    'theta_s', [media_params(10:12)], ...     % Saturated water content
    'S_s',     [media_params(13:15)], ...     % Specific storage (1/m)
    'Ks',      [media_params(16:18)]);        % Saturated hydraulic conductivity (m/s)


media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];
n_layers = length(media_thicknesses);

media_id = zeros(params.Nz, 1);
for i = 1:params.Nz
    zi = params.z(i);
    for j = 1:n_layers
        if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
            media_id(i) = j;
            break;
        elseif zi == media_interfaces(end)
            media_id(i) = n_layers;
        end
    end
end

params.alpha   = media_props.alpha(media_id);
params.n       = media_props.n(media_id);
params.m       = 1 - 1 ./ params.n;
params.theta_r = media_props.theta_r(media_id);
params.theta_s = media_props.theta_s(media_id);
params.S_s     = media_props.S_s(media_id);
params.Ks      = media_props.Ks(media_id);


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

% === TIME LOOP ===========================================================
while t <= t_end

    %% === 💾 Store Previous Step ==========================================
    h_old = h;
    h_new = h_old;
    tstep = tstep + 1;
    index_failure = 0;

    
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
            lambda_min = 1e-6;
            beta       = 0.5;  % Backtracking factor
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
                if norm(F_trial) < (1 - eta * lambda) * res_norm_0
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
            [mb_error, ~, ~] = mass_balance_check( ...
                h_old, h_new, ponding_prev, ponding_depth, q_prev, q, t, params.dt, mb_error_cumulative, params, cumulative_net_flux_prev, Q_orifice, Q_spillway, Q_orifice_prev, Q_spillway_prev, top_bc_type_used, print_error);
            
            if tstep == 1
                mb_error = 0;
            end
            mb_tol = 1e-3;  % Set mass balance tolerance [adjust as needed]

            if max(h_old) > 10
                ttt  = 1;
            end

            if norm(delta) < params.tol && norm(F) < params.tol && abs(mb_error) < mb_tol
                converged = true;
                break;
            end
        end

        %% === 📊 Mass Balance ================================================
        if tstep > 1 & converged
            print_error = 0;
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
                % warning('❌ Simulation failed: minimum dt reached without convergence.');
                error('❌ Simulation failed: minimum dt reached without convergence.');
                break
            end
        end
    end

    %% === 📈 Plotting ====================================================
    % if t >= plot_times(plot_index) || tstep == 1
    %     if tstep == 1
    %         figure('Color','w', 'Units','inches', 'Position',[1 1 7.5 7]);
    %         tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
    %     end
    %     clf;
    %     plot_soil_profiles(h_new, h_old, params, t, tstep, @theta_vgm, @K_vgm, q);
    %     drawnow;
    %     plot_index = plot_index + 1;
    % end

    % clf;
    % plot_soil_profiles(h_new, h_old, params, t, tstep, @theta_vgm, @K_vgm, @compute_residual);
    % drawnow;
    % pause(0.0001)
    %% === 💾 Save Outputs ===============================================
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
        top_flux(1,save_count)      = q_now(end);
        time_series(save_count)     = t;

        if params.bottom_bc_type == "noflow"
            seepage_flux(save_count) = max(0, q_now(1));
        else
            seepage_flux(save_count) = q_now(1);
        end

        Q_orifice_total(save_count)  = sum(Q_orifice(:));   % Total orifice flow [m³/m²/s]
        Q_spillway_total(save_count) = sum(Q_spillway(:));  % Total spillway flow [m³/m²/s]


        save_index = save_index + 1;
    end



    %% === ✅ Accept Time Step ===========================================
    h = h_new;

    % === ⏳ Snap Time Step to Next Save Interval ============================
    next_save_time = params.save_interval * save_count;

    % Check if current step would overshoot the next save time
    if t < next_save_time && (t + params.dt) > next_save_time
        params.dt = next_save_time - t;  % Reduce dt to land exactly on save point
    end

    t = t + params.dt;
    q_prev = q;
end