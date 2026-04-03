function [head_out, theta_out, flux_out, ponding_series, ...
          outlet_flux, time_series, max_ponding_depth, ...
          Se_top, Se_mid, Se_bottom, success] = Main_Solver_Function(params)
% =========================================================================
% 📂 File Location : Numerical_Solver/Main_Solver.m
%
% 🧠 MIXED-FORM RICHARDS SOLVER — Time-Stepping Engine (FUNCTION VERSION)
% -------------------------------------------------------------------------
% Description :
%   Solves the transient, one-dimensional mixed-form Richards Equation
%   using an implicit Newton-Raphson method with adaptive time stepping
%   and robust line search for convergence control.
%
% Governing Equation:
%   dθ/dt + S_s·dh/dt = d/dz [K(h)·(dh/dz + 1)] - S(z,t)
%
% Inputs:
%   params : struct with all model, soil, BC and numerical parameters.
%            Required fields (as used below) include, e.g.:
%               Nz, Tmax, dt, dt_min, dt_max, adapt_up, adapt_down,
%               n_up, n_down, max_iters, tol,
%               theta_r, theta_s, alpha, n, m, Ks, S_s, dz,
%               top_bc_type, bottom_bc_type, ...
%               surface_flux_time, surface_flux_vals,
%               bottom_flux_time,  bottom_flux_vals,
%               save_interval_min or save_interval, etc.
%
% Outputs:
%   head_out        [Nz x n_save]   pressure head [m]
%   theta_out       [Nz x n_save]   water content [-]
%   flux_out        [Nz+1 x n_save] vertical fluxes [m/s]
%   ponding_series  [1 x n_save]    surface ponding depth [m]
%   outlet_flux     [1 x n_save]    bottom boundary flux [m/s]
%   time_series     [1 x n_save]    time stamps [s]
%   seepage_flux    [1 x n_save]    seepage flux at bottom [m/s]
%   Q_orifice_total [1 x n_save]    total orifice flow [m³/m²/s]
%   Q_spillway_total[1 x n_save]    total spillway flow [m³/m²/s]
%
% Author      : Marcus Nóbrega, Ph.D.
% Updated     : Function wrapper, 2025
% =========================================================================

%% === Basic Setup & Preallocation =======================================
success = true;   % assume success unless proven otherwise

Nz   = params.Nz;
Tmax = params.Tmax;

% Save interval (if not precomputed outside)
if isfield(params, "save_interval")
    save_interval = params.save_interval;
else
    % Fallback: use save_interval_min (in minutes) or dt
    if isfield(params, "save_interval_min")
        save_interval = max(params.save_interval_min * 60, params.dt);
    else
        save_interval = params.dt;   % last resort
    end
    params.save_interval = save_interval;
end

% Time series length for saved outputs
n_save = floor(Tmax / save_interval) + 1;

% Allocate outputs
head_out        = nan(Nz,     n_save);
theta_out       = nan(Nz,     n_save);
flux_out        = nan(Nz + 1, n_save);
ponding_series  = zeros(1,    n_save);
outlet_flux     = zeros(1,    n_save);
time_series     = zeros(1,    n_save);

%%% NEW: Effective saturation at top / mid / bottom node
Se_top    = nan(1, n_save);   % top node (= surface)
Se_mid    = nan(1, n_save);   % middle node
Se_bottom = nan(1, n_save);   % bottom node

seepage_flux     = zeros(1, n_save);
Q_orifice_total  = zeros(1, n_save);
Q_spillway_total = zeros(1, n_save);

% Plotting schedule
n_plots    = 0;
plot_times = linspace(0, Tmax, n_plots);
plot_index = 1;

% Optional Delta/Gamma (if used by get_boundary_values)
if isfield(params, "Delta")
    Delta = params.Delta;
else
    Delta = [];
end
if isfield(params, "Gamma")
    Gamma = params.Gamma;
else
    Gamma = [];
end
% === Initial pressure head profile (hydrostatic by default) =============
if isfield(params, "h_init")
    % Highest priority: user explicitly provided full initial profile
    h = params.h_init;
else
    % Use soil-dependent field-capacity suction if available,
    % otherwise fall back to a generic value
    if isfield(params, "initial_suction") && ~isnan(params.initial_suction)
        p0 = params.initial_suction;    % [m], e.g. -1 m (coarse) or -3.3 m (medium/fine)
    else
        p0 = -3;                        % [m] fallback uniform suction
    end

    h = p0 * ones(1, Nz);              % uniform profile at FC (approx.)

    % Respect Dirichlet boundary conditions if specified
    if params.bottom_bc_type == "dirichlet"
        h(1) = params.bottom_bc_value;
    end
    if params.top_bc_type == "dirichlet"
        h(end) = params.top_bc_value;
    end
end


% Initial ponding depth
if params.top_bc_type == "dirichlet" && isfield(params, "top_bc_value")
    ponding_depth = max(params.top_bc_value, 0);
else
    ponding_depth = 0;
end

% 🔹 Track maximum ponding depth over the whole simulation
max_ponding_depth = ponding_depth;
%% === ⏱ Initialization ==================================================
t_end = Tmax;
nonconv_dtmin_steps = 0;   % number of consecutive time steps that failed at dt_min
t     = params.dt;   % time starts at first step
tstep = 0;
save_count = 0;
save_index = 1;
q_prev = zeros(Nz + 1, 1);
cumulative_net_flux_prev = 0;
Q_orifice  = 0;
Q_spillway = 0;
inflow_vol      = 0;
outflow_vol     = 0;
seepage_vol     = 0;
evaporation_vol = 0;
infiltration_vol = 0;
mb_error_cumulative = 0;
current_storage     = 0;   % placeholder initial storage

%% === TIME LOOP =========================================================
while t <= t_end
    %% === 💾 Store Previous Step ========================================
    h_old = h;
    h_new = h_old;
    tstep = tstep + 1;
    index_failure = 0; %#ok<NASGU>  % kept for compatibility

    ponding_prev = ponding_depth;
    [top_val, bottom_val, ponding_depth, top_bc_type_used] = ...
        get_boundary_values(h_old, params, t, ponding_prev, Delta, Gamma);

    %% === 🔁 Newton-Raphson + Line Search ===============================
    for adapt_attempt = 1:20
        converged   = false;
        % Source term interpolation (if given)
        if isfield(params, "source_times") && ~isempty(params.source_times)
            % Interpolate source term in time, then reshape to column and divide by dz
            source_term_row = interp1(params.source_times, params.source_profile', ...
                t, 'linear', 'extrap');      % 1 x Nz
            source_term     = source_term_row ./ params.dz;  % 1 x Nz  [1/s]
        else
            source_term = zeros(1, Nz);  % no distributed source
        end

        for k = 1:params.max_iters
            theta     = theta_vgm(h_new, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
            theta_old = theta_vgm(h_old, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
            K         = K_vgm(h_new, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);

            Q_orifice_prev  = Q_orifice;
            Q_spillway_prev = Q_spillway;

            % === Apply structural drainage sinks before residual =========
            [source_drainage, ~, ~, Q_orifice, Q_spillway] = drainage_sinks( ...
                h_new, (theta - params.theta_r), params.dz, params.dt, ...
                params.K_orifice, params.exp_orifice, ...
                params.spillway_enabled, params.c_spillway, ...
                params.h_spill, params.exp_spillway);

            % Combine total source term
            source_term_value = source_term + source_drainage;

            [F, q] = compute_residual(h_new, h_old, theta, theta_old, K, ...
                params, top_val, bottom_val, source_term_value, top_bc_type_used);

            J = compute_jacobian(h_new, h_old, params, top_val, bottom_val, ...
                source_term, top_bc_type_used);

            delta = -(J \ F')';   % Newton step

            %% === 🔍 Line Search: Robust Convergence =====================
            lambda     = 1.0;   % full Newton step
            lambda_min = 1e-6;
            beta       = 0.5;   % backtracking factor
            eta        = 1e-4;  % Armijo criterion
            res_norm_0 = norm(F);

            success = false;
            for ls_iter = 1:100
                h_trial = h_new + lambda * delta;

                theta_trial = theta_vgm(h_trial, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
                K_trial     = K_vgm(h_trial, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);

                % (We reuse source_term_value from h_new; if you want full
                % consistency, you could recompute drainage here for h_trial.)

                F_trial = compute_residual(h_trial, h_old, theta_trial, theta_old, ...
                    K_trial, params, top_val, bottom_val, source_term_value, top_bc_type_used);

                if norm(F_trial) < (1 - eta * lambda) * res_norm_0
                    h_new  = h_trial;
                    success = true;
                    break;
                else
                    lambda = lambda * beta;
                    if lambda < lambda_min
                        break;
                    end
                end
            end
            
            %% === ✅ Check Convergence + Mass Balance ====================
            print_error = 0;
            [mb_error, ~, ~, current_storage] = mass_balance_check( ...
                h_old, h_new, ponding_prev, ponding_depth, q_prev, q, t, ...
                params.dt, mb_error_cumulative, params, ...
                cumulative_net_flux_prev, Q_orifice, Q_spillway, ...
                Q_orifice_prev, Q_spillway_prev, top_bc_type_used, print_error);

            mb_tol = 1e-3;

            % if norm(delta) < params.tol && norm(F) < params.tol && abs(mb_error) < mb_tol
            %     converged = true;
            %     break;
            % end
            if norm(delta) < params.tol && norm(F) < params.tol 
                converged = true;
                break;
            end
        end

        %% === 📊 Mass Balance (verbose) =================================
        if tstep > 1 && converged
            print_error = 0;
            [mb_error, mb_error_cumulative, cumulative_net_flux] = mass_balance_check( ...
                h_old, h_new, ponding_prev, ponding_depth, q_prev, q, t, ...
                params.dt, mb_error_cumulative, params, ...
                cumulative_net_flux_prev, Q_orifice, Q_spillway, ...
                Q_orifice_prev, Q_spillway_prev, top_bc_type_used, print_error);
            cumulative_net_flux_prev = cumulative_net_flux;
        else
            mb_error            = 0;
            mb_error_cumulative = 0;
        end

        %% === ⏱ Adaptive Timestep Control ==============================
        if converged
            % Normal adaptive behaviour when we *did* converge
            if k <= params.n_up
                params.dt = min(params.dt * params.adapt_up, params.dt_max);
            elseif k >= params.n_down
                params.dt = max(params.dt * params.adapt_down, params.dt_min);
            end

            if h_new(end) > 1 % Larger than 100 cm, we break
                warning('Ponding depth larger than 100 cm')
                success = false;
                return
            end

            % On success, reset the non-convergence counter
            nonconv_dtmin_steps = 0;

            break;  % leave adapt_attempt loop (time step accepted)

        else
            % Not converged for this Newton iteration: shrink dt
            params.dt = max(params.dt * params.adapt_down, params.dt_min);

            % If this was the last allowed adapt_attempt for this time step,
            % we give up on this time step and handle it *after* the for-loop.
            if adapt_attempt == 20
                break;
            end
        end
    end  % <-- end of adapt_attempt loop

    % ---------------------------------------------------------------------
    % If we exit the adapt_attempt loop WITHOUT convergence, decide what to do
    % ---------------------------------------------------------------------
    if ~converged
        if params.dt < params.dt_min
            % We are at the minimum time step and still did not converge
            nonconv_dtmin_steps = nonconv_dtmin_steps + 1;
            warning('Time step at t = %.2f s did not converge with dt_min (attempt %d of 10).', ...
                t, nonconv_dtmin_steps);

            if nonconv_dtmin_steps > 2
                warning('❌ Simulation stopped: more than 10 consecutive non-convergent steps at dt_min.');
                success = false;
                head_out        = [];
                theta_out       = [];
                flux_out        = [];
                ponding_series  = [];
                outlet_flux     = [];
                time_series     = [];
                max_ponding_depth = [];
                Se_top    = [];    
                Se_mid    = [];    
                Se_bottom = [];   
                return;
            else
                % Try this *same* time step again in the next while-iteration
                % (t is not advanced, h remains h_old)
                continue;   % go to next iteration of the while t <= t_end loop
            end

        else
            % We failed before even reaching dt_min: treat as hard failure
            warning('❌ Simulation failed: unable to converge before reaching dt_min.');
            success = false;
            head_out        = [];
            theta_out       = [];
            flux_out        = [];
            ponding_series  = [];
            outlet_flux     = [];
            time_series     = [];
            max_ponding_depth = []; 
            Se_top    = [];         
            Se_mid    = [];          
            Se_bottom = [];         
            return;
        end
    end
    %% === 📈 Plotting (optional) ========================================
    if (plot_index <= numel(plot_times)) && (t >= plot_times(plot_index) || tstep == 1)
        % You can wrap this in a plotting flag if needed:
        if tstep == 1
            figure('Color','w', 'Units','inches', 'Position',[1 1 7.5 7]);
            tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
        end
        clf;

        sim_percentage = 100 * t / Tmax;
        plot_soil_profiles(h_new, h_old, params, t, tstep, @theta_vgm, @K_vgm, q);

        sgtitle(sprintf(['Time = %.2f min — Step = %d — Simulation = %.2f days ' ...
            '(%.1f%% Complete)'], ...
            t / 60, tstep, t / 86400, sim_percentage), ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'FontName', 'Helvetica', 'Interpreter', 'none');

        drawnow;
        plot_index = plot_index + 1;
    end

    %% === 💾 Save Outputs ===============================================

    % Accumulate volumes (signs depend on q convention; kept as in your code)
    inflow_vol       = -q(end) * params.dt + inflow_vol;       % [m]
    outflow_vol      = -q(1)   * params.dt + outflow_vol;      % [m]
    seepage_vol      = -min(q(1), 0) * params.dt + seepage_vol;
    infiltration_vol = -min(q(end), 0) * params.dt + infiltration_vol;
    evaporation_vol  =  max(q(end), 0) * params.dt + evaporation_vol;
    final_storage    = current_storage; %#ok<NASGU>

    current_save_time = params.save_interval * (save_count);

    if t >= current_save_time || tstep == 1
        save_count = save_count + 1;

        theta_now = theta_vgm(h_new, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
        K_now     = K_vgm(h_new, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
        [~, q_now] = compute_residual(h_new, h_old, theta_now, theta_now, K_now, ...
            params, top_val, bottom_val, source_term, top_bc_type_used);

        head_out(:, save_count)    = h_new(:);
        theta_out(:, save_count)   = theta_now(:);
        flux_out(:, save_count)    = q_now(:);
        ponding_series(save_count) = ponding_depth;
        outlet_flux(save_count)    = q_now(1);
        time_series(save_count)    = t;
        %%% NEW: Effective saturation profile and extract top/mid/bottom
        % Se = (theta - theta_r) / (theta_s - theta_r)
        Se_profile = (theta_now - params.theta_r) ./ (params.theta_s - params.theta_r);
    
        idx_bottom = 1;                  % z(1) is bottom (deepest node)
        idx_top    = Nz;                 % z(end) is surface node
        idx_mid    = round((Nz + 1)/2);  % middle node
    
        Se_bottom(save_count) = Se_profile(idx_bottom);
        Se_mid(save_count)    = Se_profile(idx_mid);
        Se_top(save_count)    = Se_profile(idx_top);

        if params.bottom_bc_type == "noflow"
            seepage_flux(save_count) = max(0, q_now(1));
        else
            seepage_flux(save_count) = q_now(1);
        end

        Q_orifice_total(save_count)  = sum(Q_orifice(:));
        Q_spillway_total(save_count) = sum(Q_spillway(:));

        % 🔹 Update maximum ponding depth
        if ponding_depth > max_ponding_depth
            max_ponding_depth = ponding_depth;
        end

        save_index = save_index + 1;
    end

    %% === ✅ Accept Time Step ===========================================
    h = h_new;

    % Snap dt to hit the next save time exactly (optional but nice)
    next_save_time = params.save_interval * save_count;
    if t < next_save_time && (t + params.dt) > next_save_time
        params.dt = next_save_time - t;
    end

    t      = t + params.dt;
    q_prev = q;
end

%% ✅ Completion Message ==================================================
if success
    disp('✅ Simulation completed.');
end

end
