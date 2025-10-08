function [params, clog] = update_clogging(params, clog, q, dt, t)
% UPDATE_CLOGGING — evolve θs and Ks from cumulative loading & maintenance
% q sign convention in your solver:
%   In your main, inflow_vol = -q(end)*dt;  % so negative q(end) = downward/infiltration
% We’ll define: downward_flux = max(-q_cellface, 0)

    % ---- 1) Maintenance check ----
    if ~isempty(clog.next_maint_t) && t >= clog.next_maint_t
        % Partial recovery toward clean state
        restored_theta_s = clog.phi_min + clog.eta .* (clog.theta_s0 - clog.phi_min);
        params.theta_s   = max(restored_theta_s, params.theta_r + 1e-6);

        % Ks from porosity–permeability link
        ratio            = params.theta_s ./ clog.theta_s0;
        params.Ks        = clog.Ks0 .* (max(ratio,1e-6)).^clog.mK;

        % Reset the loading “clock”
        clog.Ic(:) = 0;

        % Schedule next maintenance
        if ~isempty(clog.maintenance_period_days)
            clog.next_maint_t = t + clog.maintenance_period_days*86400;
        else
            % consume one explicit date
            if ~isempty(clog.maintenance_dates_sec)
                idx = find(clog.maintenance_dates_sec > t, 1, 'first');
                if isempty(idx), clog.next_maint_t = []; else, clog.next_maint_t = clog.maintenance_dates_sec(idx); end
            else
                clog.next_maint_t = [];
            end
        end
    end

    % ---- 2) Accumulate loading Ic (only where mask=1) ----
    % Choose what flux drives clogging. For pavements, the controlling layer is near the top.
    % Use the face above each masked cell. We can map from node fluxes to cells simply:
    %   q_top_face = q(end); (top boundary). If masking a few top cells, apply same driver.
    q_down = max(-q(end), 0);  % downward-only (m/s)
    clog.Ic(clog.mask) = clog.Ic(clog.mask) + q_down * dt;  % add meters of load

    % ---- 3) Update θs and Ks from the paper’s relations ----
    % θs(Ic) = φ_min + (φ0 − φ_min)*exp(−γ Ic)
    theta_s_new = clog.phi_min + (clog.theta_s0 - clog.phi_min) .* exp(-clog.gamma .* clog.Ic);

    % Apply only where mask=1, leave other layers unchanged
    theta_s_eff = params.theta_s;
    theta_s_eff(clog.mask) = max(theta_s_new(clog.mask), params.theta_r(clog.mask) + 1e-6);
    params.theta_s = theta_s_eff;

    % Ks = Ks0 * (θs/φ0)^{mK}
    ratio = params.theta_s ./ clog.theta_s0;
    Ks_new = clog.Ks0 .* (max(ratio,1e-6)).^clog.mK;
    params.Ks(clog.mask) = Ks_new(clog.mask);
end
