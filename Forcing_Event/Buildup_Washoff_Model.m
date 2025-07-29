function [tint, C, L, Bw, Ww, B, W, washed_tot] = Buildup_Washoff_Model(C1, C2, C3, C4, ADD, A_km2, Qin, tint)
% =========================================================================
% 🚰 BUILDUP-WASHOFF MODEL — Urban Water Quality Simulation
% =========================================================================
% Description:
%   Simulates time-varying pollutant buildup and washoff in a watershed
%   receiving stormwater runoff. Considers antecedent dry days, hydrologic
%   inputs, and pollutant transport dynamics.
%
% Inputs:
%   C1      - Buildup coefficient [kg/ha]
%   C2      - Buildup exponent [1/days]
%   C3      - Washoff coefficient [(mm/h)^(-C4)(h^-1)]
%   C4      - Washoff exponent [-]
%   ADD     - Antecedent dry days [days]
%   A_km2   - Drainage area [km^2]
%   Qin     - Inflow hydrograph [m^3/s]
%   tint    - Time vector [min]
%
% Outputs:
%   tint       - Time vector [min]
%   C          - Pollutant concentration [mg/L]
%   L          - Pollutant load [kg/s]
%   Bw         - Total buildup mass [kg]
%   Ww         - Washoff mass rate [kg/h]
%   B          - Buildup per area [kg/ha]
%   W          - Washoff per area [kg/ha/h]
%   washed_tot - Total pollutant washed off [kg]
% =========================================================================

% === Preprocessing =======================================================
Qin  = Qin(:);      % Ensure column vector
nt   = length(tint);
tint = tint(:);     % Time [min]
dt_h = (tint(2)-tint(1))/60;  % Time step [h]

% === Derived Parameters ==================================================
A_m2 = A_km2 * 1e6;           % [m^2]
A_ha = A_m2 / 1e4;            % [ha]
q    = Qin / A_m2 * 1000 * 3600;  % Flow rate [mm/h]

% === Initialization ======================================================
B  = zeros(nt,1);   % Buildup [kg/ha]
W  = zeros(nt,1);   % Washoff rate [kg/ha/h]
Bw = zeros(nt,1);   % Buildup mass [kg]
Ww = zeros(nt,1);   % Washoff mass rate [kg/h]
C  = zeros(nt,1);   % Concentration [mg/L]
L  = zeros(nt,1);   % Load [kg/s]

% Initial buildup based on ADD
B(1)  = C1 * (1 - exp(-C2 * ADD));
W(1)  = max(C3 * q(1)^C4 * B(1), B(1)/dt_h);
Bw(1) = B(1) * A_ha;
Ww(1) = W(1) * A_ha;
C(1)  = (Qin(1) > 0) * Ww(1)/3600/Qin(1)*1e6;
L(1)  = C(1) * Qin(1) / 1000;  % [kg/s]

% === Dry Weather Buildup Initialization ==================================
dry_days = 0;
B_dry_previous = B(1);
dB = 0;

% === Cumulative Values for Mass Balance ==================================
Bsum = zeros(nt,1);
Wsum = zeros(nt,1);
Bsum(1) = B(1);

% === Main Loop ==========================================================
for i = 2:nt
    % Adaptive buildup time step
    max_attempts = 10;
    min_dt_h = dt_h / (2^max_attempts);
    dt_try = dt_h;
    stable = false;

    while ~stable && dt_try >= min_dt_h
        B_trial = B(i-1) - W(i-1)*dt_try + dB;
        if B_trial >= 0
            B(i) = B_trial;
            stable = true;
        else
            dt_try = dt_try / 2;
        end
    end

    % Washoff rate calculation
    W(i) = max(C3 * q(i)^C4 * B(i), B(i)/dt_try);
    Wsum(i) = Wsum(i-1) + W(i) * dt_try;
    
    Bw(i) = B(i) * A_ha;
    Ww(i) = W(i) * A_ha;

    % Dry or wet condition handling
    if Qin(i) == 0
        C(i) = 0;
        dry_days = dry_days + dt_try / 24;
        B_new = C1 * (1 - exp(-C2 * dry_days));
        dB = B_new - B_dry_previous;
        B_dry_previous = B_new;
        Bsum(i) = Bsum(i-1) + max(dB,0);
    else
        C(i) = Ww(i)/3600/Qin(i)*1e6;
        dry_days = 0;
        dB = 0;
        B_dry_previous = B(i);
        Bsum(i) = Bsum(i-1);
    end

    % Load computation
    L(i) = C(i) * Qin(i) / 1000;
end

% === Final Mass Balance ================================================
B_tot = Bsum(end)*A_ha;
W_tot = Wsum(end)*A_ha;
mass_error = B_tot - W_tot;
washed_tot = W_tot;

fprintf('\n--- MASS BALANCE CHECK ---\n');
fprintf('Total Buildup: %.2f kg\n', B_tot);
fprintf('Total Washoff: %.2f kg\n', W_tot);
fprintf('Error        : %.2f kg (%.2f%%)\n', mass_error, 100*mass_error/max(B_tot,1e-6));

% === Plotting ===========================================================
plot_buildup_washoff(tint, Qin, C, Bsum, Wsum, L, A_ha);
end

function plot_buildup_washoff(t, Qin, C, Bsum, Wsum, L, A_ha)
    figure('Units','inches','Position',[2, 0, 6.5, 7]);
    set(groot,'defaultTextInterpreter','latex');
    fs = 13;

    % --- Subplot 1: Hydrograph and Concentration ---
    subplot(3,1,1); hold on; grid on;
    plot(t, Qin, 'b', 'LineWidth', 1.8);
    yyaxis right;
    plot(t, C, 'r', 'LineWidth', 1.8);
    ylabel('$C$ [mg/L]', 'FontSize', fs);
    yyaxis left;
    ylabel('$Q$ [m$^3$/s]', 'FontSize', fs);
    xlabel('Time [min]', 'FontSize', fs);
    legend('Discharge','Concentration','Location','best','Interpreter','latex');

    % --- Subplot 2: Cumulative Masses ---
    subplot(3,1,2); hold on; grid on;
    plot(t, Bsum*A_ha, 'k', 'LineWidth', 2);
    plot(t, Wsum*A_ha, 'g', 'LineWidth', 2);
    ylabel('Cumulative Mass [kg]', 'FontSize', fs);
    xlabel('Time [min]', 'FontSize', fs);
    legend('Buildup','Washoff','Location','best','Interpreter','latex');

    % --- Subplot 3: Load and Concentration ---
    subplot(3,1,3); hold on; grid on;
    plot(t, L, 'm', 'LineWidth', 1.8);
    yyaxis right;
    plot(t, C, 'r', 'LineWidth', 1.8);
    ylabel('$C$ [mg/L]', 'FontSize', fs);
    yyaxis left;
    ylabel('Load [kg/s]', 'FontSize', fs);
    xlabel('Time [min]', 'FontSize', fs);
    legend('Load','Concentration','Location','best','Interpreter','latex');

    exportgraphics(gcf,'Water_Quality_Model.pdf','ContentType','vector');
end
