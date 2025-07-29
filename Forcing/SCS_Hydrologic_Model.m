% =========================================================================
% 🌧️ SCS-CN Hydrologic Model Coupled with Nonlinear Reservoir Routing
% =========================================================================
% Description:
%   - Simulates rainfall-runoff using SCS-CN method and nonlinear reservoir routing
%   - Supports Green Infrastructure (GI) and impervious surfaces
%   - Can read user-defined rainfall, ETP, and inflow hydrographs from Excel
%   - Outputs: hydrograph, water depth, infiltration, cumulative rainfall
% Author:
%   Marcus Nóbrega, Ph.D. | Optimized: May 2025
% =========================================================================
function [time_min, Q_total, H, infil_rate, inflow_t_min, inflow_q] = SCS_Hydrologic_Model( ...
dt, duration_min, width, length, CN, h0, n_mann, Aimp, slope, baseflow, Ks, Kr, ...
    A_GI, CN_GI, imp_drain_to_GI, forcing_path)

%% === 1. Basic Setup and Check ==========================================
nsteps = duration_min * 60 / dt;    % Total number of time steps
A = width * length;                 % Watershed area [m^2]
A_ha = A / 1e4;                     % Watershed area [ha]
nw = 1;                             % Only one watershed supported now

if Aimp + A_GI > 1
    error('Impervious and GI fractions exceed watershed area (must be <= 1).');
end

%% === 2. Load Rainfall, ETP, and Optional Inflow ========================
data = readtable(forcing_path);

rain = -table2array(data(2:end,2));           % [mm/h]
etp  = +table2array(data(2:end,3));           % [mm/h]
time_vec = table2array(data(2:end,1));       % [min]
manual_flag = table2array(data(1,6));

% DELETE
manual_flag = 0;
% Inflow Hydrograph (optional)
if manual_flag == 1
    inflow_t_min = time_vec;
    inflow_q = table2array(data(2:end,4));    % [mm/h]
else
    inflow_t_min = [];
    inflow_q = [];
end

% Ensure time steps are uniform
dt_rain = (time_vec(2) - time_vec(1));
n_rain_steps = round(time_vec(end) / dt_rain);

% Disaggregate Rainfall/ETP to Model dt
if min(rain) < 0
    rain = (-1) * rain; % Rain was entered as negative
end
if min(etp) < 0
    etp = (-1) * etp; % Rain was entered as negative
end
rain_interp = interp1(time_vec, rain, linspace(0, time_vec(end), nsteps), 'previous')';
etp_interp  = interp1(time_vec, etp,  linspace(0, time_vec(end), nsteps), 'previous')';

%% === 3. Soil Parameters and Effective Rainfall ==========================
lambda = 0.2;  % Initial abstraction ratio (e.g. 0.05 - 0.2 for SCS)
S_CN   = 25400 / CN - 254;    % SCS retention param [mm]
S_GI   = 25400 / CN_GI - 254; % GI CN [mm]
f_per  = 1 - Aimp - A_GI;
f_GI   = imp_drain_to_GI .* Aimp + A_GI;
f_imp  = Aimp;
f_GI_inflow = max(f_GI / max(A_GI, 1e-6), 0);

% GI cumulative rainfall
Pcum_GI = cumsum(rain_interp * dt / 3600 * f_GI_inflow);
EffRain_GI = (max(Pcum_GI - lambda*S_GI, 0)).^2 ./ (Pcum_GI + (1 - lambda)*S_GI);
EffRain_GI = EffRain_GI / 1000; % [m]

%% === 4. Preallocate Vectors ============================================
H  = zeros(nsteps,1);     % Water depth [m]
f  = zeros(nsteps,1);     % Infiltration rate [mm/h]
Q  = zeros(nsteps,1);     % Runoff [m^3/s]
P  = 0;                   % Precip. storage for SCS
F  = 0;                   % Cumulative infiltration
Tdry = 0;                 % Dry duration tracker
S   = S_CN;               % Current retention
emptT = 4.5 / sqrt(Ks / 25.4); % Emptying time [hr]

%% === 5. Main Loop ======================================================
for i = 1:nsteps
    % --- Infiltration Logic (with ponded water allowance) ---
    if rain_interp(i) > 0 
        % Allow infiltration if it is raining 
        P = P + rain_interp(i) * dt / 3600;  % Add rainfall to cumulative P        
        F1 = P - P^2 / (P + S);                 % SCS cumulative infiltration
        f(i) = max((F1 - F) / (dt/3600), 0);    % Instantaneous infiltration rate
        F = F1;
    elseif H(i) > 0 %  Ponding
        % No need to update P
        f(i) = f(i - 1); % Take previous infilration rate
        F1 = F1 + f(i) * dt / 3600; % SCS cumulative infiltration
    else
        f(i) = 0; % No infiltration
    end
    
    % --- Soil Retention Recovery (Kr) ---
    if rain_interp(i) == 0 && H(i) == 0
        Tdry = Tdry + dt/3600;
        if Tdry > emptT
            % Full soil recovery
            P = 0; F = 0; S = S_CN;
        else
            % Gradual recovery of retention capacity
            S = S + Kr * (S_CN - S) * dt/3600;
        end
    else
        Tdry = 0;
    end

    % --- Mass Balance for Water Depth ---
    qimp = rain_interp(i) / 1000 / 3600 * f_imp; % m/s
    runoff = (1/n_mann) * width * slope^0.5 * max(H(i) - h0, 0)^(5/3); % Pervious Areas Runoff [m3/s]
    available_flow = H(i) / dt * A; % [m3/s]
    runoff = min(runoff, available_flow); % [m3/s]

    dH = (rain_interp(i) - etp_interp(i) - f(i)) * dt / 3600 / 1000 ...
        - runoff / A  * dt ...
        + qimp * dt + EffRain_GI(min(i,end)) * f_GI;


    if i < nsteps
        H(i+1) = max(H(i) + dH,0); % Ensuring no negative ponding (ETP > ET)
    end

    Q(i) = runoff + baseflow;  % [m^3/s]
end

%% === 6. Output Formatting and Plot ======================================
time_min = (1:nsteps) * dt / 60;   % Time [min]
Q_total = Q;                       % [m^3/s]
infil_rate = f;                    % [mm/h]

% === Load Color Palettes from coloramps() =================================
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();

% Assign palette colors for consistency
c1 = pallete.red_colors(2,:);   % Flow - Red
c2 = pallete.green_colors(2,:);    % Infiltration - Green
c3 = pallete.blue_colors(1,:);  % Rainfall -Blue

% === Professional Hydrology Plot with Rainfall (Half-A4) =================
figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 21, 14.85]);  % Half-A4
tiledlayout(3,1, 'TileSpacing','compact', 'Padding','compact');

% === LaTeX + Nature Style Config ==========================================
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesFontName', 'Times New Roman'); 
set(groot, 'defaultAxesFontSize', 11); 
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% === Subplot 1: Flow and Rainfall ========================================
nexttile
yyaxis left
set(gca,'YColor','black')
plot(time_min / 1440, Q, '-', 'LineWidth', 1.8, 'Color', c1); 
ylabel('$Q$ [m$^3$/s]', 'FontSize', 12);
ylim([0, max(Q)*1.8])
yyaxis right
set(gca,'YColor', c3)
bar(time_min / 1440, rain_interp, 1, 'FaceColor', c3, 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
ylim([0,max(rain_interp)*4])
ylabel('Rainfall [mm/h]', 'FontSize', 12);

title('Hydrograph with Rainfall Forcing', 'FontWeight', 'normal', 'FontSize', 13);
grid on; box on;
set(gca, 'TickDir', 'out', 'FontSize', 11, 'LineWidth', 2);

% === Subplot 2: Infiltration Rate ========================================
nexttile
plot(time_min / 1440, infil_rate, '-', 'LineWidth', 1.8, 'Color', c2);
ylabel('Infiltration [mm/h]', 'FontSize', 12);
grid on; box on;
set(gca, 'TickDir', 'out', 'FontSize', 11, 'LineWidth', 2);

% === Subplot 3: Cumulative Values ========================================
nexttile
plot(time_min / 1440, cumsum(rain_interp)*dt/3600, '--', 'Color', c3, 'LineWidth', 2); hold on;
plot(time_min / 1440, cumsum(f)*dt/3600, '-', 'Color', c2, 'LineWidth', 2);
legend({'Cumulative Rain', 'Cumulative Infiltration'}, 'Location', 'best', 'FontSize', 10);
xlabel('Time [days]', 'FontSize', 12);
ylabel('Cumulative [mm]', 'FontSize', 12);
grid on; box on;
set(gca, 'TickDir', 'out', 'FontSize', 11, 'LineWidth', 2);

% === Export ============================================================== 
exportgraphics(gcf, 'SCS_CN_Hydrology.pdf', 'ContentType', 'vector');


end
