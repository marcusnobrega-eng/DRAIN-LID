%% 3 x 4 surface plots of performance vs. Ld and area ratio for 3 soils
clear; clc;

%% ===================== USER SETTINGS ====================================
% Database files (one per soil)
db_files   = { ...
    'DRAIN_LID_database_Example2_Sand.mat', ...
    'DRAIN_LID_database_Example2_LoamySand.mat', ...
    'DRAIN_LID_database_Example2_SandyLoam.mat'};

% Number of discrete colors per variable
nColors_eta  = 10;   % Peak flow reduction
nColors_dt   = 10;   % Time-to-peak delay
nColors_det  = 10;   % Detention time
nColors_hmax = 10;   % Maximum ponding depth (used only to define ticks)

% Colorbar limits (min, max) for each variable
eta_clim  = [ 0 100];   % [%]       Peak flow reduction
dt_clim   = [ 1 900];   % [min]     Time-to-peak delay (start at 1 for log scale)
det_clim  = [10 40];    % [h]       Detention time
hmax_clim = [ 0 0.15];   % [m]       Maximum ponding depth

% Number of ticks for Ld axis (smooth)
nTicks_Ld = 5;
% ========================================================================

% --- LaTeX + font defaults ---
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultTextFontSize',14);

nSoils = numel(db_files);

%% --------- Pre-scan to get global x/y ranges and common ticks ----------
Ld_all = [];
AR_all = [];

for iSoil = 1:nSoils
    S = load(db_files{iSoil}, 'Metrics');
    Metrics = S.Metrics;
    Ld_all = [Ld_all, [Metrics.Ld_m]];
    AR_all = [AR_all, [Metrics.area_ratio]];
end

Ld_all = unique(Ld_all);
AR_all = unique(AR_all);

Ld_min_global = min(Ld_all);
Ld_max_global = max(Ld_all);

AR_min_global = min(AR_all);
AR_max_global = max(AR_all);

% X ticks (use tested area ratios, but expressed in %)
xTicks      = AR_all * 100;
xTickLabels = arrayfun(@(v) sprintf('%.0f', v), xTicks, 'UniformOutput', false);

% Y ticks (smooth over global Ld range)
yTicks      = linspace(Ld_min_global, Ld_max_global, nTicks_Ld);
yTickLabels = arrayfun(@(v) sprintf('%.1f', v), yTicks, 'UniformOutput', false);

%% --- Figure setup: 3 x 4, full A4 width in landscape ---
fig = figure('Units','centimeters', ...
             'Position',[1 1 29.7 18], ...
             'PaperUnits','centimeters', ...
             'PaperSize',[29.7 21], ...
             'PaperPosition',[0 0 29.7 18]);

for iSoil = 1:nSoils

    % ---------- Load database for this soil ----------
    S = load(db_files{iSoil}, 'Metrics');
    Metrics = S.Metrics;

    % ---------- Extract design variables (local, but ranges are global) ----------
    Ld_vals = sort(unique([Metrics.Ld_m]).');          % [m]
    AR_vals = sort(unique([Metrics.area_ratio]).');    % [-]

    nLd = numel(Ld_vals);
    nAR = numel(AR_vals);

    % ---------- Build metric surfaces ----------
    Z_eta_p   = nan(nLd, nAR);   % peak flow reduction [%]
    Z_dt_peak = nan(nLd, nAR);   % time-to-peak delay [min]
    Z_det_hr  = nan(nLd, nAR);   % detention time [h]
    Z_hmax    = nan(nLd, nAR);   % maximum ponding depth [m]

    for iLd = 1:nLd
        for jAR = 1:nAR
            idx = find( abs([Metrics.Ld_m]       - Ld_vals(iLd))  < 1e-6 & ...
                        abs([Metrics.area_ratio] - AR_vals(jAR)) < 1e-6, ...
                        1, 'first');
            if ~isempty(idx)
                Z_eta_p(iLd, jAR)   = Metrics(idx).eta_p_pct;        % [%]
                Z_dt_peak(iLd, jAR) = Metrics(idx).Delta_tp_min;     % [min]
                Z_det_hr(iLd, jAR)  = Metrics(idx).DetTime_min/60;   % [h]
                Z_hmax(iLd, jAR)    = Metrics(idx).max_ponding;      % [m]
            end
        end
    end

    % Ensure strictly positive values for log-scale time-to-peak delay
    Z_dt_peak(Z_dt_peak <= 0) = NaN;

    % ---------- Meshgrid ----------
    [AR_grid, Ld_grid] = meshgrid(AR_vals, Ld_vals);
    AR_grid_pct = AR_grid * 100;

    %% Column 1: Peak flow reduction
    ax = subplot(nSoils, 4, (iSoil-1)*4 + 1);

    surf(AR_grid_pct, Ld_grid, Z_eta_p);
    view(2); shading interp;
    xlabel('Area ratio [\%]');
    ylabel('$L_d\;[\mathrm{m}]$');
    xlim([AR_min_global AR_max_global] * 100);
    ylim([Ld_min_global Ld_max_global]);

    if iSoil == 1
        title('Peak flow reduction');
    end

    set(ax, 'TickDir','in', ...
            'Box','on', ...
            'Layer','top', ...
            'FontSize',14, ...
            'LineWidth',2.0, ...
            'XTick', xTicks, ...
            'XTickLabel', xTickLabels, ...
            'YTick', yTicks, ...
            'YTickLabel', yTickLabels);
    grid(ax, 'off');

    caxis(ax, eta_clim);
    colormap(ax, parula(nColors_eta));

    cb = colorbar(ax);
    cb.Label.Interpreter = 'latex';
    cb.Label.String      = '$\eta_p\;[\%]$';
    cb.LineWidth         = 2.0;
    cb.TickDirection     = 'in';
    cb.TickLength        = 0.025;
    cb.FontSize          = 14;

    eta_ticks = round(linspace(eta_clim(1), eta_clim(2), nColors_eta+1));
    cb.Ticks  = unique(eta_ticks);

    %% Column 2: Time-to-peak delay (log-scale color)
    ax = subplot(nSoils, 4, (iSoil-1)*4 + 2);

    surf(AR_grid_pct, Ld_grid, Z_dt_peak);
    view(2); shading interp;
    xlabel('Area ratio [\%]');
    ylabel('$L_d\;[\mathrm{m}]$');
    xlim([AR_min_global AR_max_global] * 100);
    ylim([Ld_min_global Ld_max_global]);

    if iSoil == 1
        title('Time-to-peak delay');
    end

    set(ax, 'TickDir','in', ...
            'Box','on', ...
            'Layer','top', ...
            'FontSize',14, ...
            'LineWidth',2.0, ...
            'XTick', xTicks, ...
            'XTickLabel', xTickLabels, ...
            'YTick', yTicks, ...
            'YTickLabel', yTickLabels, ...
            'ColorScale','log');   % log-scale for color mapping
    grid(ax, 'off');

    caxis(ax, dt_clim);
    % colormap(ax, spring(nColors_dt));
    colormap(ax,spring(6))

    cb = colorbar(ax);
    cb.Label.Interpreter = 'latex';
    cb.Label.String      = '$\Delta t_p\;[\mathrm{min}]$';
    cb.LineWidth         = 2.0;
    cb.TickDirection     = 'in';
    cb.TickLength        = 0.025;
    cb.FontSize          = 14;

    % Log-spaced ticks within dt_clim
    dt_ticks_all = [1 3 10 30 100 300 900];
    dt_ticks = dt_ticks_all(dt_ticks_all >= dt_clim(1) & dt_ticks_all <= dt_clim(2));
    cb.Ticks = dt_ticks;

    %% Column 3: Detention time (spring colormap)
    ax = subplot(nSoils, 4, (iSoil-1)*4 + 3);

    surf(AR_grid_pct, Ld_grid, Z_det_hr);
    view(2); shading interp;
    xlabel('Area ratio [\%]');
    ylabel('$L_d\;[\mathrm{m}]$');
    xlim([AR_min_global AR_max_global] * 100);
    ylim([Ld_min_global Ld_max_global]);

    if iSoil == 1
        title('Detention time');
    end

    set(ax, 'TickDir','in', ...
            'Box','on', ...
            'Layer','top', ...
            'FontSize',14, ...
            'LineWidth',2.0, ...
            'XTick', xTicks, ...
            'XTickLabel', xTickLabels, ...
            'YTick', yTicks, ...
            'YTickLabel', yTickLabels);
    grid(ax, 'off');

    caxis(ax, det_clim);
    colormap(ax, autumn(nColors_det));

    cb = colorbar(ax);
    cb.Label.Interpreter = 'latex';
    cb.Label.String      = '$t_d\;[\mathrm{h}]$';
    cb.LineWidth         = 2.0;
    cb.TickDirection     = 'in';
    cb.TickLength        = 0.025;
    cb.FontSize          = 14;

    det_ticks = round(linspace(det_clim(1), det_clim(2), nColors_det+1));
    cb.Ticks  = unique(det_ticks);

    %% Column 4: Maximum ponding depth (continuous hot scale)
    ax = subplot(nSoils, 4, (iSoil-1)*4 + 4);

    surf(AR_grid_pct, Ld_grid, Z_hmax);
    view(2); shading interp;
    xlabel('Area ratio [\%]');
    ylabel('$L_d\;[\mathrm{m}]$');
    xlim([AR_min_global AR_max_global] * 100);
    ylim([Ld_min_global Ld_max_global]);

    if iSoil == 1
        title('Maximum ponding depth');
    end

    set(ax, 'TickDir','in', ...
            'Box','on', ...
            'Layer','top', ...
            'FontSize',14, ...
            'LineWidth',2.0, ...
            'XTick', xTicks, ...
            'XTickLabel', xTickLabels, ...
            'YTick', yTicks, ...
            'YTickLabel', yTickLabels);
    grid(ax, 'off');

    caxis(ax, hmax_clim);
    % Continuous hot color scale
    colormap(ax, hot);

    cb = colorbar(ax);
    cb.Label.Interpreter = 'latex';
    cb.Label.String      = '$h_{\max}\;[\mathrm{m}]$';
    cb.LineWidth         = 2.0;
    cb.TickDirection     = 'in';
    cb.TickLength        = 0.025;
    cb.FontSize          = 14;

    % "Nice" ticks for ponding depth
    cb.Ticks = linspace(hmax_clim(1), hmax_clim(2), nColors_hmax+1);

end

drawnow;
