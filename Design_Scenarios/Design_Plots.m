%% DRAIN-LID design maps: ATC vs Ld (colored performance maps)
% Rows  = Peff
% Cols  = metrics (eta_p, DetTime, Delta_tp)
% Color = metric value (averaged over all tp values)
%
% Each panel shows a colored map of the metric in ATC–Ld space for a given Peff.
% For each (metric, Peff) panel, values are averaged across all available tp.

clear; clc;
load('DRAIN_LID_database.mat', 'Metrics');

%% Global style settings
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

metricFields = {'eta_p_pct', 'DetTime_min', 'Delta_tp_min'};
metricTitles = { ...
    'Peak flow reduction, $\eta_{\mathrm{p}}$ [\%]', ...
    'Detention time [min]', ...
    'Time-to-peak delay, $\Delta t_{\mathrm{p}}$ [min]'};
nMetrics = numel(metricFields);

%% Unique soils
soilNames = unique({Metrics.soilName}, 'stable');

for iSoil = 1:numel(soilNames)

    soilName = soilNames{iSoil};
    idxSoil  = strcmp({Metrics.soilName}, soilName);
    Msoil    = Metrics(idxSoil);

    % Unique design / forcing values for this soil
    ATC_vals  = unique([Msoil.ATC_m2]);
    Ld_vals   = unique([Msoil.Ld_m]);
    Peff_vals = unique([Msoil.Peff_catch_mm]);
    Tp_vals   = unique([Msoil.tp_min]);   % [min]

    nA   = numel(ATC_vals);
    nLd  = numel(Ld_vals);
    nP   = numel(Peff_vals);
    nTp  = numel(Tp_vals);

    % Precompute metric grids: (nLd x nA x nP x nTp)
    metricGrids = cell(nMetrics, 1);
    for im = 1:nMetrics
        metricGrids{im} = nan(nLd, nA, nP, nTp);
    end

    % Fill metricGrids
    for ip = 1:nP
        Peff_now = Peff_vals(ip);
        for it = 1:nTp
            Tp_now = Tp_vals(it);

            idxPT = ([Msoil.Peff_catch_mm] == Peff_now) & ...
                    ([Msoil.tp_min]        == Tp_now);

            MsoilPT = Msoil(idxPT);
            if isempty(MsoilPT), continue; end

            for iA = 1:nA
                ATC_now = ATC_vals(iA);
                idxA    = [MsoilPT.ATC_m2] == ATC_now;
                MsoilPTA = MsoilPT(idxA);
                if isempty(MsoilPTA), continue; end

                for iL = 1:nLd
                    Ld_now = Ld_vals(iL);
                    idxL   = [MsoilPTA.Ld_m] == Ld_now;
                    Mcell  = MsoilPTA(idxL);
                    if isempty(Mcell), continue; end

                    for im = 1:nMetrics
                        vals = [Mcell.(metricFields{im})];
                        metricGrids{im}(iL, iA, ip, it) = mean(vals, 'omitnan');
                    end
                end
            end
        end
    end

    % --------- Figure for this soil: rows = Peff, cols = metrics ----------
    fig = figure('Units', 'inches', ...
                 'Position', [1 1 8.27 11.69], ... % A4
                 'PaperUnits', 'inches', ...
                 'PaperSize', [8.27 11.69], ...
                 'PaperPosition', [0 0 8.27 11.69], ...
                 'Color', 'w');

    tlay = tiledlayout(fig, nP, nMetrics, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    sgtitle(tlay, sprintf('%s -- Metrics in $A_{\\mathrm{TC}}$--$L_{\\mathrm{d}}$ space', soilName), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % ------------------------ Plotting loop -------------------------------
    for ip = 1:nP          % rows: Peff
        for im = 1:nMetrics % cols: metric

            nexttile(tlay, (ip-1)*nMetrics + im);
            hold on; box on;

            Zall = metricGrids{im};     % [nLd x nA x nP x nTp] for this metric

            % Extract all tp surfaces for this Peff
            Z_Peff_allTp = squeeze(Zall(:,:,ip,:));  % [nLd x nA x nTp]

            % Aggregate over tp: here we use mean over tp (omit NaNs)
            Z_mean = mean(Z_Peff_allTp, 3, 'omitnan');  % [nLd x nA]

            if all(isnan(Z_mean), 'all')
                % Nothing simulated for this panel
                text(0.5,0.5,'No data', ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'FontSize',10);
                axis off;
                continue;
            end

            % Local min/max for this panel (for color scaling)
            Zflat = Z_mean(:);
            Zflat = Zflat(~isnan(Zflat));
            zMin = min(Zflat);
            zMax = max(Zflat);

            % Plot colored map in ATC–Ld space
            % Note: imagesc needs 'axis xy' to keep y increasing upwards
            hImg = imagesc(ATC_vals, Ld_vals, Z_mean);
            set(hImg, 'AlphaData', ~isnan(Z_mean)); % make NaNs transparent
            set(gca, 'YDir', 'normal');             % same as axis xy
            colormap(gca, parula);                  % or any colormap you like

            % Color limits
            if zMin == zMax
                zMax = zMin + 1;
            end
            caxis([zMin zMax]);

            % Colorbar (one per panel; if too busy, you can move this outside)
            cb = colorbar;
            cb.TickLabelInterpreter = 'latex';
            cb.FontSize = 10;

            % Axes styling
            set(gca, 'TickDir', 'out', ...
                     'FontSize', 12, ...
                     'LineWidth', 1.0, ...
                     'Layer', 'top');

            xlim([min(ATC_vals), max(ATC_vals)]);
            ylim([min(Ld_vals),  max(Ld_vals)]);

            % Column titles (metrics)
            if ip == 1
                title(metricTitles{im}, 'FontSize', 14);
            end

            % X label only on bottom row
            if ip == nP
                xlabel('$A_{\mathrm{TC}}~[\mathrm{m}^2]$', 'FontSize', 14);
            end

            % Y label only on first column
            if im == 1
                ylabel('$L_{\mathrm{d}}~[\mathrm{m}]$', 'FontSize', 14);
            end

            % Row label: Peff (top-left inside axes)
            xl = xlim; yl = ylim;
            txt = sprintf('$P_{\\mathrm{eff}} = %g~\\mathrm{mm}$', Peff_vals(ip));
            text(xl(1) + 0.02*(xl(2)-xl(1)), ...
                 yl(2) - 0.05*(yl(2)-yl(1)), ...
                 txt, ...
                 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'top', ...
                 'FontSize', 11);

            hold off;
        end
    end

    % Optional: save figure per soil
    % fname = sprintf('DesignMaps_%s.pdf', strrep(soilName,' ','_'));
    % print(fig, fname, '-dpdf', '-painters');

end
