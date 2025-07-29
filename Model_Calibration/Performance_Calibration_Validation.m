%% Run Calibration vs Validation Overview
filename_cal = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Model_Calibration\Rainfall_ETP_Discharge_Obs_Discharge_Timeseries_Events_new.xlsx';
sheetname_cal = 'Sheet1';

filename_val = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Model_Calibration\Rainfall_ETP_Discharge_Obs_Discharge_Timeseries_new.xlsx';
sheetname_val = 'Sheet1';

% Event definitions
event_calib(1).start = datetime(2024,1,21,22,0,0);
event_calib(1).end   = datetime(2024,1,23);
event_calib(2).start = datetime(2024,1,2,12,0,0);
event_calib(2).end   = datetime(2024,1,3,12,0,0);

event_valid(1).start = datetime(2024,4,28,6,0,0);
event_valid(1).end   = datetime(2024,4,28,16,0,0);
event_valid(2).start = datetime(2024,5,31);
event_valid(2).end   = datetime(2024,6,1);

close all
analyze_discharge_performance(filename_cal, sheetname_cal, ...
    filename_val, sheetname_val, ...
    event_calib, event_valid);



function analyze_discharge_performance(filename_cal, sheetname_cal, ...
    filename_val, sheetname_val, ...
    event_calib, event_valid)
% filename: Excel file with timeseries data
% sheet_calib: sheet name for calibration
% sheet_valid: sheet name for validation
% event_calib, event_valid: two-element struct arrays with fields "start" and "end" (datetime)

% Load both datasets
data_cal = readtable(filename_cal, 'Sheet', sheetname_cal);
data_val = readtable(filename_val, 'Sheet', sheetname_val);


% Convert datetime and extract variables
[time_cal, rain_cal, etp_cal, Qobs_cal, Qsim_cal] = extract_data(data_cal);
[time_val, rain_val, etp_val, Qobs_val, Qsim_val] = extract_data(data_val);

dt_hr = 15/60;

% Performance metrics
metrics_cal = compute_metrics(Qobs_cal, Qsim_cal, rain_cal, etp_cal, dt_hr);
metrics_val = compute_metrics(Qobs_val, Qsim_val, rain_val, etp_val, dt_hr);

% Create Figure
figure('Color','w','Position',[100 100 1600 1200]);
tiledlayout(3,4, 'Padding', 'compact', 'TileSpacing','compact');

%% Row 1: Calibration
nexttile(1); plot_hydrograph(time_cal, rain_cal, Qobs_cal, Qsim_cal, 'Calibration');
format_axes()
nexttile(2); plot_event(time_cal, rain_cal, Qobs_cal, Qsim_cal, event_calib(1), 'Cal Event 1');
format_axes()
nexttile(3); plot_event(time_cal, rain_cal, Qobs_cal, Qsim_cal, event_calib(2), 'Cal Event 2');
format_axes()
nexttile(4); plot_scatter(Qobs_cal, Qsim_cal, metrics_cal, 'Calibration');

%% Row 2: Validation
nexttile(5); plot_hydrograph(time_val, rain_val, Qobs_val, Qsim_val, 'Validation');
nexttile(6); plot_event(time_val, rain_val, Qobs_val, Qsim_val, event_valid(1), 'Val Event 1');
nexttile(7); plot_event(time_val, rain_val, Qobs_val, Qsim_val, event_valid(2), 'Val Event 2');
nexttile(8); plot_scatter(Qobs_val, Qsim_val, metrics_val, 'Validation');

%% Row 3: FDCs and Cumulative
nexttile(9); plot_fdc(Qobs_cal, Qsim_cal, 'FDC Calibration');
nexttile(10); plot_fdc(Qobs_val, Qsim_val, 'FDC Validation');
nexttile(11); plot_cumulative(time_cal, rain_cal, etp_cal, Qobs_cal, Qsim_cal, 'Cumulative (Calibration)');
nexttile(12); plot_cumulative(time_val, rain_val, etp_val, Qobs_val, Qsim_val, 'Cumulative (Validation)');

letters = 'abcdefghijkl';
for i = 1:12
    ax = nexttile(i);
    text(ax.XLim(1), ax.YLim(2), ['(', letters(i), ')'], ...
        'FontWeight', 'bold', 'FontSize', 16, ...
        'FontName', 'Montserrat', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');
end


sgtitle('Calibration and Validation Overview', 'FontWeight', 'bold', 'FontSize', 16);

% At the end of the function
exportgraphics(gcf, 'Calibration_Validation_Overview.png', 'Resolution', 300);

end

%% === Helper Functions ===

function format_axes()
ax = gca;

% Standard formatting
set(ax, ...
    'FontName','Montserrat', ...
    'FontSize',14, ...
    'LineWidth',2.5, ...
    'TickDir','in', ...
    'TickLength',[0.02 0.01], ...
    'GridLineStyle','-', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

% Set black axis color without forcing creation of right axis
yyaxis left
ax.YColor = 'k';

% If this plot already has a second y-axis
if strcmp(ax.YAxis(2).Visible, 'on')
    yyaxis right
    ax.YColor = 'k';
end

grid on;
end

function format_axes_scatter()
ax = gca;

% Standard formatting
set(ax, ...
    'FontName','Montserrat', ...
    'FontSize',14, ...
    'LineWidth',2.5, ...
    'TickDir','in', ...
    'TickLength',[0.02 0.01], ...
    'GridLineStyle','-', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

grid on;
end


function tf = hasTwoYAxes(ax)
% Check if the axes has yyaxis enabled
yy = findall(ax, 'Type', 'line');
yAxes = get(yy, 'YData');
if iscell(yAxes)
    tf = any(cellfun(@(y) any(~isnan(y)), yAxes));
else
    tf = any(~isnan(yAxes));
end
end

function [t, r, e, qobs, qsim] = extract_data(data)
t = datetime(data{:,1});
r = data{:,2}; e = data{:,3};
qobs = data{:,4}; qsim = data{:,5};
end

function metrics = compute_metrics(qobs, qsim, rain, etp, dt)
dt_hr = 15/60;
dt_break = round(77.75/dt_hr); % 1 day
qobs(1:dt_break) = [];
qsim(1:dt_break) = [];
rain(1:dt_break) = [];
etp(1:dt_break) = [];
NSE = 1 - sum((qobs - qsim).^2) / sum((qobs - mean(qobs)).^2);
RMSE = sqrt(mean((qobs - qsim).^2));
r = corr(qobs, qsim, 'rows','complete');
alpha = std(qsim) / std(qobs);
beta = mean(qsim) / mean(qobs);
KGE = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
KGE_lim = 1 - sqrt((r - 1)^2 + (beta - 1)^2);
PBIAS = 100 * (sum(qsim - qobs) / sum(qobs));
rain_total = sum(rain) * dt;
etp_total = sum(etp) * dt;
Qobs_total = sum(qobs) * dt;
Qsim_total = sum(qsim) * dt;
RC_obs = Qobs_total / rain_total;
RC_sim = Qsim_total / rain_total;
metrics = struct('NSE', NSE, 'RMSE', RMSE, 'KGE', KGE, 'KGE_lim', KGE_lim, ...
    'PBIAS', PBIAS, 'R2', r^2, 'RC_obs', RC_obs, 'RC_sim', RC_sim, ...
    'rain_total', rain_total, 'etp_total', etp_total);
end

function plot_hydrograph(t, r, qobs, qsim, label)
bright = get_bright_colors();  % use helper function
yyaxis left
plot(t, qobs, 'o', 'Color', 'k', 'MarkerSize', 6, 'DisplayName', 'Observed'); hold on
plot(t, qsim, '-', 'Color', bright.red, 'LineWidth', 2.5, 'DisplayName', 'Simulated');
ylabel('Discharge [mm/h]')

yyaxis right
bar(t, r, 1, 'FaceColor', bright.dark_blue, 'EdgeColor', bright.dark_blue, 'BarWidth', 1);
ylabel('Rainfall [mm/h]');
set(gca, 'YDir', 'reverse');

title(['Hydrograph - ', label])
format_axes();
ylim([0,150])
end

function plot_event(t, r, qobs, qsim, ev, title_str)
bright = get_bright_colors();  % Use same color source as hydrograph

% Extract event time window
idx = t >= ev.start & t <= ev.end;
t_ev = t(idx); r_ev = r(idx);
qobs_ev = qobs(idx); qsim_ev = qsim(idx);

% Compute NSE (only if data exists)
if all(isnan(qobs_ev)) || all(isnan(qsim_ev))
    NSE = NaN;
else
    NSE = 1 - sum((qobs_ev - qsim_ev).^2) / sum((qobs_ev - mean(qobs_ev)).^2);
end

% Left Y-axis: Discharge
yyaxis left
plot(t_ev, qsim_ev, '-', 'Color', bright.red, 'LineWidth', 2.5, 'DisplayName', 'Simulated'); hold on    
plot(t_ev, qobs_ev, '.', 'Color', 'k', 'MarkerSize', 16, 'DisplayName', 'Observed'); 
ylabel('Discharge [mm/h]');

% Right Y-axis: Rainfall
yyaxis right
bar(t_ev, r_ev, 1, 'FaceColor', bright.dark_blue, ...
    'EdgeColor', bright.dark_blue, 'BarWidth', 1);

ylabel('Rainfall [mm/h]');
set(gca, 'YDir', 'reverse');

% Title
title(title_str);

% Format axes consistently
format_axes();

% Annotate NSE value (top-right corner of left axis)
yyaxis left
ylims = ylim;
xlims = xlim;
x_pos = xlims(1) + 0.75 * (xlims(2) - xlims(1));
y_pos = ylims(1) + 0.9 * (ylims(2) - ylims(1));
text(x_pos, y_pos, sprintf('NSE = %.3f', NSE), ...
    'FontSize', 10, ...
    'FontName', 'Montserrat', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'Margin', 4);
axis tight
yyaxis right
ylim([0,150])
end




function plot_scatter(qobs, qsim, metrics, label)
bright = get_bright_colors();  % Returns struct with hex strings
color_used = hex2rgb(bright.purple);  % Convert hex to RGB

% Plot scatter
scatter(qobs, qsim, 25, color_used, 'filled', 'DisplayName', 'Data Points'); hold on;
line([0 max(qobs)], [0 max(qobs)], 'Color', 'k', 'LineStyle', '--','linewidth',2);
axis equal; grid on;
xlabel('Observed [mm/h]'); ylabel('Simulated [mm/h]');
title(['Scatter - ', label]);
txt = sprintf(['\\bfMetrics\\rm\nNSE=%.3f\nRMSE=%.2f\nKGE=%.3f\nKGE_{lim}=%.3f\nPBIAS=%.1f%%\nR^2=%.3f\nRC_{obs}=%.2f\nRC_{sim}=%.2f'], ...
    metrics.NSE, metrics.RMSE, metrics.KGE, metrics.KGE_lim, ...
    metrics.PBIAS, metrics.R2, metrics.RC_obs, metrics.RC_sim);
xlim([0 max(qobs)]); ylim([0 max(qsim)]);
xpos = max(qobs)*0.05; ypos = max(qsim)*0.95;
text(xpos, ypos, txt, 'FontSize', 9, 'BackgroundColor','white','EdgeColor',[0.4 0.4 0.4],'Margin',5);
format_axes_scatter()
end



function plot_fdc(qobs, qsim, title_str)
    bright = get_bright_colors();  % Load bright palette

    % Threshold: remove values ≤ 0.001
    threshold = 0.001;
    valid_idx = (qobs > threshold) & (qsim > threshold);
    qobs_valid = qobs(valid_idx);
    qsim_valid = qsim(valid_idx);

    % Sort descending for FDC
    qobs_sorted = sort(qobs_valid, 'descend');
    qsim_sorted = sort(qsim_valid, 'descend');
    p_obs = (1:length(qobs_sorted)) / length(qobs_sorted) * 100;
    p_sim = (1:length(qsim_sorted)) / length(qsim_sorted) * 100;

    % Trim to same length (for NSE)
    min_len = min(length(qobs_sorted), length(qsim_sorted));
    qobs_sorted = qobs_sorted(1:min_len);
    qsim_sorted = qsim_sorted(1:min_len);
    p_obs = p_obs(1:min_len);
    p_sim = p_sim(1:min_len);

    % Plot log-scale FDC
    semilogy(p_obs, qobs_sorted, '-', ...
        'Color',  bright.med_blue, 'LineWidth', 2.5, 'DisplayName', 'Observed'); hold on;
    semilogy(p_sim, qsim_sorted, '--', ...
        'Color', bright.green, 'LineWidth', 2.5, 'DisplayName', 'Simulated');

    xlabel('Exceedance Probability [%]');
    ylabel('Discharge [mm/h]');
    title(title_str);
    legend('Location', 'northeast');

    % Custom Y ticks and labels
    ytick_vals = [0.001, 0.01, 0.1, 1, 10];
    yticks(ytick_vals);
    yticklabels({'0.001','0.01','0.1','1','10'});
    ylim([0.001, 10]);
    grid off

    % --- Compute NSE between FDC curves ---
    NSE_fdc = 1 - sum((qobs_sorted - qsim_sorted).^2) / sum((qobs_sorted - mean(qobs_sorted)).^2);

    % Annotate NSE
    xlims = xlim; ylims = ylim;
    x_pos = xlims(1) + 0.6 * (xlims(2) - xlims(1));
    y_pos = 10^(log10(ylims(1)) + 0.8 * (log10(ylims(2)) - log10(ylims(1))));
    text(x_pos, y_pos, sprintf('NSE = %.3f', NSE_fdc), ...
        'FontSize', 10, ...
        'FontName', 'Montserrat', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', [0.4 0.4 0.4], ...
        'Margin', 5);

    format_axes_scatter();
end


function plot_cumulative(t, r, e, qobs, qsim, title_str)
    dt = 15/60;

    % Compute cumulative values
    cum_r     = cumsum(r)     * dt;
    cum_e     = cumsum(e)     * dt;
    cum_qobs  = cumsum(qobs)  * dt;
    cum_qsim  = cumsum(qsim)  * dt;

    % Load qualitative color palette
    qual = get_qual_colors();

    % Plot using defined colors and styles
    plot(t, cum_r, '-',  'Color', qual.dark_blue, 'DisplayName','Rain',  'LineWidth',2); hold on;
    plot(t, cum_e, '-',  'Color', qual.red,       'DisplayName','ETP',   'LineWidth',2);
    plot(t, cum_qobs, '-', 'Color', qual.teal,    'DisplayName','Q obs', 'LineWidth',2);
    plot(t, cum_qsim, '--', 'Color', qual.purple, 'DisplayName','Q sim', 'LineWidth',2);

    ylabel('Cumulative Depth [mm]');
    title(title_str);
    legend show;
    
    format_axes_scatter();  % consistent styling
end


function c = get_bright_colors()
c = struct( ...
    'dark_blue', '#003a7d', ...
    'med_blue', '#008dff', ...
    'pink', '#ff73b6', ...
    'purple', '#c701ff', ...
    'green', '#4ecb8d', ...
    'orange', '#ffd93a', ...
    'yellow', '#f9e858', ...
    'red', '#d83034');
end

function c = get_qual_colors()
c = struct( ...
    'gray', '#cecece', ...
    'purple', '#a559aa', ...
    'teal', '#59a89c', ...
    'gold', '#f0c571', ...
    'red', '#e02b35', ...
    'dark_blue', '#082a54');
end
