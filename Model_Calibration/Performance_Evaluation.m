%% Run Hydrograph Analysis
filename = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Model_Calibration\Rainfall_ETP_Discharge_Obs_Discharge_Timeseries_Events_new.xlsx';
sheetname = 'Sheet1';
close all
analyze_discharge_performance(filename, sheetname)

function analyze_discharge_performance(filename, sheetname)
% Analyze hydrologic performance between observed and modeled discharge (mm/h)
%
% Inputs:
%   filename   – Excel file with timeseries (datetime, rain, etp, Qobs, Qsim)
%   sheetname  – Sheet name in Excel file

%% === Load Data ===
data = readtable(filename, 'Sheet', sheetname);
time        = datetime(data{:,1});
rain_mm_h   = data{:,2};
etp_mm_h    = data{:,3};
Q_obs_mm_h  = data{:,4};
Q_sim_mm_h  = data{:,5};


dt_hr = 15 / 60;  % h

%% === Break Data === %
dt_break = round(77.75/dt_hr); % 1 day
time(1:dt_break) = [];
rain_mm_h(1:dt_break) = [];
etp_mm_h(1:dt_break) = [];
Q_obs_mm_h(1:dt_break) = [];
Q_sim_mm_h(1:dt_break) = [];

% === Colors and markers ===
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();
colors = [pallete.blue_colors(1,:); pallete.red_colors(1,:); [192,192,192]/255];

% Font
set(groot, 'defaultAxesFontName', 'Montserrat');
set(groot, 'defaultTextFontName', 'Montserrat');
%% === Performance Metrics (in mm/h) ===
NSE     = 1 - sum((Q_obs_mm_h - Q_sim_mm_h).^2) / sum((Q_obs_mm_h - mean(Q_obs_mm_h)).^2);
NSElog  = 1 - sum((log(Q_obs_mm_h + eps) - log(Q_sim_mm_h + eps)).^2) / ...
              sum((log(Q_obs_mm_h + eps) - mean(log(Q_obs_mm_h + eps))).^2);
RMSE    = sqrt(mean((Q_obs_mm_h - Q_sim_mm_h).^2));
PBIAS   = 100 * (sum(Q_sim_mm_h - Q_obs_mm_h) / sum(Q_obs_mm_h));
r       = corr(Q_sim_mm_h, Q_obs_mm_h, 'rows','complete');
alpha   = std(Q_sim_mm_h) / std(Q_obs_mm_h);
beta    = mean(Q_sim_mm_h) / mean(Q_obs_mm_h);
KGE     = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
KGE_lim = 1 - sqrt((r - 1)^2 + (beta - 1)^2);

%% === Runoff Coefficients (unitless, all in mm) ===
rain_total = sum(rain_mm_h) * dt_hr;
Qobs_total = sum(Q_obs_mm_h) * dt_hr;
Qsim_total = sum(Q_sim_mm_h) * dt_hr;
etp_total  = sum(etp_mm_h) * dt_hr;

RC_obs = Qobs_total / rain_total;
RC_sim = Qsim_total / rain_total;

%% === Metric Text Box ===
% === Build plain TeX string with line breaks ===
metrics_text = sprintf(['\\bfPerformance Metrics\\rm\n' ...
    'NSE = %.3f\n' ...
    'RMSE = %.3f mm/h\n' ...
    'KGE = %.3f\n' ...
    'KGE_{lim} = %.3f\n' ...
    'PBIAS = %.2f%%%%\n' ...
    'R^2 = %.3f\n' ...
    'RC_{obs} = %.2f\n' ...
    'RC_{sim} = %.2f\n' ...
    'Total Rain = %.1f mm\n' ...
    'Total ETP = %.1f mm'], ...
    NSE, RMSE, KGE, KGE_lim, PBIAS, r^2, RC_obs, RC_sim, rain_total, etp_total);

% === Place the text inside the plot using text() ===
x_pos = time(round(length(time) * 0.75));  % 75% along the time axis
y_pos = max(Q_obs_mm_h) * 0.9;             % 90% of the max discharge

text(x_pos, y_pos, metrics_text, ...
    'Interpreter', 'tex', ...
    'FontSize', 11, ...
    'FontName', 'Montserrat', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'Margin', 5);


figure('Color','w','Position',[100 100 1200 500]);  % Wide figure for two subplots

%% === Subplot 1: Hydrograph and Rainfall ===
% === TOGGLE: Use full time series or a specific interval ===
use_full_series = true;  % Set to true to use full time series

if use_full_series
    t_start = min(time);
    t_end   = max(time);
else
    t_start = datetime(2024, 1, 21, 20, 0, 0);
    t_end   = datetime(2024, 1, 23, 0, 0, 0);
end

% === Filter time window ===
interval_idx = time >= t_start & time <= t_end;

% Use filtered series for NSE and plotting
Q_obs_sel = Q_obs_mm_h(interval_idx);
Q_sim_sel = Q_sim_mm_h(interval_idx);
time_sel  = time(interval_idx);
rain_sel  = rain_mm_h(interval_idx);

% === Compute NSE for interval ===
NSE_interval = 1 - sum((Q_obs_sel - Q_sim_sel).^2) / sum((Q_obs_sel - mean(Q_obs_sel)).^2);

% === Subplot 1: Hydrograph with Rainfall ===
subplot(1,2,1)

yyaxis left
set(gca, 'YColor', 'black');

% Observed flow: red circles only
plot(time, Q_obs_mm_h, 'o', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', ...
    'LineStyle', 'none', ...
    'MarkerSize', 6, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Observed Q'); hold on;

% Simulated flow: thick continuous line
plot(time, Q_sim_mm_h, '-', ...
    'Color', colors(2,:), ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Simulated Q');
ylabel('Discharge [mm/h]', 'FontName', 'Montserrat', 'FontSize', 12);
set(gca,'ycolor','black')

yyaxis right
area(time, rain_mm_h, ...
    'FaceColor', colors(1,:), ...
    'FaceAlpha', 0.5, ...
    'DisplayName', 'Rainfall');
set(gca, 'YColor', 'black')
ylabel('Net Rainfall [mm/h]','FontName','Montserrat');
ylim([0, max(rain_mm_h)*4]);
set(gca, 'YDir', 'reverse');
xlim([t_start, t_end]);


xlabel('Time', 'FontName', 'Montserrat', 'FontSize', 12);
title('Observed vs Simulated Discharge with Rainfall', ...
    'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Montserrat');

legend('Location', 'northwest', 'Box', 'on', 'FontName', 'Montserrat', 'FontSize', 11);
grid on;

set(gca, 'TickDir', 'in', 'Box', 'on', ...
    'FontName', 'Montserrat', 'FontSize', 12, ...
    'LineWidth', 2.5, 'XMinorTick', 'on', 'YMinorTick', 'on');

% === Plot NSE as text inside chart ===
x_nse = t_start + minutes(30);
y_nse = max(Q_obs_mm_h(interval_idx)) * 0.9;

text(x_nse, y_nse, ...
    sprintf('NSE_{interval} = %.3f', NSE_interval), ...
    'FontSize', 11, ...
    'FontName', 'Montserrat', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'Margin', 5);

%% === Subplot 2: Scatter Plot ===
subplot(1,2,2)

% Scatter
scatter(Q_obs_mm_h, Q_sim_mm_h, 25, 'filled', ...
    'MarkerFaceColor', [0.7 0.7 0.7], ...
    'DisplayName', 'Data Points'); hold on;

% 45° reference line
max_val = max([Q_obs_mm_h; Q_sim_mm_h]);
line([0 max_val], [0 max_val], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, ...
    'DisplayName', '1:1 Line');

axis equal
xlim([0, max_val])
ylim([0, max_val])
xlabel('Observed Discharge [mm/h]', 'FontName', 'Montserrat', 'FontSize', 12)
ylabel('Simulated Discharge [mm/h]', 'FontName', 'Montserrat', 'FontSize', 12)
title('Observed vs Simulated Discharge (Scatter)', ...
    'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Montserrat')
grid on

set(gca, 'TickDir', 'in', 'Box', 'on', ...
    'FontName', 'Montserrat', 'FontSize', 12, ...
    'LineWidth', 2.5, 'XMinorTick', 'on', 'YMinorTick', 'on');

legend('Location', 'northwest', 'Box', 'on', 'FontSize', 11)

% === Metrics box using text (subplot-relative) ===
x_pos = 0.05 * max_val;
y_pos = 0.95 * max_val;
text(x_pos, y_pos, metrics_text, ...
    'Interpreter', 'tex', ...
    'FontSize', 11, ...
    'FontName', 'Montserrat', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'Margin', 5);



%% === Plot 2: Cumulative Volumes ===
cum_rain = cumsum(rain_mm_h) * dt_hr;
cum_etp  = cumsum(etp_mm_h)  * dt_hr;
cum_Qobs = cumsum(Q_obs_mm_h) * dt_hr;
cum_Qsim = cumsum(Q_sim_mm_h) * dt_hr;

figure('Color','w','Position',[100 600 1000 400]); hold on;
plot(time, cum_rain, '-', 'Color', pallete.blue_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Rainfall');
plot(time, cum_etp,  '-', 'Color', pallete.red_colors(1,:), 'LineWidth', 2, 'DisplayName', 'ETP');
plot(time, cum_Qobs, '-','Color', pallete.green_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Observed Q');
plot(time, cum_Qsim, '--','Color', pallete.red_colors(2,:), 'LineWidth', 2, 'DisplayName', 'Simulated Q');

xlabel('Time', 'FontName', 'Montserrat', 'FontSize', 12);
ylabel('Cumulative Depth [mm]');
title('Cumulative Rainfall, ETP, and Discharge');
legend('Location', 'northwest', 'Box', 'on', 'FontName', 'Montserrat', 'FontSize', 11);
grid on;
set(gca, 'TickDir', 'out', 'Box', 'on', 'FontName', 'Montserrat', 'FontSize', 12, ...
    'LineWidth', 2.5, 'XMinorTick', 'on', 'YMinorTick', 'on');

%% === Identify Events Based on Inter-Event Period ===
inter_event_hrs = 6;
inter_event_steps = inter_event_hrs / dt_hr;

event_id = zeros(size(rain_mm_h));
event_counter = 0;
i = 1;
while i <= length(rain_mm_h)
    if rain_mm_h(i) > 0
        % Start new event
        event_counter = event_counter + 1;
        start_idx = i;
        i = i + 1;
        while i <= length(rain_mm_h)
            if rain_mm_h(i) == 0
                remaining = length(rain_mm_h) - i;
                if remaining >= inter_event_steps
                    if all(rain_mm_h(i+1:i+inter_event_steps) == 0)
                        break;
                    end
                else
                    break;
                end
            end
            i = i + 1;
        end
        end_idx = i;
        event_id(start_idx:end_idx) = event_counter;
    end
    i = i + 1;
end

unique_events = unique(event_id(event_id > 0));
n_event = numel(unique_events);

fprintf('📌 Identified %d rainfall events.\n', n_event);

%% === IDENTIFY AND PLOT TOP N EVENTS ===
n_top_events = 4;  % Change as desired
inter_event_steps = round(6 / dt_hr);  % 6-hour inter-event dry period

% Identify event blocks
event_id = zeros(size(rain_mm_h));
event_counter = 0;
i = 1;
while i <= length(rain_mm_h)
    if rain_mm_h(i) > 0
        event_counter = event_counter + 1;
        start_idx = i;
        i = i + 1;
        while i <= length(rain_mm_h)
            if rain_mm_h(i) == 0
                end_i = min(i + inter_event_steps, length(rain_mm_h));
                if all(rain_mm_h(i+1:end_i) == 0)
                    break;
                end
            end
            i = i + 1;
        end
        end_idx = min(i, length(rain_mm_h));
        event_id(start_idx:end_idx) = event_counter;
    end
    i = i + 1;
end

% Rank events by total rainfall volume
unique_events = unique(event_id(event_id > 0));
event_volumes = zeros(numel(unique_events), 1);
for e = 1:numel(unique_events)
    idx = find(event_id == unique_events(e));
    if all(idx <= length(rain_mm_h))
        event_volumes(e) = sum(rain_mm_h(idx)) * dt_hr;
    end
end

[~, sorted_idx] = sort(event_volumes, 'descend');
top_events = unique_events(sorted_idx(1:min(n_top_events, end)));

% Plot subplots
figure('Color','w','Position',[100 100 1200 300*n_top_events]);

for p = 1:length(top_events)
    idx = find(event_id == top_events(p));
    if isempty(idx), continue; end

    t_ev = time(idx);
    r_ev = rain_mm_h(idx);
    q_obs_ev = Q_obs_mm_h(idx);
    q_sim_ev = Q_sim_mm_h(idx);
    r_max = max(r_ev);

    subplot(n_top_events, 1, p);

    yyaxis right
    area(t_ev, r_ev, 'FaceColor', colors(1,:), 'FaceAlpha', 1, 'DisplayName', 'Rainfall');
    set(gca, 'YDir', 'reverse', 'YColor', 'k');
    ylim([0, 3 * r_max]);
    ylabel('Rainfall [mm/h]', 'FontName','Montserrat');

    yyaxis left
    plot(t_ev, q_obs_ev, 'o', 'Color', colors(3,:), 'MarkerSize', 5, 'DisplayName', 'Obs'); hold on;
    plot(t_ev, q_sim_ev, 's', 'Color', colors(2,:), 'MarkerSize', 5, 'DisplayName', 'Sim');
    ylabel('Q [mm/h]', 'FontName','Montserrat');

    title(sprintf('Top Event #%d (%.1f mm Rain)', p, event_volumes(sorted_idx(p))), ...
          'FontName', 'Montserrat', 'FontSize', 13, 'FontWeight', 'bold');
    
    grid on;
    legend('Location','northeast');
    set(gca, ...
        'TickDir', 'in', ...
        'Box', 'on', ...
        'FontName', 'Montserrat', ...
        'FontSize', 11, ...
        'LineWidth', 2.5, ...
        'XColor', 'k', ...
        'YColor', 'k');

end

sgtitle('Top Rainfall-Runoff Events (Hydrographs)', ...
        'FontWeight', 'bold', ...
        'FontSize', 16, ...
        'FontName', 'Montserrat');


%% === Unified Plot with Hydrographs, Scatter, and Cumulative Volumes ===
figure('Color','w','Position',[100 100 1200 800]);

% === Bright Color Palette ===
bright_colors = struct( ...
    'dark_blue', [0, 58, 125]/255, ...
    'med_blue',  [0, 141, 255]/255, ...
    'pink',      [255, 115, 182]/255, ...
    'purple',    [199, 1, 255]/255, ...
    'green',     [78, 203, 141]/255, ...
    'orange',    [255, 157, 58]/255, ...
    'yellow',    [249, 232, 88]/255, ...
    'red',       [216, 48, 52]/255 ...
);

%% --- Subplot (2,1:2): Hydrograph ---
subplot(2,2,[1 2]) 
yyaxis left
plot(time, Q_obs_mm_h, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', ...
    'LineStyle', 'none', 'MarkerSize', 6, 'LineWidth', 2.5, ...
    'DisplayName', 'Observed Q'); hold on;
plot(time, Q_sim_mm_h, '-', 'Color', bright_colors.purple, ...
    'LineWidth', 2.5, 'DisplayName', 'Simulated Q');
ylabel('Discharge [mm/h]', 'FontName', 'Montserrat', 'FontSize', 12);

yyaxis right
area(time, rain_mm_h, 'FaceColor', bright_colors.med_blue, ...
    'FaceAlpha', 0.5, 'DisplayName', 'Rainfall');
ylabel('Rainfall [mm/h]', 'FontName', 'Montserrat', 'FontSize', 12);
set(gca, 'YDir', 'reverse', 'YColor', 'k');

title('Observed vs Simulated Discharge with Rainfall', ...
    'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Montserrat');
xlabel('Time', 'FontName', 'Montserrat');
legend('Location', 'northwest');
grid on;
axis tight

set(gca, 'TickDir', 'out', 'LineWidth', 2.5, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.02 0.01]);

%% --- Subplot (2,1,2): Scatter Plot ---
subplot(2,2,3)
scatter(Q_obs_mm_h, Q_sim_mm_h, 25, 'filled', ...
    'MarkerFaceColor', bright_colors.orange, 'DisplayName', 'Data Points'); hold on;
max_val = max([Q_obs_mm_h; Q_sim_mm_h]);
line([0 max_val], [0 max_val], 'Color', 'k', 'LineStyle', '--', ...
    'LineWidth', 2.5, 'DisplayName', '1:1 Line');
axis equal
xlim([0, max_val]); ylim([0, max_val]);
xlabel('Observed [mm/h]','FontName','Montserrat');
ylabel('Simulated [mm/h]','FontName','Montserrat');
title('Observed vs Simulated Discharge (Scatter)', ...
    'FontWeight','bold','FontSize',14);
grid on;
legend('Location','northwest');
set(gca, 'TickDir','out','LineWidth',2.5, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.01]);
axis tight
% === Metrics box using text (subplot-relative) ===
x_pos = 0.05 * max_val;
y_pos = 0.95 * max_val;
text(x_pos, y_pos, metrics_text, ...
    'Interpreter', 'tex', ...
    'FontSize', 11, ...
    'FontName', 'Montserrat', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'Margin', 5);


%% --- Subplot (2,1,4): Cumulative Volumes ---
subplot(2,2,4)
cum_rain = cumsum(rain_mm_h)*dt_hr;
cum_etp  = cumsum(etp_mm_h)*dt_hr;
cum_Qobs = cumsum(Q_obs_mm_h)*dt_hr;
cum_Qsim = cumsum(Q_sim_mm_h)*dt_hr;
plot(time, cum_rain, '-', 'Color', bright_colors.med_blue, ...
    'LineWidth', 2.5, 'DisplayName', 'Rainfall'); hold on;
plot(time, cum_etp, '-', 'Color', bright_colors.orange, ...
    'LineWidth', 2.5, 'DisplayName', 'ETP');
plot(time, cum_Qobs, '-', 'Color', bright_colors.green, ...
    'LineWidth', 2.5, 'DisplayName', 'Observed Q');
plot(time, cum_Qsim, '--', 'Color', bright_colors.red, ...
    'LineWidth', 2.5, 'DisplayName', 'Simulated Q');
xlabel('Time','FontName','Montserrat');
ylabel('Cumulative Depth [mm]','FontName','Montserrat');
title('Cumulative Rainfall, ETP, and Discharge');
legend('Location','northwest');
grid on;
set(gca, 'TickDir','out','LineWidth',2.5, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.01]);
axis tight


end

