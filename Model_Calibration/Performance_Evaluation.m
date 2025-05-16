%% Run Hydrograph Analysis
filename = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Course\Model_Calibration\Rainfall_ETP_Discharge_Obs_Discharge_Timeseries.xlsx';
sheetname = 'Sheet1';

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

dt_hr = 5 / 60;  % 5-min in hours

% === Colors and markers ===
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();
colors = [pallete.blue_colors(1,:); pallete.red_colors(1,:); [192,192,192]/255];

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
metrics_text = sprintf(['\\bfPerformance Metrics\\rm\n' ...
    'NSE       = %.3f\n' ...
    'NSE$_{log}$ = %.3f\n' ...
    'RMSE      = %.3f mm/h\n' ...
    'KGE       = %.3f\n' ...
    'KGE$_{lim}$ = %.3f\n' ...
    'PBIAS     = %.2f %%\n' ...
    'RC$_{obs}$ = %.2f\n' ...
    'RC$_{sim}$ = %.2f\n' ...
    'Total Rain = %.1f mm\n' ...
    'Total ETP  = %.1f mm'], ...
    NSE, NSElog, RMSE, KGE, KGE_lim, PBIAS, RC_obs, RC_sim, rain_total, etp_total);

%% === Plot 1: Hydrographs and Rainfall ===
figure('Color','w','Position',[100 100 1000 400]);

yyaxis left
set(gca,'Ycolor','black');
plot(time, Q_obs_mm_h, 'Color',colors(3,:), 'MarkerSize', 4, 'DisplayName', 'Observed Q'); hold on;
plot(time, Q_sim_mm_h, 'Color',colors(2,:), 'MarkerSize', 4, 'DisplayName', 'Simulated Q');
ylabel('Discharge [mm/h]','Interpreter','latex');
yyaxis right
area(time, rain_mm_h, 'FaceColor', colors(1,:), 'FaceAlpha', 0.5, 'DisplayName', 'Rainfall');
set(gca,'Ycolor','black')
ylim([0, max(rain_mm_h)*1.2]);
set(gca, 'YDir', 'reverse');
xlabel('Time','Interpreter','latex');
title('Observed vs Simulated Discharge with Rainfall');
legend('Location', 'northwest');
grid on;
set(gca,'TickDir', 'out', 'Box', 'on', 'FontName', 'Montserrat', 'FontSize', 12);

% Metrics box
annotation('textbox', [0.72 0.45 0.2 0.4], ...
    'String', metrics_text, 'FitBoxToText','on', ...
    'Interpreter','latex', 'FontSize', 11, ...
    'BackgroundColor','white', 'EdgeColor',[0.4 0.4 0.4]);

%% === Plot 2: Cumulative Volumes ===
cum_rain = cumsum(rain_mm_h) * dt_hr;
cum_etp  = cumsum(etp_mm_h)  * dt_hr;
cum_Qobs = cumsum(Q_obs_mm_h) * dt_hr;
cum_Qsim = cumsum(Q_sim_mm_h) * dt_hr;

figure('Color','w','Position',[100 600 1000 400]); hold on;
plot(time, cum_rain, '-', 'Color', pallete.blue_colors(1,:), 'LineWidth', 1.6, 'DisplayName', 'Rainfall');
plot(time, cum_etp,  '-', 'Color', pallete.red_colors(1,:), 'LineWidth', 1.6, 'DisplayName', 'ETP');
plot(time, cum_Qobs, '-','Color', pallete.green_colors(1,:), 'DisplayName', 'Observed Q');
plot(time, cum_Qsim, '--','Color', pallete.red_colors(2,:), 'DisplayName', 'Simulated Q');

xlabel('Time');
ylabel('Cumulative Depth [mm]');
title('Cumulative Rainfall, ETP, and Discharge');
legend('Location', 'northwest');
grid on;
set(gca, 'TickDir', 'out', 'Box', 'on', 'FontName', 'Montserrat', 'FontSize', 12);

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
    area(t_ev, r_ev, 'FaceColor', colors(1,:), 'FaceAlpha', 0.35, 'DisplayName', 'Rainfall');
    set(gca, 'YDir', 'reverse', 'YColor', 'k');
    ylim([0, 3 * r_max]);
    ylabel('Rainfall [mm/h]', 'FontName','Montserrat', 'Interpreter','latex');

    yyaxis left
    plot(t_ev, q_obs_ev, 'o', 'Color', colors(3,:), 'MarkerSize', 5, 'DisplayName', 'Obs'); hold on;
    plot(t_ev, q_sim_ev, 's', 'Color', colors(2,:), 'MarkerSize', 5, 'DisplayName', 'Sim');
    ylabel('Q [mm/h]', 'FontName','Montserrat', 'Interpreter','latex');

    title(sprintf('Top Event #%d (%.1f mm Rain)', p, event_volumes(sorted_idx(p))), ...
          'FontName', 'Montserrat', 'FontSize', 13, 'FontWeight', 'bold');
    
    grid on;
    legend('Location','northeast');
    set(gca, ...
        'TickDir', 'out', ...
        'Box', 'on', ...
        'FontName', 'Montserrat', ...
        'FontSize', 11, ...
        'XColor', 'k', ...
        'YColor', 'k');
end

sgtitle('Top Rainfall-Runoff Events (Hydrographs)', ...
        'FontWeight', 'bold', ...
        'FontSize', 16, ...
        'FontName', 'Montserrat');

end