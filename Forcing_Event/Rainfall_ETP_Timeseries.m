clear; clc;

%% === USER SETTINGS ===
dt_target_minutes = 5;  % 🔁 Change to resample to a different resolution

%% === DAILY ETP TIMESERIES ===
data_ETP = readtable('Catchment_Meteorological_Timeseries.xlsx');
date_ETP = datetime(table2array(data_ETP(:,2)));  % Ensure datetime
julian_day = day(date_ETP, 'dayofyear');
% Tavg = table2array(data_ETP(:,7));
Tmax = table2array(data_ETP(:,8));
Tmin = table2array(data_ETP(:,9));
Tavg = 1/2 * (Tmax + Tmin);
latitude = 29.5365;  % degrees
ET0_daily = max(hargreaves_ET0(Tavg, Tmax, Tmin, latitude, julian_day),0);  % [mm/day]

%% === RAINFALL TIMESERIES (with Gaps Filled) ===
data_rainfall = readtable('Observed_Rainfall.xlsx');
time_rain = datetime(data_rainfall{:,1});
rain_vals = data_rainfall{:,2};  % [mm over interval]

% 🧭 Define the original time step of the rainfall data (e.g., 5 minutes)
dt_input = minutes(5);  % ⏱ Adjust this to your actual native resolution

% ⏳ Create a continuous time vector from first to last timestamp at dt_input intervals
t_start = min(time_rain);
t_end   = max(time_rain);
time_filled = (t_start:dt_input:t_end)';

% 🔎 Fill missing timestamps with zero rainfall
% Create a timetable for known data
rain_tt = timetable(time_rain, rain_vals);

% Re-time to regular interval and fill missing with zeros
rain_tt_filled = retime(rain_tt, time_filled, 'fillwithconstant', 'Constant', 0);

% Convert to regular variables
time_rain_filled = rain_tt_filled.time_rain;
rain_mm_per_hr   = rain_tt_filled.rain_vals / (minutes(dt_input)/60);  % Convert to mm/hr

% === ETP Interpolation ===
ET0_hourly_interp = interp1(date_ETP, ET0_daily / 24, time_rain_filled, 'linear', 'extrap');

%% === RESAMPLING TO TARGET TIMESTEP ===
% 🕓 Define target timestep (e.g., 1440 = daily, 60 = hourly)
dt_target = minutes(dt_target_minutes);
time_uniform = (t_start:dt_target:t_end)';

% 🌧 Resample rainfall to uniform time grid (integrating over interval)
rain_tt_uniform = retime(timetable(time_rain_filled, rain_mm_per_hr), ...
                         time_uniform, 'mean');
rain_interp_mm_per_hr = rain_tt_uniform.rain_mm_per_hr;

% ☀️ Resample ET0 to average per interval
etp_tt_uniform = retime(timetable(time_rain_filled, ET0_hourly_interp), ...
                        time_uniform, 'mean');
ET0_interp_mm_per_hr = etp_tt_uniform.ET0_hourly_interp;

% === Convert to mm per timestep
dt_hr = dt_target_minutes / 60;
rain = rain_interp_mm_per_hr * dt_hr;
etp  = ET0_interp_mm_per_hr  * dt_hr;


%% === STATISTICS ===
total_rain = sum(rain, 'omitnan');
total_etp = sum(etp, 'omitnan');
mean_rain = mean(rain, 'omitnan');
mean_etp = mean(etp, 'omitnan');
max_rain = max(rain);
min_rain = min(rain);
max_etp = max(etp);
std_rain = std(rain, 'omitnan');
std_etp = std(etp, 'omitnan');

years = year(time_uniform);
unique_years = unique(years);
n_years = numel(unique_years);
annual_rain = zeros(n_years,1);
annual_etp = zeros(n_years,1);

for i = 1:n_years
    idx = years == unique_years(i);
    if idx == 5
        ttt = 1;
    end
    annual_rain(i) = sum(rain(idx), 'omitnan');
    if annual_rain(i) == 0
        ttt = 1;
    end
    annual_etp(i) = sum(etp(idx), 'omitnan');
end

mean_annual_rain = mean(annual_rain, 'omitnan');
mean_annual_etp = mean(annual_etp, 'omitnan');
aridity_index = mean_annual_etp / mean_annual_rain;

%% === DISPLAY RESULTS ===
fprintf('\n--- Summary Statistics ---\n');
fprintf('Total Rainfall:              %.2f mm\n', total_rain);
fprintf('Total ETP:                   %.2f mm\n', total_etp);
fprintf('Mean Rainfall per step:      %.4f mm\n', mean_rain);
fprintf('Mean ETP per step:           %.4f mm\n', mean_etp);
fprintf('Max Rainfall per step:       %.2f mm\n', max_rain);
fprintf('Min Rainfall per step:       %.2f mm\n', min_rain);
fprintf('Max ETP per step:            %.2f mm\n', max_etp);
fprintf('Std Rainfall per step:       %.4f mm\n', std_rain);
fprintf('Std ETP per step:            %.4f mm\n', std_etp);
fprintf('Mean Annual Rainfall:        %.2f mm\n', mean_annual_rain);
fprintf('Mean Annual ETP:             %.2f mm\n', mean_annual_etp);
fprintf('Aridity Index (ETP/Rain):    %.2f\n', aridity_index);

%% === PLOT ===
figure;

% Global figure settings
set(groot, 'defaultAxesFontName', 'Montserrat');
set(groot, 'defaultTextFontName', 'Montserrat');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% Plot Rainfall and ET0
yyaxis left
plot(time_uniform, rain_interp_mm_per_hr, 'b-', 'LineWidth', 1.5);
ylabel('Rainfall [mm/h]', 'Interpreter', 'latex');

yyaxis right
plot(time_uniform, ET0_interp_mm_per_hr, 'r-', 'LineWidth', 1.5);
ylabel('Interpolated $ET_0$ [mm/h]', 'Interpreter', 'latex');

xlabel('Time', 'Interpreter', 'latex');
title('Rainfall and Downscaled $ET_0$', 'Interpreter', 'latex');
legend({'Rainfall', '$ET_0$'}, 'Location', 'best');

% Axes styling
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Montserrat';
ax.Box = 'on';
ax.LineWidth = 2;
ax.TickDir = 'out';
grid on;

% Optional export
% print(gcf, 'Rainfall_ET0_Plot', '-dpdf', '-r300');

%% === PLOT: Annual Rainfall and ETP with Std Dev ===

figure('Units', 'inches', 'Position', [1, 1, 6.7, 5.8]);  % [left bottom width height]

% === Set paper properties for perfect export ===
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.7 5.8]);
set(gcf, 'PaperSize', [6.7 5.8]);

% Global figure settings
set(groot, 'defaultAxesFontName', 'Montserrat');
set(groot, 'defaultTextFontName', 'Montserrat');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% --- Calculate Means and Standard Deviations ---
mean_annual_rain = mean(annual_rain, 'omitnan');
std_annual_rain  = std(annual_rain, 'omitnan');

mean_annual_etp  = mean(annual_etp, 'omitnan');
std_annual_etp   = std(annual_etp, 'omitnan');

% === PLOT: Annual Rainfall and ETP without Error Bars ===
subplot(2,1,1)
hold on;

% Rainfall
plot(unique_years, annual_rain, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);

% ETP
plot(unique_years, annual_etp, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);

% Labels and title
xlabel('Year', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Annual Volume [mm]', 'Interpreter', 'latex', 'FontSize', 14);
title('Annual Rainfall and Potential ET', 'Interpreter', 'latex', 'FontSize', 16);

% Create annotation text
ai_text = sprintf('Aridity Index: %.2f', aridity_index);

% Add annotation box
annotation('textbox', [0.15, 0.69, 0.2, 0.1], ...
           'String', ai_text, ...
           'Interpreter', 'latex', ...
           'FontSize', 12, ...
           'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', ...
           'EdgeColor', [0 0 0]);

% Add a legend
legend({'Annual Rainfall', 'Annual ETP'}, 'Location', 'northwest', 'FontSize', 12);

% Formatting the axes
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Montserrat';
ax.Box = 'on';
ax.LineWidth = 2;
ax.TickDir = 'out';
grid on;


% === PLOT: Monthly Cumulative Rainfall and ETP (Mean and Std) ===
subplot(2,1,2)
% Extract year and month
year_vec = year(time_uniform);
month_vec = month(time_uniform);

% Identify unique years
unique_years = unique(year_vec);
n_years = numel(unique_years);

% Initialize matrices [year x month]
monthly_rain_totals = NaN(n_years,12);
monthly_etp_totals = NaN(n_years,12);

for i = 1:n_years
    for m = 1:12
        idx = (year_vec == unique_years(i)) & (month_vec == m);
        monthly_rain_totals(i,m) = sum(rain(idx), 'omitnan');
        monthly_etp_totals(i,m)  = sum(etp(idx), 'omitnan');
    end
end

% Now calculate mean and std across years
monthly_rain_mean = mean(monthly_rain_totals, 1, 'omitnan');
monthly_rain_std  = std(monthly_rain_totals, 0, 1, 'omitnan');
monthly_etp_mean  = mean(monthly_etp_totals, 1, 'omitnan');
monthly_etp_std   = std(monthly_etp_totals, 0, 1, 'omitnan');

% Global figure settings
set(groot, 'defaultAxesFontName', 'Montserrat');
set(groot, 'defaultTextFontName', 'Montserrat');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

hold on;

% Plot Rainfall
errorbar(1:12, monthly_rain_mean, monthly_rain_std, 'bo-', 'LineWidth', 1.8, 'MarkerSize', 7, 'CapSize', 10);
ylabel('Monthly Volume [mm]', 'Interpreter', 'latex');

% Plot ETP
errorbar(1:12, monthly_etp_mean, monthly_etp_std, 'rs-', 'LineWidth', 1.8, 'MarkerSize', 7, 'CapSize', 10);

% Common settings
xlabel('Month', 'Interpreter', 'latex');
xlim([0.5 12.5]);
xticks(1:12);
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'});
title('Monthly Cumulative Rainfall and Potential ET', 'Interpreter', 'latex');

legend({'Rainfall', 'Potential ET'}, 'Location', 'northwest');

% Format axes
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.TickDir = 'out';
grid on;

% Optional: Save the plot
% print(gcf, 'Monthly_Cumulative_Rainfall_ETP_StdDev', '-dpdf', '-r300');

%% === EXPORT TIMESERIES TO EXCEL ===

% Create table
timeseries_table = table(time_uniform, rain_interp_mm_per_hr, ET0_interp_mm_per_hr, ...
    'VariableNames', {'Time', 'Rain_mm_per_h', 'ETP_mm_per_h'});

% Ensure the Modeling_Results folder exists
output_folder = 'Modeling_Results';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Define output filename
output_filename = fullfile(output_folder, 'Rainfall_ETP_Timeseries.xlsx');

% Write table to Excel file
writetable(timeseries_table, output_filename);

fprintf('Timeseries successfully exported to %s\n', output_filename);

%% === DETECT RAINFALL EVENTS AND EXPORT TO EXCEL WITH INTENSITY ===

% === USER SETTINGS ===
min_event_volume = 0.5;           % Minimum event volume [mm]
inter_event_hrs = 6;              % Minimum inter-event dry period [hours]
lookback_days = 5;                % Antecedent period for rainfall summary [days]
dt_hr = dt_target_minutes / 60;   % Target timestep in hours
inter_event_steps = inter_event_hrs * 60 / dt_target_minutes;
margin = inter_event_steps;       % Margin before and after event for plotting/context

% === STEP 1: IDENTIFY POTENTIAL RAINFALL EVENTS ===
is_rain = rain > 0;
rain_event_idx = bwlabel(is_rain);  % Label continuous rain periods
n_events = max(rain_event_idx);     % Number of detected rain clusters

% === Prepare Excel Output ===
excel_filename = fullfile('Modeling_Results', 'Rainfall_ETP_Events.xlsx');
if exist(excel_filename, 'file'), delete(excel_filename); end  % Overwrite old file

% === Initialize Storage for Summary and Time Series ===
summary = [];                    % Summary statistics for each event
events_data = struct();          % Struct to store individual event time series
event_counter = 0;               % Track number of valid events exported

% === STEP 2: LOOP THROUGH EVENTS AND FILTER VALID ONES ===
for i = 1:n_events
    idx = find(rain_event_idx == i);  % Index of current rain cluster
    if isempty(idx), continue; end

    % Check if there's enough dry time before this event
    if i > 1
        prev_idx = find(rain_event_idx == i-1, 1, 'last');
        gap = idx(1) - prev_idx;
        if gap < inter_event_steps, continue; end
    end

    % Check if total rainfall volume meets threshold
    event_rain_volume = sum(rain(idx));
    if event_rain_volume < min_event_volume, continue; end

    % === Event is valid ===
    event_counter = event_counter + 1;

    % Define event duration and context window
    start_idx = idx(1);
    end_idx = idx(end);
    start_time = time_uniform(start_idx);
    end_time = time_uniform(end_idx);
    duration_hrs = hours(end_time - start_time);

    % Define extended window (margin before and after)
    i1 = max(1, start_idx - margin);
    i2 = min(length(rain), end_idx + margin);
    t_event = time_uniform(i1:i2);
    r_event = rain(i1:i2);
    e_event = etp(i1:i2);

    % === Compute Rainfall and ETP Stats ===
    R_total = sum(r_event);
    R_max   = max(r_event);
    R_mean  = mean(r_event);
    R_min   = min(r_event);
    R_std   = std(r_event);

    E_total = sum(e_event);
    E_max   = max(e_event);
    E_mean  = mean(e_event);
    E_min   = min(e_event);
    E_std   = std(e_event);

    rain_intensity_mean = R_total / duration_hrs;
    rain_etp_ratio = R_total / max(E_total, eps);

    % === Antecedent Dry Days Calculation ===
    dry_days = 0;
    check_idx = start_idx - 1;
    while check_idx > 0 && rain(check_idx) == 0
        dry_days = dry_days + dt_target_minutes / (60 * 24);  % Increment in days
        check_idx = check_idx - 1;
    end
    dry_days = round(dry_days, 2);  % Round for display

    % === 5-Day Antecedent Rainfall Volume ===
    lookback_steps = lookback_days * 24 * 60 / dt_target_minutes;
    start_5day_idx = max(1, start_idx - lookback_steps);
    antecedent_rain_5d = sum(rain(start_5day_idx:start_idx - 1), 'omitnan');

    % === Store Summary Row ===
    summary = [summary; table(event_counter, start_time, end_time, duration_hrs, ...
        R_total, R_max, R_mean, R_min, R_std, ...
        E_total, E_max, E_mean, E_min, E_std, ...
        rain_intensity_mean, rain_etp_ratio, dry_days, antecedent_rain_5d)];

    % === Prepare Time Series Table for Export ===
    rain_intensity = r_event / dt_hr;  % [mm/h]
    etp_intensity  = e_event / dt_hr;  % [mm/h]
    elapsed_min = minutes(t_event - t_event(1));  % Time since window start

    T_event = table(t_event, elapsed_min, rain_intensity, etp_intensity, ...
        'VariableNames', {'Time', 'Elapsed_Time_min', ...
                          'Rain_Intensity_mm_per_h', 'ETP_Intensity_mm_per_h'});

    % Create safe sheet name for Excel
    sheet_name = sprintf('Ev%d_%s', event_counter, datestr(start_time, 'mmmdd_HHMM'));
    sheet_name = regexprep(sheet_name, '[^a-zA-Z0-9]', '_');

    % Store in struct for delayed writing
    events_data(event_counter).sheet_name = sheet_name;
    events_data(event_counter).table = T_event;
end

% === STEP 3: WRITE SUMMARY TO FIRST SHEET ===
summary.Properties.VariableNames = {'EventID', 'StartTime', 'EndTime', 'Duration_hr', ...
    'Rain_Total_mm', 'Rain_Max_mm', 'Rain_Mean_mm', 'Rain_Min_mm', 'Rain_Std_mm', ...
    'ETP_Total_mm', 'ETP_Max_mm', 'ETP_Mean_mm', 'ETP_Min_mm', 'ETP_Std_mm', ...
    'Rain_Intensity_mm_per_hr', 'Rain_ETP_Ratio', ...
    'Antecedent_Dry_Days', 'Antecedent_Rain_5d_mm'};

writetable(summary, excel_filename, 'Sheet', 1, 'WriteRowNames', false);

% === STEP 4: WRITE EACH EVENT TO ITS OWN SHEET ===
for k = 1:numel(events_data)
    writetable(events_data(k).table, excel_filename, ...
        'Sheet', events_data(k).sheet_name);
end

fprintf('✅ %d rainfall events exported with summary to "%s"\n', event_counter, excel_filename);

