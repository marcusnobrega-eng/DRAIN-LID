clear; clc;

%% === USER SETTINGS ===
dt_target_minutes = 15;  % 🔁 Change to resample to a different resolution

output_folder = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results';

%% === DAILY ETP TIMESERIES ===
data_ETP = readtable('Catchment_Meteorological_Timeseries_Monitored_PP.xlsx');
date_ETP = datetime(table2array(data_ETP(:,2)));  % Ensure datetime
julian_day = day(date_ETP, 'dayofyear');
% Tavg = table2array(data_ETP(:,7));
Tmax = table2array(data_ETP(:,8));
Tmin = table2array(data_ETP(:,9));
Tavg = 1/2 * (Tmax + Tmin);
latitude = 29.5365;  % degrees
ET0_daily = max(hargreaves_ET0(Tavg, Tmax, Tmin, latitude, julian_day),0);  % [mm/day]

%% === RAINFALL TIMESERIES (with Gaps Filled) ===
% Load Rainfall Data
data_rainfall = readtable('Observed_Rainfall_PP.xlsx');
time_rain = datetime(data_rainfall{:,1});
rain_vals = data_rainfall{:,2};  % [mm over interval]

% Define Input Time Step
dt_input = minutes(15);  % 15-minute resolution

% Create Regular Time Vector
t_start = min(time_rain);
t_end = max(time_rain);
time_filled = (t_start:dt_input:t_end)';

% Build Timetable
rain_tt = timetable(time_rain, rain_vals);

% Re-time using 'nearest' to avoid losing data
rain_tt_filled = retime(rain_tt, time_filled, 'nearest');

% Optional: If you prefer to keep gaps and fill missing manually:
% rain_tt_filled = synchronize(rain_tt, timetable(time_filled), 'first', 'fillwithconstant', 'Constant', 0);

% Extract Variables
time_rain_filled = rain_tt_filled.time_rain;
rain_mm_per_hr = rain_tt_filled.rain_vals / (minutes(dt_input)/60);  % Convert mm/interval to mm/hr

% Load ET0 Data
% Assuming date_ETP is datetime and ET0_daily is in mm/day
% date_ETP = datetime(...);
% ET0_daily = ...;

ET0_hourly = ET0_daily / 24;  % mm/hr assumption
ET0_hourly_interp = interp1(date_ETP, ET0_hourly, time_rain_filled, 'linear', 'extrap');


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

%% === OBSERVED DISCHARGE (cfs ➜ m³/s, irregular ➜ regular) ===

% ⚠️ Make sure your discharge file is formatted as:
%  • Column 1: datetime (Excel or datetime format)
%  • Column 2: discharge in cfs

% === Load and parse observed discharge ===
data_discharge = readtable('Observed_Discharge_PP.xlsx');
time_discharge = datetime(data_discharge{:,1});        % time
Q_obs_cfs = data_discharge{:,2};                       % discharge in cfs

% === Convert to SI units (1 cfs = 0.0283168 m³/s)
Q_obs_m3s = Q_obs_cfs * 0.0283168;

% === Interpolate to the same uniform time vector as rainfall/ETP
Q_obs_interp = interp1(time_discharge, Q_obs_m3s, time_uniform, 'linear', 'extrap');

% === Optional NaN handling for extrapolation beyond actual data range
Q_obs_interp(isnan(Q_obs_interp)) = 0;  % or use fillmissing(...)

% === PLOT: Original vs Interpolated Discharge ===

figure('Color','w', 'Position', [300, 300, 800, 400]);

% Plot original discharge data
plot(time_discharge, Q_obs_m3s, 'ko', 'MarkerSize', 4, 'DisplayName', 'Original Data'); hold on;

% Plot interpolated discharge
plot(time_uniform, Q_obs_interp, 'r-', 'LineWidth', 1.4, 'DisplayName', 'Interpolated Data');

xlabel('Time');
ylabel('Discharge [m³/s]');
title('Observed Discharge: Original vs Interpolated');
legend('Location', 'best');
grid on;



% === Append to existing exported table
timeseries_table = table(time_uniform, ...
                         rain_interp_mm_per_hr, ...
                         ET0_interp_mm_per_hr, ...
                         Q_obs_interp, ...
    'VariableNames', {'Time', 'Rain_mm_per_h', 'ETP_mm_per_h', 'Discharge_m3_per_s'});

% === Define output path
output_filename = fullfile(output_folder, 'Rainfall_ETP_Discharge_Timeseries.xlsx');

% === Export unified timeseries
writetable(timeseries_table, output_filename);

fprintf('📤 Timeseries with discharge exported to: %s\n', output_filename);

%% === RUNOFF COEFFICIENT AND TOTAL VOLUMES ===
catchment_area_m2 = 195.1;  % 📐 Surface area of the watershed in m²

% Total rainfall volume [m³] = sum(mm) * area [m²] * 1e-3
total_rainfall_volume_m3 = sum(rain, 'omitnan') * catchment_area_m2 * 1e-3;

% Total discharge volume [m³] = sum(flow [m³/s] * Δt [s])
dt_sec = dt_target_minutes * 60;
total_discharge_volume_m3 = sum(Q_obs_interp, 'omitnan') * dt_sec;

% Runoff coefficient = discharge volume / rainfall volume
runoff_coefficient = total_discharge_volume_m3 / total_rainfall_volume_m3;

% Display new metrics
fprintf('--- Hydrologic Volume Stats ---\n');
fprintf('Total Rainfall Volume     : %.2f m³\n', total_rainfall_volume_m3);
fprintf('Total Discharge Volume    : %.2f m³\n', total_discharge_volume_m3);
fprintf('Runoff Coefficient (Q/P)  : %.4f\n', runoff_coefficient);


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

% Create combined time series table with discharge included
timeseries_table = table(time_uniform, ...
                         rain_interp_mm_per_hr, ...
                         ET0_interp_mm_per_hr, ...
                         Q_obs_interp, ...
    'VariableNames', {'Time', 'Rain_mm_per_h', 'ETP_mm_per_h', 'Discharge_m3_per_s'});


% Ensure the Modeling_Results folder exists
output_folder = 'Modeling_Results';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Define output filename
output_filename = fullfile(output_folder, 'Rainfall_ETP_Discharge_Timeseries.xlsx');

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

% === STEP 1: IDENTIFY RAINFALL EVENTS WITH PROPER INTER-EVENT DRY PERIOD ===
event_id = zeros(size(rain));  % Vector to store event IDs
event_counter = 0;
i = 1;
while i <= length(rain)
    if rain(i) > 0
        % Start new event
        event_counter = event_counter + 1;
        event_start = i;
        i = i + 1;
        while i <= length(rain)
            if rain(i) == 0
                remaining_steps = length(rain) - i;
                if (i + inter_event_steps) <= length(rain)
                    if all(rain(i+1:i+inter_event_steps) == 0)
                        break;
                    end
                else
                    break;
                end
            end
            i = i + 1;
        end
    event_end = min(i, length(rain));  % Ensure we stay within bounds
    event_id(event_start:event_end) = event_counter;

    end
    i = i + 1;
end

% Get list of valid event IDs
unique_events = unique(event_id(event_id > 0));


% === Prepare Excel Output ===
excel_filename = fullfile('Modeling_Results', 'Rainfall_ETP_Discharge_Events.xlsx');
if exist(excel_filename, 'file'), delete(excel_filename); end  % Overwrite old file

% === Initialize Storage for Summary and Time Series ===
summary = [];                    % Summary statistics for each event
events_data = struct();          % Struct to store individual event time series
event_counter = 0;               % Track number of valid events exported

% === STEP 2: LOOP THROUGH EVENTS AND FILTER VALID ONES ===
for i = 1:numel(unique_events)
    idx = find(event_id == unique_events(i));
    if isempty(idx), continue; end

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
    duration_hrs = (end_idx - start_idx + 1) * dt_hr;


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

    % === Get observed discharge in m³/s for the event window ===
    Q_event_m3s = Q_obs_interp(i1:i2);
    
    % === Convert to discharge per unit area in m/s ===
    Q_event_ms = Q_event_m3s / catchment_area_m2;

    % === Convert discharge from m³/s to mm/h over the watershed
    Q_event_mm_per_h = Q_event_m3s / catchment_area_m2 * 3600 * 1000;  % [mm/h]
    
    % === Compute rainfall and discharge volume over the event
    rain_volume_m3 = R_total * catchment_area_m2 * 1e-3;               % [m³]
    Q_volume_m3 = sum(Q_event_m3s) * dt_hr * 3600;                     % [m³]
    
    % === Runoff Coef
    runoff_coef = Q_volume_m3 / rain_volume_m3;
    
    % === Store Summary Row with added Q stats and rain/Q volume ratio
    summary = [summary; table(event_counter, start_time, end_time, duration_hrs, ...
        R_total, R_max, R_mean, R_min, R_std, ...
        E_total, E_max, E_mean, E_min, E_std, ...
        rain_intensity_mean, rain_etp_ratio, ...
        dry_days, antecedent_rain_5d, ...
        max(Q_event_mm_per_h), mean(Q_event_mm_per_h), ...
        min(Q_event_mm_per_h), std(Q_event_mm_per_h), ...
        runoff_coef)];


    % === Prepare Time Series Table for Export ===
    rain_intensity = r_event / dt_hr;  % [mm/h]
    etp_intensity  = e_event / dt_hr;  % [mm/h]
    elapsed_min = minutes(t_event - t_event(1));  % Time since window start

    
    % === Build event export table ===
    T_event = table(t_event, elapsed_min, ...
        rain_intensity, etp_intensity, ...
        Q_event_ms * 1000 * 3600, Q_event_m3s, ...
        'VariableNames', {'Time', 'Elapsed_Time_min', ...
                          'Rain_Intensity_mm_per_h', ...
                          'ETP_Intensity_mm_per_h', ...
                          'Discharge_Observed_mm_per_h', ...
                          'Discharge_Observed_m3_per_s'});



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
    'Antecedent_Dry_Days', 'Antecedent_Rain_5d_mm', ...
    'Q_Max_mm_per_h', 'Q_Mean_mm_per_h', 'Q_Min_mm_per_h', 'Q_Std_mm_per_h', ...
    'runoff_coef'};


writetable(summary, excel_filename, 'Sheet', 1, 'WriteRowNames', false);

% === STEP 4: WRITE EACH EVENT TO ITS OWN SHEET ===
for k = 1:numel(events_data)
    writetable(events_data(k).table, excel_filename, ...
        'Sheet', events_data(k).sheet_name);
end

fprintf('✅ %d rainfall events exported with summary to "%s"\n', event_counter, excel_filename);

