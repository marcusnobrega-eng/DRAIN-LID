% Example parameters:
a = 0.006;       % Coefficient for rating curve
b = 2.5532;      % Exponent
evap = 0;        % Evaporation losses in m³
storage = 0;     % Remaining storage in m³
catchment_area_m2 = 195.1;

optimize_flow_bias('pressure_and_rainfall_data.xlsx', a, b, catchment_area_m2, evap, storage)

function [bias_all, fval] = optimize_flow_bias(filename, a, b, catchment_area_m2, evaporation_loss, remaining_storage);
    %% === READ INPUT DATA ===
    raw = readtable(filename);
    t_level = datetime(raw{:,1});
    h_raw = raw{:,2}; % bubbler level in inches
        
    t_rain = datetime(raw{:,3});
    r_raw = raw{:,4}; % rainfall in mm/h
    
    % Limiting maximum level to 13 inches
    h_raw(h_raw > 13) = 12;

    rain_5min = r_raw * (5/60);  % Convert to mm per 5 min

    %% === CLEAN & RESAMPLE ===
    % Clean non-finite
    valid_h = isfinite(h_raw) & isfinite(datenum(t_level));
    valid_r = isfinite(rain_5min) & isfinite(datenum(t_rain));
    t_level_clean = t_level(valid_h);
    h_clean = h_raw(valid_h);
    t_rain_clean = t_rain(valid_r);
    r_clean = rain_5min(valid_r);

    % Create uniform time grid
    t_start = min([min(t_level_clean), min(t_rain_clean)]);
    t_end   = max([max(t_level_clean), max(t_rain_clean)]);
    t_uniform = (t_start:minutes(5):t_end)';

    % Interpolate to 5-min intervals, fill missing with 0
    tt_h = timetable(t_level_clean, h_clean);
    tt_r = timetable(t_rain_clean, r_clean);
    tt_h_uniform = retime(tt_h, t_uniform, 'fillwithconstant', 'Constant', 0);
    tt_r_uniform = retime(tt_r, t_uniform, 'fillwithconstant', 'Constant', 0);
    h_interp = tt_h_uniform.h_clean;
    r_interp = tt_r_uniform.r_clean;

    % Build unified timetable
    ts = timetable(t_uniform, h_interp, r_interp);
    ts.Properties.VariableNames = {'Level_in', 'Rain_mm'};

    %% === EVENT DETECTION ===
    inter_event_hours = 6;
    min_event_depth_mm = 1;
    dt_minutes = minutes(median(diff(t_uniform)));
    inter_event_steps = round(inter_event_hours * 60 / dt_minutes);
    drainage_extension_steps = round(6 * 60 / dt_minutes);  % ±3h extension

    event_id = zeros(height(ts), 1);
    ev_count = 0;
    i = 1;

    while i <= height(ts)
        if ts.Rain_mm(i) > 0
            ev_count = ev_count + 1;
            ev_start = i;
            i = i + 1;
            while i <= height(ts)
                if ts.Rain_mm(i) == 0
                    if (i + inter_event_steps <= height(ts)) && all(ts.Rain_mm(i+1:i+inter_event_steps) == 0)
                        break;
                    end
                end
                i = i + 1;
            end
            ev_end = min(i, height(ts));
            event_id(ev_start:ev_end) = ev_count;
        end
        i = i + 1;
    end

    % === Filter events by min rainfall ===
    unique_ids = unique(event_id);
    unique_ids(unique_ids == 0) = [];
    valid_event_id = zeros(size(event_id));
    new_ev_count = 0;

    for k = 1:length(unique_ids)
        idx = find(event_id == unique_ids(k));
        total_rain = sum(ts.Rain_mm(idx), 'omitnan');
        if total_rain >= min_event_depth_mm
            new_ev_count = new_ev_count + 1;
            valid_event_id(idx) = new_ev_count;
        end
    end

    event_id = valid_event_id;
    ev_count = new_ev_count;

    %% === OPTIMIZATION PER EVENT ===
    bias_all = zeros(ev_count,1);
    Q_mm_per_hr = zeros(height(ts),1);
    rain_vol = zeros(ev_count,1);

    dt_sec = seconds(median(diff(t_uniform)));

    for k = 1:ev_count
        idx_core = find(event_id == k);
        if isempty(idx_core), continue; end

        % Extended window
        idx_start = max(1, min(idx_core) - drainage_extension_steps);
        idx_end   = min(height(ts), max(idx_core) + drainage_extension_steps);
        idx_full  = idx_start:idx_end;

        subTS = ts(idx_full, :);

        % Rain volume (only over core event)
        rain_vol(k) = sum(ts.Rain_mm(idx_core), 'omitnan') / 1000 * catchment_area_m2;

        if k == 10
            ttt = 1;
        end

        % Objective function for bias (mass balance)
        obj_fun = @(bias) volume_difference(bias, subTS.Level_in, a, b, catchment_area_m2, rain_vol(k), evaporation_loss, remaining_storage, dt_sec);
        [bias_all(k), fval(k)]= fminbnd(obj_fun, -30, 5);  % Bias in inches

        % Apply bias
        h_m = (subTS.Level_in + bias_all(k));  % still in inches
        h_m(h_m < 0) = 0;
        Q = (0.3048^3) * (a * h_m .^ b);  % cfs to m³/s

        Q_event_mm_per_hr = Q / catchment_area_m2 * 3600 * 1000;  % mm/h
        Q_mm_per_hr(idx_full) = Q_event_mm_per_hr;
    end

    %% === OUTPUT ===
    result_ts = timetable(t_uniform, Q_mm_per_hr);
    result_ts.Properties.VariableNames = {'Discharge_Estimated_mm_per_hr'};
    writetable(timetable2table(result_ts), 'Corrected_Discharge_mm_per_hr.xlsx');
    fprintf('✅ Exported corrected discharge timeseries to Corrected_Discharge_mm_per_hr.xlsx\n');

    %% === RUNOFF RATIO AND PLOTS ===
    Q_m3s = Q_mm_per_hr / 1000 / 3600 * catchment_area_m2;
    total_discharge_m3 = sum(Q_m3s) * dt_sec;
    total_rainfall_m3 = sum(ts.Rain_mm, 'omitnan') / 1000 * catchment_area_m2;
    runoff_ratio = total_discharge_m3 / total_rainfall_m3;

    fprintf('--- Runoff Ratio Summary ---\n');
    fprintf('Total Rainfall Volume  : %.2f m³\n', total_rainfall_m3);
    fprintf('Total Discharge Volume : %.2f m³\n', total_discharge_m3);
    fprintf('Runoff Ratio (Q / P)   : %.4f\n', runoff_ratio);

    % Plot cumulative rainfall and runoff
    cumulative_rain_mm = cumsum(ts.Rain_mm, 'omitnan');
    cumulative_Q_mm = cumsum(Q_mm_per_hr, 'omitnan') * dt_minutes / 60;  % mm/h * hours

    figure('Color','w','Position',[200 300 800 400]);
    plot(t_uniform, cumulative_rain_mm, 'b-', 'LineWidth', 2); hold on;
    plot(t_uniform, cumulative_Q_mm, 'r-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Cumulative Depth [mm]');
    title('Cumulative Rainfall and Estimated Discharge');
    legend('Cumulative Rainfall', 'Cumulative Discharge', 'Location', 'best');
    grid on;
end

function diff = volume_difference(bias, level_in, a, b, area, rain_vol, evap, storage, dt_sec)
    h = (level_in + bias);  % inches
    h(h < 0) = 0;
    Q = a * h .^ b;         % cfs
    Q = (0.3048^3) * Q;     % m³/s
    outflow_vol = sum(Q) * dt_sec;
    diff = abs((0.90 * rain_vol - evap - storage) - outflow_vol);
end
