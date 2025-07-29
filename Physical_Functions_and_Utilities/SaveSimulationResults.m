%% =========================================================================
% 📦 EXPORT SCRIPT: SaveSimulationResults.m
%
% Purpose  : Saves simulation outputs (pressure head, moisture, fluxes, etc.)
%            into a single Excel file with proper time axes and physical units.
%
% Location : Add to the end of your Main Solver
%
% Output   : Results_<sim_name>/Data/SimulationResults.xlsx
%
% Includes :
%   - Head at all nodes        [m]
%   - Moisture θ at all nodes  [m³/m³]
%   - Darcy fluxes at all interfaces [m/s]
%   - Ponding depth            [m]
%   - Outlet flux              [m/s]
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
%% =========================================================================

% === Prepare Time Columns ===============================================
time_seconds = time_series';         % [s]
time_minutes = time_seconds / 60;    % [min]
time_days    = time_seconds / 86400; % [days]
data_length  = length(time_seconds); % [length]

% === Create Output Excel File Path ======================================
filename = fullfile(data_dir, 'SimulationResults.xlsx');

% === Delete existing file if it exists ==================================
if exist(filename, 'file')
    delete(filename);
end

%% === 1. Pressure Head h(z,t) ============================================
T_h = array2table([time_seconds, time_minutes, time_days, head_out(:,1:data_length)'], ...
    'VariableNames', ...
    ['Time_s', 'Time_min', 'Time_day', ...
     strcat("h_Node_m_", string(1:params.Nz))]);

% Optional metadata (not displayed, but readable with readtable)
T_h.Properties.VariableUnits = [ ...
    {'s', 'min', 'day'}, repmat({'m'}, 1, size(head_out,1))];

% Write to Excel
writetable(T_h, filename, 'Sheet', 'Head', 'WriteVariableNames', true);

%% === 2. Volumetric Moisture θ(z,t) ======================================
T_theta = array2table([time_seconds, time_minutes, time_days, theta_out(:,1:data_length)'], ...
    'VariableNames', ...
    ['Time_s', 'Time_min', 'Time_day', ...
     strcat("theta_Node_", string(1:params.Nz))]);

% Optional metadata (not displayed, but readable with readtable)
T_theta.Properties.VariableUnits = [ ...
    {'s', 'min', 'day'}, repmat({'m3/m3'}, 1, size(theta_out,1))];


writetable(T_theta, filename, 'Sheet', 'Moisture', 'WriteVariableNames', true);

%% === 3. Vertical Flux q(z,t) ============================================
T_q = array2table([time_seconds, time_minutes, time_days, flux_out(:,1:data_length)'], ...
    'VariableNames', ...
    ['Time_s', 'Time_min', 'Time_day', ...
     strcat("q_Interface_m_per_sec_", string(1:params.Nz+1))]);

% Optional metadata (not displayed, but readable with readtable)
T_q.Properties.VariableUnits = [ ...
    {'s', 'min', 'day'}, repmat({'m/s'}, 1, params.Nz + 1)];

writetable(T_q, filename, 'Sheet', 'Flux', 'WriteVariableNames', true);

%% === 4. Surface and Outlet Variables ====================================
T_extra = table(time_seconds, time_minutes, time_days, ...
                ponding_series(1:data_length)', seepage_flux(1:data_length)', top_flux(1:data_length)', ...
                Q_orifice_total(1:data_length)', Q_spillway_total(1:data_length)', ...
                'VariableNames', {'Time_s', 'Time_min', 'Time_day', ...
                                  'PondingDepth_m', 'Seepage_Flux_m_per_s', 'Top_Flux_m_per_s','Q_orifice_m_per_s','Q_spillway_m_per_s'});

T_extra.Properties.VariableUnits = ...
    {'s', 'min', 'day', 'm', 'm/s', 'm/s', 'm/s', 'm/s'};

writetable(T_extra, filename, 'Sheet', 'SurfaceOutlet', 'WriteVariableNames', true);

%% === 5. Volume Balance Summary (Scalar Values) ==========================
T_vol = table(inflow_vol, outflow_vol, seepage_vol, evaporation_vol, final_storage, ...
    'VariableNames', {'InflowVol_m', 'OutflowVol_m', 'SeepageVol_m', 'EvaporationVol_m', 'final_storage_m'});

T_vol.Properties.VariableUnits = {'m', 'm', 'm', 'm', 'm'};

% Write to Excel
writetable(T_vol, filename, 'Sheet', 'VolumeBalance', 'WriteVariableNames', true);

%% ✅ Done
fprintf('\n📁 Results saved to: %s\n', filename);

%% === 6. Save Simulation Info to Text File ================================

info_txt = fullfile(data_dir, 'SimulationInfo.txt');
fid = fopen(info_txt, 'w');

fprintf(fid, '============================================\n');
fprintf(fid, '   MIXED-FORM RICHARDS SIMULATION SUMMARY\n');
fprintf(fid, '============================================\n');
fprintf(fid, 'Simulation Name       : %s\n', sim_name);
fprintf(fid, 'Date                  : %s\n', datestr(now));
fprintf(fid, '\n--- Domain & Discretization ---\n');
fprintf(fid, 'Number of Nodes (Nz)  : %d\n', params.Nz);
fprintf(fid, 'Domain Depth (L) [m]  : %.3f\n', params.L);
fprintf(fid, 'Grid Refinement (n)   : %.2f\n', nonlin_factor);

fprintf(fid, '\n--- Time Settings ---\n');
fprintf(fid, 'Initial Δt [s]        : %.2f\n', params.dt);
fprintf(fid, 'Min Δt [s]            : %.2f\n', params.dt_min);
fprintf(fid, 'Max Δt [s]            : %.2f\n', params.dt_max);
fprintf(fid, 'Total Time [s]        : %.2f\n', params.Tmax);
fprintf(fid, 'Save Interval [s]     : %.2f\n', params.save_interval);
fprintf(fid, 'Max Newton Iterations : %d\n', params.max_iters);

fprintf(fid, '\n--- Boundary Conditions ---\n');
fprintf(fid, 'Top BC Type           : %s\n', params.top_bc_type);
fprintf(fid, 'Top BC Value          : %.3f\n', params.top_bc_value);
fprintf(fid, 'Bottom BC Type        : %s\n', params.bottom_bc_type);
fprintf(fid, 'Bottom BC Value       : %.3f\n', params.bottom_bc_value);

fprintf(fid, '\n--- Output Paths ---\n');
fprintf(fid, 'Base Folder           : %s\n', base_output_dir);
fprintf(fid, 'Data Folder           : %s\n', data_dir);
fprintf(fid, 'Figures Folder        : %s\n', figures_dir);

fprintf(fid, '\n--- Notes ---\n');
fprintf(fid, '(Model Developed by Marcus Nobrega, Ph.D. Feel free to contact me at marcusnobrega.engcivil@gmail.com.)\n');

fclose(fid);
fprintf('📝 Metadata saved to: %s\n', info_txt);
