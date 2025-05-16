%% =========================================================================
% 📘 EXAMPLE 1: Infiltration Experiment — Dirichlet–Free Drainage
% =========================================================================
clear; clc;
sim_name = 'Infiltration_Experiment';

% Define subdirectories using sim_name
base_output_dir   = fullfile('Modeling_Results', sim_name);
figures_dir       = fullfile(base_output_dir, 'Figures');
mesh_dir          = fullfile(figures_dir, 'Mesh');
data_dir         = fullfile(base_output_dir, 'Data');
profiles_dir      = fullfile(figures_dir, 'Profiles');
time_series_dir   = fullfile(figures_dir, 'TimeSeries');
diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
fdc_dir           = fullfile(figures_dir, 'FlowDuration');
wetting_dir       = fullfile(figures_dir, 'WettingFront');

% Create each directory only if it does not already exist
dirs_to_create = {base_output_dir, figures_dir, mesh_dir, profiles_dir, ...
                  time_series_dir, diagnostics_dir, data_dir, fdc_dir, wetting_dir};

for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end


% 🧱 DOMAIN SETUP
params.Nz = 40;
params.L = 1.0;
nonlin_factor = 1.5;
[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor,mesh_dir);

% ⏱️ TIME SETTINGS
params.Tmax = 3600;
params.dt = 1;
params.dt_min = 1;
params.dt_max = 10;
params.adapt_down = 0.5;
params.adapt_up = 1.2;
params.n_up = 3;
params.n_down = 10;
params.max_iters = 30;

% Saving/Plotting setup
params.save_interval_min = 10;
params.save_interval = max(params.save_interval_min * 60, params.dt);
params.save_steps = round(0 : params.save_interval / params.dt : params.Tmax / params.dt);
n_save = length(params.save_steps);
params.plot_interval = round(0.05 * length(params.save_steps));
params.Nt = round(params.Tmax / params.dt);

% SOIL PROFILE: Two-layer system
media_thicknesses = [0.5, 0.5];
media_props = struct( ...
    'alpha',   [2.0, 2.0], ...
    'n',       [1.5, 1.5], ...
    'theta_r', [0.04, 0.04], ...
    'theta_s', [0.42, 0.42], ...
    'S_s',     [1e-4, 1e-4], ...
    'Ks',      [5e-4, 5e-5]);

media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];
media_id = zeros(params.Nz, 1);
for i = 1:params.Nz
    zi = params.z(i);
    for j = 1:length(media_thicknesses)
        if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
            media_id(i) = j;
            break;
        elseif zi == media_interfaces(end)
            media_id(i) = length(media_thicknesses);
        end
    end
end

params.alpha   = media_props.alpha(media_id);
params.n       = media_props.n(media_id);
params.m       = 1 - 1 ./ params.n;
params.theta_r = media_props.theta_r(media_id);
params.theta_s = media_props.theta_s(media_id);
params.S_s     = media_props.S_s(media_id);
params.Ks      = media_props.Ks(media_id);

% BOUNDARY CONDITIONS
params.top_bc_type = "dirichlet";      % Constant head (ponding)
params.top_bc_value = 0.0;
params.bottom_bc_type = "free";        % Gravity drainage
params.bottom_bc_value = 0;

% SOURCES
params.surface_flux_series = zeros(params.Nt, 1);  % No rain
params.bottom_flux_series  = zeros(params.Nt, 1);  % No recharge
params.source_times = linspace(0, params.Tmax, params.Nt);
params.source_profile = zeros(params.Nz, params.Nt);

% INITIAL CONDITION: hydrostatic suction
h = -1 * ones(params.Nz, 1);
if params.bottom_bc_type == "dirichlet"
    h(1) = params.bottom_bc_value;
end
if params.top_bc_type == "dirichlet"
    h(end) = params.top_bc_value;
end

% === OUTPUT PREALLOCATION ====================================================
head_out         = nan(params.Nz, n_save);
theta_out        = nan(params.Nz, n_save);
flux_out         = nan(params.Nz + 1, n_save);
ponding_series   = zeros(1, n_save);
seepage_flux     = zeros(1, n_save);
top_flux         = zeros(1, n_save);
save_count       = 0;
mb_error_cumulative = 0;
ponding_depth = 0;
params.tol = 1e-7;

% Plotting time vector
n_plots = 20;
plot_times = linspace(0, params.Tmax, n_plots);
plot_index = 1;

% Pressure Limiter for Evaporation (Feddes)
params.h_lim_upper = -10; % [m]
params.h_lim_down = -0.1; % [m]

% === Set top_val and bottom_val based on BC type ============================
if params.top_bc_type == "dirichlet"
    top_val = params.top_bc_value;
else
    top_val = nan;
end

if params.bottom_bc_type == "dirichlet"
    bottom_val = params.bottom_bc_value;
else
    bottom_val = nan;
end

% === DRAINAGE CONFIGURATION =====================================
params.K_orifice    = zeros(params.Nz, 1);   % [m^(e_o)/s] Orifice coefficients (default: 0 [not activated])
params.exp_orifice    = 0.5 * ones(params.Nz, 1);    % [–] Exponent for each orifice (usually 0.5)

% Example: Place orifice at node node_idx with K = 0.0216 m^(0.5)/s, assuming a 10
% cm orifice
node_idx = [];
params.K_orifice(node_idx) = 0.0216 * pi() * (0.10) ^ 2 / 4 * sqrt(2 * 9.81); % Cd * Aef * sqrt(2 * g)
params.exp_orifice(node_idx) = 0.5;

% Spillway at top (only used if top BC is Neumann)
params.spillway_enabled   = true;
params.c_spillway         = 1.80 * 1.5;   % [m^(e)/s] -> Cd * Lef
params.h_spill            = 0.05;         % [m] (height above soil where spillway activates)
params.exp_spillway       = 1.5;          % [–]

Q_orifice_total = zeros(1, n_save);   % [m³/m²/s]
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s]


% === 📄 Save Simulation Log/Metadata =====================================
logfile = fullfile(base_output_dir, 'Log.txt');
fid = fopen(logfile, 'w');

fprintf(fid, '=====================================================\n');
fprintf(fid, '🔧 Mixed-Form Richards Model — Simulation Log\n');
fprintf(fid, '=====================================================\n\n');

% === Simulation Info
fprintf(fid, '🟢 Simulation Name     : %s\n', sim_name);
fprintf(fid, '📅 Date and Time       : %s\n', datestr(now));
fprintf(fid, '⏱️ Duration            : %.2f hr\n', params.Tmax / 3600);
fprintf(fid, '💾 Save Interval       : %.2f min\n', params.save_interval / 60);
fprintf(fid, '\n');

% === Mesh Info
fprintf(fid, '📐 Mesh Information\n');
fprintf(fid, '    • Number of Nodes      : %d\n', params.Nz);
fprintf(fid, '    • Total Depth          : %.2f m\n', abs(params.z(1)));
fprintf(fid, '    • Top Node Elevation   : %.2f m\n', params.z(end));
fprintf(fid, '    • Min Cell Thickness   : %.4f m\n', min(params.dz));
fprintf(fid, '    • Max Cell Thickness   : %.4f m\n', max(params.dz));
fprintf(fid, '\n');

% === Soil Properties
fprintf(fid, '🌱 Soil Hydraulic Properties (Van Genuchten)\n');
fprintf(fid, '    • α       : %.4f 1/m\n', mean(params.alpha));
fprintf(fid, '    • n       : %.4f [-]\n', mean(params.n));
fprintf(fid, '    • m       : %.4f [-]\n', mean(params.m));
fprintf(fid, '    • θ_s     : %.4f [m³/m³]\n', mean(params.theta_s));
fprintf(fid, '    • θ_r     : %.4f [m³/m³]\n', mean(params.theta_r));
fprintf(fid, '    • K_s     : %.2e [m/s]\n', mean(params.Ks));
fprintf(fid, '    • S_s     : %.2e [1/m]\n', mean(params.S_s));
fprintf(fid, '\n');

% === Boundary Conditions
fprintf(fid, '🔲 Boundary Conditions\n');
fprintf(fid, '    • Top BC Type     : %s\n', params.top_bc_type);
if isfield(params, 'top_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.top_bc_value);
end
fprintf(fid, '    • Bottom BC Type  : %s\n', params.bottom_bc_type);
if isfield(params, 'bottom_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.bottom_bc_value);
end
fprintf(fid, '\n');

% === Solver and Time-Stepping
fprintf(fid, '⚙️ Solver Configuration\n');
fprintf(fid, '    • Initial dt            : %.4f s\n', params.dt);
fprintf(fid, '    • Min/Max dt            : [%.4f, %.4f] s\n', params.dt_min, params.dt_max);
fprintf(fid, '    • Max Newton Iterations : %d\n', params.max_iters);
fprintf(fid, '    • Adapt Up/Down         : [%.2f, %.2f]\n', params.adapt_up, params.adapt_down);
fprintf(fid, '    • Convergence Threshold : %.1e\n', params.tol);
fprintf(fid, '\n');

% === Source/Sink Terms
if isfield(params, 'source_profile')
    total_source = sum(params.source_profile(:));
    fprintf(fid, '🛠️  Source Term Included  : Yes\n');
    fprintf(fid, '    • Total Integrated Source [m³/m²]: %.4f\n', total_source * sum(params.dz));
else
    fprintf(fid, '🛠️  Source Term Included  : No\n');
end
fprintf(fid, '\n');

% === Structural Drainage (Orifices and Spillways)
fprintf(fid, '🚰 Structural Drainage Sinks\n');
if isfield(params, 'K_orifice') && any(params.K_orifice > 0)
    idxs = find(params.K_orifice > 0);
    for i = 1:length(idxs)
        fprintf(fid, '    • Orifice at Node %d: K = %.2e, e = %.2f\n', ...
            idxs(i), params.K_orifice(idxs(i)), params.exp_orifice(idxs(i)));
    end
else
    fprintf(fid, '    • No orifices enabled.\n');
end

if isfield(params, 'spillway_enabled') && params.spillway_enabled
    fprintf(fid, '    • Spillway at Surface: c = %.2e, h_s = %.3f m, e = %.2f\n', ...
        params.c_spillway, params.h_spill, params.exp_spillway);
else
    fprintf(fid, '    • No spillway enabled.\n');
end
fprintf(fid, '\n');

fprintf(fid, '=====================================================\n');
fprintf(fid, '📍 End of Log\n');
fprintf(fid, '=====================================================\n');

fclose(fid);


% SAVE
save('Examples/Example1_Infiltration.mat')

%% =========================================================================
% 📘 EXAMPLE 2: Capillary Rise from Groundwater — Neumann–Dirichlet
% =========================================================================
clear; clc;
sim_name = 'Capillary_Rise_from_Groundwater';

% Define subdirectories using sim_name
base_output_dir   = fullfile('Modeling_Results', sim_name);
figures_dir       = fullfile(base_output_dir, 'Figures');
data_dir         = fullfile(base_output_dir, 'Data');
mesh_dir          = fullfile(figures_dir, 'Mesh');
profiles_dir      = fullfile(figures_dir, 'Profiles');
time_series_dir   = fullfile(figures_dir, 'TimeSeries');
diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
fdc_dir           = fullfile(figures_dir, 'FlowDuration');
wetting_dir       = fullfile(figures_dir, 'WettingFront');

% Create each directory only if it does not already exist
dirs_to_create = {base_output_dir, figures_dir, mesh_dir, profiles_dir, ...
                  time_series_dir, diagnostics_dir, fdc_dir, data_dir, wetting_dir};

for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end

% === DOMAIN SETUP =======================================================
params.Nz = 40;
params.L = 1.0;                         % Total depth [m]
nonlin_factor = 1.5;
[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);

% === TIME SETTINGS ======================================================
params.Tmax = 3600;                     % [s]
params.dt = 1;                          % [s]
params.dt_min = 1;
params.dt_max = 10;
params.adapt_down = 0.5;
params.adapt_up = 1.2;
params.n_up = 3;
params.n_down = 10;
params.max_iters = 30;

% Save and plot intervals
params.save_interval_min = 10;
params.save_interval = max(params.save_interval_min * 60, params.dt);
params.save_steps = round(0 : params.save_interval / params.dt : params.Tmax / params.dt);
n_save = length(params.save_steps);
params.plot_interval = round(0.05 * length(params.save_steps));
params.Nt = round(params.Tmax / params.dt);

% === SOIL PROFILE — 4-layer stratified column ===========================
media_thicknesses = [0.3, 0.2, 0.1, 0.4];  % [m]
media_props = struct( ...
    'alpha',   [2.0, 1.6, 1.2, 0.8], ...
    'n',       [1.5, 1.4, 1.35, 1.3], ...
    'theta_r', [0.04, 0.05, 0.06, 0.07], ...
    'theta_s', [0.42, 0.45, 0.47, 0.48], ...
    'S_s',     [1e-4, 5e-5, 3e-5, 1e-5], ...
    'Ks',      [5e-4, 2e-4, 1e-4, 5e-5]);

media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];
media_id = zeros(params.Nz, 1);
for i = 1:params.Nz
    zi = params.z(i);
    for j = 1:length(media_thicknesses)
        if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
            media_id(i) = j;
            break;
        elseif zi == media_interfaces(end)
            media_id(i) = length(media_thicknesses);
        end
    end
end

params.alpha   = media_props.alpha(media_id);
params.n       = media_props.n(media_id);
params.m       = 1 - 1 ./ params.n;
params.theta_r = media_props.theta_r(media_id);
params.theta_s = media_props.theta_s(media_id);
params.S_s     = media_props.S_s(media_id);
params.Ks      = media_props.Ks(media_id);

% === BOUNDARY CONDITIONS ================================================
params.top_bc_type     = "neumann";     % No surface flux (sealed)
params.top_bc_value    = 0.0;
params.bottom_bc_type  = "dirichlet";   % Water table at bottom
params.bottom_bc_value = 0.0;

% === SOURCES (none in this example) =====================================
params.surface_flux_series = zeros(params.Nt, 1);
params.bottom_flux_series  = zeros(params.Nt, 1);
params.source_times        = linspace(0, params.Tmax, params.Nt);
params.source_profile      = zeros(params.Nz, params.Nt);

% === INITIAL CONDITION: Hydrostatic =====================================
h = -1 * ones(params.Nz, 1);  % Slightly dry initial condition
if params.bottom_bc_type == "dirichlet"
    h(1) = params.bottom_bc_value;
end
if params.top_bc_type == "dirichlet"
    h(end) = params.top_bc_value;
end

% === OUTPUT PREALLOCATION ====================================================
head_out         = nan(params.Nz, n_save);
theta_out        = nan(params.Nz, n_save);
flux_out         = nan(params.Nz + 1, n_save);
ponding_series   = zeros(1, n_save);
seepage_flux     = zeros(1, n_save);
top_flux         = zeros(1, n_save);
save_count       = 0;
mb_error_cumulative = 0;
ponding_depth = 0;
params.tol = 1e-7;

% Plotting time vector
n_plots = 20;
plot_times = linspace(0, params.Tmax, n_plots);
plot_index = 1;

% Pressure Limiter for Evaporation (Feddes) 
params.h_lim_upper = -10; % [m]
params.h_lim_down = -0.1; % [m]

% === Set top_val and bottom_val based on BC type ============================
if params.top_bc_type == "dirichlet"
    top_val = params.top_bc_value;
else
    top_val = nan;
end

if params.bottom_bc_type == "dirichlet"
    bottom_val = params.bottom_bc_value;
else
    bottom_val = nan;
end

% === DRAINAGE CONFIGURATION =====================================
params.K_orifice    = zeros(params.Nz, 1);   % [m^(e_o)/s] Orifice coefficients (default: 0 [not activated])
params.exp_orifice    = 0.5 * ones(params.Nz, 1);    % [–] Exponent for each orifice (usually 0.5)

% Example: Place orifice at node node_idx with K = 0.0216 m^(0.5)/s, assuming a 10
% cm orifice
node_idx = 5;
params.K_orifice(node_idx) = 0.0216 * pi() * (0.10) ^ 2 / 4 * sqrt(2 * 9.81); % Cd * Aef * sqrt(2 * g)
params.exp_orifice(node_idx) = 0.5;

% Spillway at top (only used if top BC is Neumann)
params.spillway_enabled   = true;
params.c_spillway         = 1.80 * 1.5;   % [m^(e)/s] -> Cd * Lef
params.h_spill            = 0.05;         % [m] (height above soil where spillway activates)
params.exp_spillway       = 1.5;          % [–]

Q_orifice_total = zeros(1, n_save);   % [m³/m²/s]
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s]


% === 📄 Save Simulation Log/Metadata =====================================
logfile = fullfile(base_output_dir, 'Log.txt');
fid = fopen(logfile, 'w');

fprintf(fid, '=====================================================\n');
fprintf(fid, '🔧 Mixed-Form Richards Model — Simulation Log\n');
fprintf(fid, '=====================================================\n\n');

% === Simulation Info
fprintf(fid, '🟢 Simulation Name     : %s\n', sim_name);
fprintf(fid, '📅 Date and Time       : %s\n', datestr(now));
fprintf(fid, '⏱️ Duration            : %.2f hr\n', params.Tmax / 3600);
fprintf(fid, '💾 Save Interval       : %.2f min\n', params.save_interval / 60);
fprintf(fid, '\n');

% === Mesh Info
fprintf(fid, '📐 Mesh Information\n');
fprintf(fid, '    • Number of Nodes      : %d\n', params.Nz);
fprintf(fid, '    • Total Depth          : %.2f m\n', abs(params.z(1)));
fprintf(fid, '    • Top Node Elevation   : %.2f m\n', params.z(end));
fprintf(fid, '    • Min Cell Thickness   : %.4f m\n', min(params.dz));
fprintf(fid, '    • Max Cell Thickness   : %.4f m\n', max(params.dz));
fprintf(fid, '\n');

% === Soil Properties
fprintf(fid, '🌱 Soil Hydraulic Properties (Van Genuchten)\n');
fprintf(fid, '    • α       : %.4f 1/m\n', mean(params.alpha));
fprintf(fid, '    • n       : %.4f [-]\n', mean(params.n));
fprintf(fid, '    • m       : %.4f [-]\n', mean(params.m));
fprintf(fid, '    • θ_s     : %.4f [m³/m³]\n', mean(params.theta_s));
fprintf(fid, '    • θ_r     : %.4f [m³/m³]\n', mean(params.theta_r));
fprintf(fid, '    • K_s     : %.2e [m/s]\n', mean(params.Ks));
fprintf(fid, '    • S_s     : %.2e [1/m]\n', mean(params.S_s));
fprintf(fid, '\n');

% === Boundary Conditions
fprintf(fid, '🔲 Boundary Conditions\n');
fprintf(fid, '    • Top BC Type     : %s\n', params.top_bc_type);
if isfield(params, 'top_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.top_bc_value);
end
fprintf(fid, '    • Bottom BC Type  : %s\n', params.bottom_bc_type);
if isfield(params, 'bottom_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.bottom_bc_value);
end
fprintf(fid, '\n');

% === Solver and Time-Stepping
fprintf(fid, '⚙️ Solver Configuration\n');
fprintf(fid, '    • Initial dt            : %.4f s\n', params.dt);
fprintf(fid, '    • Min/Max dt            : [%.4f, %.4f] s\n', params.dt_min, params.dt_max);
fprintf(fid, '    • Max Newton Iterations : %d\n', params.max_iters);
fprintf(fid, '    • Adapt Up/Down         : [%.2f, %.2f]\n', params.adapt_up, params.adapt_down);
fprintf(fid, '    • Convergence Threshold : %.1e\n', params.tol);
fprintf(fid, '\n');

% === Source/Sink Terms
if isfield(params, 'source_profile')
    total_source = sum(params.source_profile(:));
    fprintf(fid, '🛠️  Source Term Included  : Yes\n');
    fprintf(fid, '    • Total Integrated Source [m³/m²]: %.4f\n', total_source * sum(params.dz));
else
    fprintf(fid, '🛠️  Source Term Included  : No\n');
end
fprintf(fid, '\n');

% === Structural Drainage (Orifices and Spillways)
fprintf(fid, '🚰 Structural Drainage Sinks\n');
if isfield(params, 'K_orifice') && any(params.K_orifice > 0)
    idxs = find(params.K_orifice > 0);
    for i = 1:length(idxs)
        fprintf(fid, '    • Orifice at Node %d: K = %.2e, e = %.2f\n', ...
            idxs(i), params.K_orifice(idxs(i)), params.exp_orifice(idxs(i)));
    end
else
    fprintf(fid, '    • No orifices enabled.\n');
end

if isfield(params, 'spillway_enabled') && params.spillway_enabled
    fprintf(fid, '    • Spillway at Surface: c = %.2e, h_s = %.3f m, e = %.2f\n', ...
        params.c_spillway, params.h_spill, params.exp_spillway);
else
    fprintf(fid, '    • No spillway enabled.\n');
end
fprintf(fid, '\n');

fprintf(fid, '=====================================================\n');
fprintf(fid, '📍 End of Log\n');
fprintf(fid, '=====================================================\n');

fclose(fid);


% === SAVE =================================================================
save('Examples/Example2_CapillaryRise.mat')

%% =========================================================================
% 📘 EXAMPLE 3: Celia et al. (1990) Benchmark — Ponded Infiltration
% =========================================================================
sim_name = 'Celia_1990';  % 📝 CHANGE THIS for each run


% Define subdirectories using sim_name
base_output_dir   = fullfile('Modeling_Results', sim_name);
figures_dir       = fullfile(base_output_dir, 'Figures');
data_dir         = fullfile(base_output_dir, 'Data');
mesh_dir          = fullfile(figures_dir, 'Mesh');
profiles_dir      = fullfile(figures_dir, 'Profiles');
time_series_dir   = fullfile(figures_dir, 'TimeSeries');
diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
fdc_dir           = fullfile(figures_dir, 'FlowDuration');
wetting_dir       = fullfile(figures_dir, 'WettingFront');

% Create each directory only if it does not already exist
dirs_to_create = {base_output_dir, figures_dir, mesh_dir, profiles_dir, ...
                  time_series_dir, diagnostics_dir, data_dir, fdc_dir, wetting_dir};

for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end


% === 1. DOMAIN & MESH DISCRETIZATION =====================================

params.Nz = 41;               % Number of vertical nodes [-]
params.L  = 1.0;              % Total soil depth [m]
nonlin_factor = 1;            % Grid refinement factor (1 = uniform)
params.LID_area = 1;       % 1D column area [m2]

% Generate refined mesh (Hydrus-style, refined near surface)
[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);

% === 2. TIME DISCRETIZATION ==============================================

params.Tmax = 24*3600;        % Total simulation time [s]
params.dt   = 5*60;           % Initial time step [s]
params.dt_min = 0.00001;           % Minimum dt [s]
params.dt_max = 15*60;        % Maximum dt [s]

% Adaptive timestep control
params.adapt_down = 0.5;      % Shrink factor
params.adapt_up   = 2.0;      % Growth factor
params.n_up       = 25;       % Threshold for fast convergence
params.n_down     = 75;       % Threshold for slow convergence

params.Nt = round(params.Tmax / params.dt);            % Max steps
params.max_iters = 100;                                % Newton max iters

% Output saving frequency
params.save_interval_min = 15;                         % [min]
params.save_interval = max(params.save_interval_min * 60, params.dt); % [s]
params.save_steps = round(0:params.save_interval / params.dt : params.Tmax / params.dt);
n_save = length(params.save_steps);

% Preallocate outputs
head_out         = nan(params.Nz, n_save);
theta_out        = nan(params.Nz, n_save);
flux_out         = nan(params.Nz + 1, n_save);
ponding_series   = zeros(1, n_save);
seepage_flux     = zeros(1, n_save);
top_flux         = zeros(1, n_save);
save_count       = 0;
mb_error_cumulative = 0;

% Plotting time vector
n_plots = 20;
plot_times = linspace(0, params.Tmax, n_plots);
plot_index = 1;

% Pressure Limiter for Evaporation (Feddes)
params.h_lim_upper = -10; % [m]
params.h_lim_down = -0.1; % [m]

% === 3. MULTILAYER SOIL HYDRAULIC PROPERTIES =============================

% Define soil layer thicknesses [m], from bottom to top
media_thicknesses = [1.0];
media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];  % Include bottom
n_layers = length(media_thicknesses);

% Van Genuchten + Ks parameters (per layer)
media_props = struct( ...
    'alpha',   [3.55], ...
    'n',       [2.0], ...
    'theta_r', [0.102], ...
    'theta_s', [0.368], ...
    'S_s',     [0], ...
    'Ks',      [9.22e-5] ...
    );

% Assign media index to each node
media_id = zeros(1 , params.Nz);
for i = 1:params.Nz
    zi = params.z(i);
    for j = 1:n_layers
        if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
            media_id(i) = j;
            break;
        elseif zi == media_interfaces(end)
            media_id(i) = n_layers;  % Top edge case
        end
    end
end

media_props.labels = {'Celia_Soil'};

% Assign per-node hydraulic parameters
params.alpha   = media_props.alpha(media_id);
params.n       = media_props.n(media_id);
params.m       = 1 - 1 ./ params.n;
params.theta_r = media_props.theta_r(media_id);
params.theta_s = media_props.theta_s(media_id);
params.S_s     = media_props.S_s(media_id);
params.Ks      = media_props.Ks(media_id);

% Water Retention Curves
plot_vg_retention_curves(media_props, params, figures_dir)

% === 4. SOLVER ITERATION ERROR TOLERANCE  ==============================================
% Newtown Raphson + Line Search Error Tolerance
params.tol = 1e-6; % [m]

% === 5. BOUNDARY CONDITIONS ==============================================

% === Top Boundary Condition Type =========================================
% Options: 'dirichlet' = fixed pressure head; 'neumann' = flux
params.top_bc_type  = "dirichlet";  % Options: 'dirichlet' or 'neumann'
params.top_bc_value = -0.75;         % Pressure head at top [m] (set NaN for Neumann)

% === Bottom Boundary Condition Type ======================================
% Options: 'dirichlet', 'neumann', 'free', or 'noflow'
params.bottom_bc_type  = "dirichlet";
params.bottom_bc_value = -10;         % [m] (ignored for 'free' and 'noflow')

% === Surface Flux (Top) — from Hydrologic Catchment Model ===============
% Catchment_Outputs returns:
%   • surface_flux_time [s]
%   • surface_flux_vals [m/s]
%   • C_top (optional water quality tracer) [mg / L]

% Example of time-varying surface fluxes and bottom fluxes, uncomment if
% you want to manually define the surface and bottom fluxes
% dt_input = 60;
% Nt_input = round(params.Tmax / dt_input);
% t_input  = (0:Nt_input - 1) * dt_input;
% surface_flux_vals   = -1e-4 * sin(linspace(0, pi, Nt_input))';   % Example rainfall or use a constant value
% bottom_flux_vals = 0 * surface_flux_vals;

if params.top_bc_type == "neumann"
    [surface_flux_time, surface_flux_vals, C_top] = Catchment_Outputs(params.dt , params.LID_area);  
    params.surface_flux_time = surface_flux_time;      % Time vector for interpolation [s]
    params.surface_flux_vals = surface_flux_vals;      % Time-varying flux [m/s]
else
    dt_input = 60; Nt_input = round(params.Tmax / dt_input); t_input  = (0:Nt_input - 1) * dt_input;
    params.surface_flux_time = Nt_input; params.surface_flux_vals = 0 * sin(linspace(0, pi, Nt_input)); % Example rainfall or use a constant value
end

% === Bottom Flux — (e.g., Recharge, Lateral Leakage, etc.) ==============
% You can assign another hydrologic model or use constant zero
params.bottom_flux_time  = params.surface_flux_time;
params.bottom_flux_vals  = 0 * params.surface_flux_vals;   % No recharge in this example

% === STRUCTURAL DRAINAGE SINKS ==========================================

% === Orifice Flow Parameters =============================================
% General form: Q_orifice = K_orifice * (max(h, 0)) ^ exp_orifice
params.K_orifice     = zeros(1, params.Nz);             % [m^(e)/s] Coefficient per node
params.exp_orifice   = 0.5 * ones(1, params.Nz);         % [-] Exponent (typically 0.5 for orifices)

% Example: Add orifice at node 1 (10 cm diameter)
node_idx = 20;                  % Vector indicating orifice nodes
n_orifices = 0;                 % Vector indicating the number of orifices in a node
Cd       = 0.6;                 % Discharge coefficient [-] (0.5 - 0.6)
D        = 0.10;                % Vector of orifice diameters [m]
Aeff     = pi * D .^2 / 4;      % Effective area [m²]
g        = 9.81;                % Gravity [m/s²]
params.K_orifice(node_idx) = n_orifices * Cd .* Aeff * sqrt(2 * g);   % Full orifice flow coefficient

% === Spillway Flow at Top Node (for Neumann BC only) =====================
% General form: Q_spillway = c_spillway * (h - h_spill) ^ exp_spillway, when h > h_spill
params.spillway_enabled   = true;              % Enable spillway only for Neumann BC
params.c_spillway         = 0 *1.8 * 1.5;      % Cd * L [m^(e)/s] for weir-type equation
params.h_spill            = 0.05;              % Activation height above soil [m]
params.exp_spillway       = 1.5;               % Exponent for free surface overflow

% === Output Arrays to Track Flow =========================================
Q_orifice_total  = zeros(1, n_save);  % [m³/m²/s] per time step
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s] per time step

% === B.C and Ponding Initial Values ======================================
if params.top_bc_type == "dirichlet"
    ponding_depth = max(params.top_bc_value, 0);
    ponding_prev  = ponding_depth;
    top_val       = params.top_bc_value;
else
    ponding_depth = 0;
    ponding_prev  = 0;
    top_val       = 0;
end

% === 6. SOURCE TERM AND INITIAL CONDITIONS ===============================

% Create source profile [Nz x Nt] – default: zero
params.source_times   = linspace(0, params.Tmax, params.Nt);
params.source_profile = zeros(params.Nz, params.Nt);  % Fill as needed

% Initial condition: hydrostatic profile
p = -10;                      % Uniform suction [m]
h = p * ones(1, params.Nz);   % Default

if params.bottom_bc_type == "dirichlet"
    h(1) = params.bottom_bc_value;
end
if params.top_bc_type == "dirichlet"
    h(end) = params.top_bc_value;
end

% === 📄 7. Save Simulation Log/Metadata =====================================
logfile = fullfile(base_output_dir, 'Log.txt');
fid = fopen(logfile, 'w');

fprintf(fid, '=====================================================\n');
fprintf(fid, '🔧 Mixed-Form Richards Model — Simulation Log\n');
fprintf(fid, '=====================================================\n\n');

% === Simulation Info
fprintf(fid, '🟢 Simulation Name     : %s\n', sim_name);
fprintf(fid, '📅 Date and Time       : %s\n', datestr(now));
fprintf(fid, '⏱️ Duration            : %.2f hr\n', params.Tmax / 3600);
fprintf(fid, '💾 Save Interval       : %.2f min\n', params.save_interval / 60);
fprintf(fid, '\n');

% === Mesh Info
fprintf(fid, '📐 Mesh Information\n');
fprintf(fid, '    • Number of Nodes      : %d\n', params.Nz);
fprintf(fid, '    • Total Depth          : %.2f m\n', abs(params.z(1)));
fprintf(fid, '    • Top Node Elevation   : %.2f m\n', params.z(end));
fprintf(fid, '    • Min Cell Thickness   : %.4f m\n', min(params.dz));
fprintf(fid, '    • Max Cell Thickness   : %.4f m\n', max(params.dz));
fprintf(fid, '\n');

% === Soil Properties
fprintf(fid, '🌱 Soil Hydraulic Properties (Van Genuchten)\n');
fprintf(fid, '    • α       : %.4f 1/m\n', mean(params.alpha));
fprintf(fid, '    • n       : %.4f [-]\n', mean(params.n));
fprintf(fid, '    • m       : %.4f [-]\n', mean(params.m));
fprintf(fid, '    • θ_s     : %.4f [m³/m³]\n', mean(params.theta_s));
fprintf(fid, '    • θ_r     : %.4f [m³/m³]\n', mean(params.theta_r));
fprintf(fid, '    • K_s     : %.2e [m/s]\n', mean(params.Ks));
fprintf(fid, '    • S_s     : %.2e [1/m]\n', mean(params.S_s));
fprintf(fid, '\n');

% === Boundary Conditions
fprintf(fid, '🔲 Boundary Conditions\n');
fprintf(fid, '    • Top BC Type     : %s\n', params.top_bc_type);
if isfield(params, 'top_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.top_bc_value);
end
fprintf(fid, '    • Bottom BC Type  : %s\n', params.bottom_bc_type);
if isfield(params, 'bottom_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.bottom_bc_value);
end
fprintf(fid, '\n');

% === Drainage Sinks Information ==========================================
fprintf(fid, '💧 Drainage Sinks Information\n');

% Orifice information
if any(params.K_orifice > 0)
    fprintf(fid, '    • Orifices Enabled     : %s\n', "Yes");
else
    fprintf(fid, '    • Orifices Enabled     : %s\n', "No");
end
if any(params.K_orifice > 0)
    fprintf(fid, '    • Orifice Exponent(s)  : %.2f (for each node)\n', mean(params.exp_orifice));
    fprintf(fid, '    • Orifice Coefficients : %.4f [m^(exp_orifice)/s]\n', mean(params.K_orifice(params.K_orifice > 0)));
end

% Spillway information
if params.spillway_enabled
    fprintf(fid, '    • Spillway Enabled     : %s\n', "Yes");
else
    fprintf(fid, '    • Spillway Enabled     : %s\n', "No");
end
if params.spillway_enabled
    fprintf(fid, '    • Spillway Coefficient : %.4f [m^(exp_spillway)/s]\n', params.c_spillway);
    fprintf(fid, '    • Spillway Exponent    : %.2f\n', params.exp_spillway);
    fprintf(fid, '    • Spillway Height      : %.2f [m]\n', params.h_spill);
end

fprintf(fid, '\n');

% === Solver and Time-Stepping
fprintf(fid, '⚙️ Solver Configuration\n');
fprintf(fid, '    • Initial dt            : %.4f s\n', params.dt);
fprintf(fid, '    • Min/Max dt            : [%.4f, %.4f] s\n', params.dt_min, params.dt_max);
fprintf(fid, '    • Max Newton Iterations : %d\n', params.max_iters);
fprintf(fid, '    • Adapt Up/Down         : [%.2f, %.2f]\n', params.adapt_up, params.adapt_down);
fprintf(fid, '    • Convergence Threshold : %.1e\n', params.tol);
fprintf(fid, '\n');

% === Source/Sink Terms
if isfield(params, 'source_profile')
    total_source = sum(params.source_profile(:));
    fprintf(fid, '🛠️  Source Term Included  : Yes\n');
    fprintf(fid, '    • Total Integrated Source [m³/m²]: %.4f\n', total_source * sum(params.dz));
else
    fprintf(fid, '🛠️  Source Term Included  : No\n');
end
fprintf(fid, '\n');

fprintf(fid, '=====================================================\n');
fprintf(fid, '📍 End of Log\n');
fprintf(fid, '=====================================================\n');

fclose(fid);

% === SAVE FILE ===========================================================
save('Examples/Celia1990.mat')

%% =========================================================================
% 📘 EXAMPLE 4: Permeable Pavement — Multi-layer Profile (San Antonio)
% =========================================================================
clear; clc;

% === 🗂️ Simulation Name and Directory Setup =============================
sim_name = 'PermeablePavement_SanAntonio';

% Define subdirectories using sim_name
base_output_dir   = fullfile('Modeling_Results', sim_name);
figures_dir       = fullfile(base_output_dir, 'Figures');
data_dir         = fullfile(base_output_dir, 'Data');
mesh_dir          = fullfile(figures_dir, 'Mesh');
profiles_dir      = fullfile(figures_dir, 'Profiles');
time_series_dir   = fullfile(figures_dir, 'TimeSeries');
diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
fdc_dir           = fullfile(figures_dir, 'FlowDuration');
wetting_dir       = fullfile(figures_dir, 'WettingFront');

% Create directories if not already existing
dirs_to_create = {base_output_dir, figures_dir, mesh_dir, ...
                  profiles_dir, time_series_dir, diagnostics_dir, ...
                  fdc_dir, data_dir, wetting_dir};

for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end

% === 1. DOMAIN & MESH ===================================================
params.Nz = 21;
params.L  = 0.203 + 0.102 + 0.0254;  % Total pavement profile depth [m]
nonlin_factor = 1.3;

[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);

% === 2. TIME SETTINGS ===================================================
params.Tmax        = 24 * 3600;   % [s] = 24 hours
params.dt          = 12 * 60;     % [s]
params.dt_min      = params.dt;
params.dt_max      = params.dt;
params.adapt_down  = 0.5;
params.adapt_up    = 1.2;
params.n_up        = 3;
params.n_down      = 10;
params.max_iters   = 1000;
params.Nt          = round(params.Tmax / params.dt);

params.save_interval_min = 1;
params.save_interval     = max(params.save_interval_min * 60, params.dt);
params.save_steps        = round(0 : params.save_interval / params.dt : params.Tmax / params.dt);
n_save                   = length(params.save_steps);

params.plot_interval = round(0.05 * params.Nt);

% === 3. MULTILAYER SOIL PROPERTIES ======================================
media_thicknesses = [0.203, 0.102, 0.0254];  % Bottom to top

media_props = struct( ...
    'alpha',   [20, 10, 0.5], ...
    'n',       [2.8, 2.5, 1.5], ...
    'theta_r', [0.02, 0.02, 0.03], ...
    'theta_s', [0.30, 0.25, 0.20], ...
    'S_s',     [1e-4, 5e-5, 2e-5], ...
    'Ks',      [1e-2, 1.5e-2, 2e-3]);

media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];
n_layers = length(media_thicknesses);

media_id = zeros(params.Nz, 1);
for i = 1:params.Nz
    zi = params.z(i);
    for j = 1:n_layers
        if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
            media_id(i) = j;
            break;
        elseif zi == media_interfaces(end)
            media_id(i) = n_layers;
        end
    end
end

params.alpha   = media_props.alpha(media_id);
params.n       = media_props.n(media_id);
params.m       = 1 - 1 ./ params.n;
params.theta_r = media_props.theta_r(media_id);
params.theta_s = media_props.theta_s(media_id);
params.S_s     = media_props.S_s(media_id);
params.Ks      = media_props.Ks(media_id);

% === 4. BOUNDARY CONDITIONS =============================================
params.top_bc_type     = "dirichlet";
params.top_bc_value    = 0.2;  % Positive ponding at surface [m]
params.bottom_bc_type  = "free";
params.bottom_bc_value = -10;  % Not used for free

if params.top_bc_type == "dirichlet"
    ponding_depth = max(params.top_bc_value, 0);
    ponding_prev  = ponding_depth;
    top_val       = params.top_bc_value;
else
    ponding_depth = 0;
    ponding_prev  = 0;
    top_val       = 0;
end

% === 5. SOURCE TERMS & INITIAL CONDITIONS ================================
params.source_times   = linspace(0, params.Tmax, params.Nt);
params.source_profile = zeros(params.Nz, params.Nt);

dt_input = 60;
Nt_input = round(params.Tmax / dt_input);
t_input  = (0:Nt_input - 1) * dt_input;

params.surface_flux_time = t_input;
params.surface_flux_vals = zeros(Nt_input, 1);
params.bottom_flux_time  = t_input;
params.bottom_flux_vals  = zeros(Nt_input, 1);

h = -10 * ones(params.Nz, 1);  % Initial pressure head

if params.top_bc_type == "dirichlet"
    h(end) = params.top_bc_value;
end
if params.bottom_bc_type == "dirichlet"
    h(1) = params.bottom_bc_value;
end


% === OUTPUT PREALLOCATION ====================================================
head_out         = nan(params.Nz, n_save);
theta_out        = nan(params.Nz, n_save);
flux_out         = nan(params.Nz + 1, n_save);
ponding_series   = zeros(1, n_save);
seepage_flux     = zeros(1, n_save);
top_flux         = zeros(1, n_save);
save_count       = 0;
mb_error_cumulative = 0;
ponding_depth = 0;
params.tol = 1e-7;

% Plotting time vector
n_plots = 20;
plot_times = linspace(0, params.Tmax, n_plots);
plot_index = 1;

% Pressure Limiter for Evaporation (Feddes)
params.h_lim_upper = -10; % [m]
params.h_lim_down = -0.1; % [m]

% === Set top_val and bottom_val based on BC type ============================
if params.top_bc_type == "dirichlet"
    top_val = params.top_bc_value;
else
    top_val = nan;
end

if params.bottom_bc_type == "dirichlet"
    bottom_val = params.bottom_bc_value;
else
    bottom_val = nan;
end

% === DRAINAGE CONFIGURATION =====================================
params.K_orifice    = zeros(params.Nz, 1);   % [m^(e_o)/s] Orifice coefficients (default: 0 [not activated])
params.exp_orifice    = 0.5 * ones(params.Nz, 1);    % [–] Exponent for each orifice (usually 0.5)

% Example: Place orifice at node node_idx with K = 0.0216 m^(0.5)/s, assuming a 10
% cm orifice
node_idx = [];
params.K_orifice(node_idx) = 0.0216 * pi() * (0.10) ^ 2 / 4 * sqrt(2 * 9.81); % Cd * Aef * sqrt(2 * g)
params.exp_orifice(node_idx) = 0.5;

% Spillway at top (only used if top BC is Neumann)
params.spillway_enabled   = true;
params.c_spillway         = 1.80 * 1.5;   % [m^(e)/s] -> Cd * Lef
params.h_spill            = 0.05;         % [m] (height above soil where spillway activates)
params.exp_spillway       = 1.5;          % [–]

Q_orifice_total = zeros(1, n_save);   % [m³/m²/s]
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s]


% === 📄 Save Simulation Log/Metadata =====================================
logfile = fullfile(base_output_dir, 'Log.txt');
fid = fopen(logfile, 'w');

fprintf(fid, '=====================================================\n');
fprintf(fid, '🔧 Mixed-Form Richards Model — Simulation Log\n');
fprintf(fid, '=====================================================\n\n');

% === Simulation Info
fprintf(fid, '🟢 Simulation Name     : %s\n', sim_name);
fprintf(fid, '📅 Date and Time       : %s\n', datestr(now));
fprintf(fid, '⏱️ Duration            : %.2f hr\n', params.Tmax / 3600);
fprintf(fid, '💾 Save Interval       : %.2f min\n', params.save_interval / 60);
fprintf(fid, '\n');

% === Mesh Info
fprintf(fid, '📐 Mesh Information\n');
fprintf(fid, '    • Number of Nodes      : %d\n', params.Nz);
fprintf(fid, '    • Total Depth          : %.2f m\n', abs(params.z(1)));
fprintf(fid, '    • Top Node Elevation   : %.2f m\n', params.z(end));
fprintf(fid, '    • Min Cell Thickness   : %.4f m\n', min(params.dz));
fprintf(fid, '    • Max Cell Thickness   : %.4f m\n', max(params.dz));
fprintf(fid, '\n');

% === Soil Properties
fprintf(fid, '🌱 Soil Hydraulic Properties (Van Genuchten)\n');
fprintf(fid, '    • α       : %.4f 1/m\n', mean(params.alpha));
fprintf(fid, '    • n       : %.4f [-]\n', mean(params.n));
fprintf(fid, '    • m       : %.4f [-]\n', mean(params.m));
fprintf(fid, '    • θ_s     : %.4f [m³/m³]\n', mean(params.theta_s));
fprintf(fid, '    • θ_r     : %.4f [m³/m³]\n', mean(params.theta_r));
fprintf(fid, '    • K_s     : %.2e [m/s]\n', mean(params.Ks));
fprintf(fid, '    • S_s     : %.2e [1/m]\n', mean(params.S_s));
fprintf(fid, '\n');

% === Boundary Conditions
fprintf(fid, '🔲 Boundary Conditions\n');
fprintf(fid, '    • Top BC Type     : %s\n', params.top_bc_type);
if isfield(params, 'top_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.top_bc_value);
end
fprintf(fid, '    • Bottom BC Type  : %s\n', params.bottom_bc_type);
if isfield(params, 'bottom_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.bottom_bc_value);
end
fprintf(fid, '\n');

% === Solver and Time-Stepping
fprintf(fid, '⚙️ Solver Configuration\n');
fprintf(fid, '    • Initial dt            : %.4f s\n', params.dt);
fprintf(fid, '    • Min/Max dt            : [%.4f, %.4f] s\n', params.dt_min, params.dt_max);
fprintf(fid, '    • Max Newton Iterations : %d\n', params.max_iters);
fprintf(fid, '    • Adapt Up/Down         : [%.2f, %.2f]\n', params.adapt_up, params.adapt_down);
fprintf(fid, '    • Convergence Threshold : %.1e\n', params.tol);
fprintf(fid, '\n');

% === Source/Sink Terms
if isfield(params, 'source_profile')
    total_source = sum(params.source_profile(:));
    fprintf(fid, '🛠️  Source Term Included  : Yes\n');
    fprintf(fid, '    • Total Integrated Source [m³/m²]: %.4f\n', total_source * sum(params.dz));
else
    fprintf(fid, '🛠️  Source Term Included  : No\n');
end
fprintf(fid, '\n');

% === Structural Drainage (Orifices and Spillways)
fprintf(fid, '🚰 Structural Drainage Sinks\n');
if isfield(params, 'K_orifice') && any(params.K_orifice > 0)
    idxs = find(params.K_orifice > 0);
    for i = 1:length(idxs)
        fprintf(fid, '    • Orifice at Node %d: K = %.2e, e = %.2f\n', ...
            idxs(i), params.K_orifice(idxs(i)), params.exp_orifice(idxs(i)));
    end
else
    fprintf(fid, '    • No orifices enabled.\n');
end

if isfield(params, 'spillway_enabled') && params.spillway_enabled
    fprintf(fid, '    • Spillway at Surface: c = %.2e, h_s = %.3f m, e = %.2f\n', ...
        params.c_spillway, params.h_spill, params.exp_spillway);
else
    fprintf(fid, '    • No spillway enabled.\n');
end
fprintf(fid, '\n');

fprintf(fid, '=====================================================\n');
fprintf(fid, '📍 End of Log\n');
fprintf(fid, '=====================================================\n');

fclose(fid);


% === 6. Save Setup to File ==============================================
save(fullfile('Examples', [sim_name, '.mat']))

close all

%% =========================================================================
% 📘 EXAMPLE 4: Permeable Pavement — Multi-layer Profile (San Antonio)
% =========================================================================
clear; clc;

% === 🗂️ Simulation Name and Directory Setup =============================
sim_name = 'Ev6_May13_0110';

% Define subdirectories using sim_name
base_output_dir   = fullfile('Modeling_Results', sim_name);
figures_dir       = fullfile(base_output_dir, 'Figures');
data_dir         = fullfile(base_output_dir, 'Data');
mesh_dir          = fullfile(figures_dir, 'Mesh');
profiles_dir      = fullfile(figures_dir, 'Profiles');
time_series_dir   = fullfile(figures_dir, 'TimeSeries');
diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
fdc_dir           = fullfile(figures_dir, 'FlowDuration');
wetting_dir       = fullfile(figures_dir, 'WettingFront');

% Create each directory only if it does not already exist
dirs_to_create = {base_output_dir, figures_dir, mesh_dir, profiles_dir, ...
                  time_series_dir, diagnostics_dir, data_dir, fdc_dir, wetting_dir};

for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end


% === 1. DOMAIN & MESH DISCRETIZATION =====================================

params.Nz = 27;                      % Number of vertical nodes [-]
params.L  = 0.203 + 0.102 + 0.0254;  % Total pavement profile depth [m]
nonlin_factor = 1.1;                 % Grid refinement factor (1 = uniform)
params.LID_area = 1;                 % 1D column area [m2]


[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);
% Generate refined mesh (Hydrus-style, refined near surface)

% === 2. TIME DISCRETIZATION ==============================================

params.Tmax = 682160*60;        % Total simulation time [s]
params.dt   = 5*60;           % Initial time step [s]
params.dt_min = 0.001;           % Minimum dt [s]
params.dt_max = 5*60;       % Maximum dt [s]

% Adaptive timestep control
params.adapt_down = 0.5;      % Shrink factor
params.adapt_up   = 2.0;      % Growth factor
params.n_up       = 5;       % Threshold for fast convergence
params.n_down     = 10;       % Threshold for slow convergence

params.Nt = round(params.Tmax / params.dt);           % Max steps
params.max_iters = 20;                                % Newton max iters

% Output saving frequency
params.save_interval_min = 5;                         % [min]
params.save_interval = max(params.save_interval_min * 60, params.dt); % [s]
params.save_steps = round(0:params.save_interval / params.dt : params.Tmax / params.dt);
n_save = length(params.save_steps);

% Preallocate outputs
head_out         = nan(params.Nz, n_save);
theta_out        = nan(params.Nz, n_save);
flux_out         = nan(params.Nz + 1, n_save);
ponding_series   = zeros(1, n_save);
seepage_flux     = zeros(1, n_save);
top_flux         = zeros(1, n_save);
save_count       = 0;
mb_error_cumulative = 0;

% Plotting time vector
n_plots = 400;
plot_times = linspace(0, params.Tmax, n_plots);
plot_index = 1;

% Pressure Limiter for Evaporation (Feddes)
params.h_lim_upper = -2; % [m]
params.h_lim_down = -0.1; % [m]

% === 3. MULTILAYER SOIL PROPERTIES ======================================
media_thicknesses = [0.203, 0.102, 0.0254];  % Bottom to top

media_props = struct( ...
    'alpha',   [7.5, 3.5, 7.5], ...        % van Genuchten alpha (1/m)
    'n',       [2.5, 1.75, 2.5], ...       % van Genuchten n
    'theta_r', [0.005, 0.005, 0.005], ...  % Residual water content
    'theta_s', [0.42, 0.42, 0.42], ...     % Saturated water content
    'S_s',     [5e-5, 1e-4, 1e-4], ...     % Specific storage (1/m)
    'Ks',      [2e-2, 2e-2, 2e-2]);        % Saturated hydraulic conductivity (m/s)

media_props.labels = {'ASTM $\#8$', ...
                      'ASTM $\#2$', ...
                      'Permeable Asphalt'};

media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];
n_layers = length(media_thicknesses);

media_id = zeros(params.Nz, 1);
for i = 1:params.Nz
    zi = params.z(i);
    for j = 1:n_layers
        if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
            media_id(i) = j;
            break;
        elseif zi == media_interfaces(end)
            media_id(i) = n_layers;
        end
    end
end

params.alpha   = media_props.alpha(media_id);
params.n       = media_props.n(media_id);
params.m       = 1 - 1 ./ params.n;
params.theta_r = media_props.theta_r(media_id);
params.theta_s = media_props.theta_s(media_id);
params.S_s     = media_props.S_s(media_id);
params.Ks      = media_props.Ks(media_id);

% Water Retention Curves
plot_vg_retention_curves(media_props, params, figures_dir)

% === 4. SOLVER ITERATION ERROR TOLERANCE  ==============================================
% Newtown Raphson + Line Search Error Tolerance
params.tol = 1e-3; % [m]

% === 5. BOUNDARY CONDITIONS ==============================================

% === Top Boundary Condition Type =========================================
% Options: 'dirichlet' = fixed pressure head; 'neumann' = flux
params.top_bc_type  = "neumann";  % Options: 'dirichlet' or 'neumann'
params.top_bc_value = 0.0;         % Pressure head at top [m] (set NaN for Neumann)

% === Bottom Boundary Condition Type ======================================
% Options: 'dirichlet', 'neumann', 'free', or 'noflow'
params.bottom_bc_type  = "free";
params.bottom_bc_value = 0;         % [m] (ignored for 'free' and 'noflow')

% === Surface Flux (Top) — from Hydrologic Catchment Model ===============
% Catchment_Outputs returns:
%   • surface_flux_time [s]
%   • surface_flux_vals [m/s]
%   • C_top (optional water quality tracer) [mg / L]

% Example of time-varying surface fluxes and bottom fluxes, uncomment if
% you want to manually define the surface and bottom fluxes
% dt_input = 60;
% Nt_input = round(params.Tmax / dt_input);
% t_input  = (0:Nt_input - 1) * dt_input;
% surface_flux_vals   = -1e-4 * sin(linspace(0, pi, Nt_input))';   % Example rainfall or use a constant value
% bottom_flux_vals = 0 * surface_flux_vals;

if params.top_bc_type == "neumann"
    [surface_flux_time, surface_flux_vals, C_top] = Catchment_Outputs(params.dt , params.LID_area);  
    params.surface_flux_time = surface_flux_time;      % Time vector for interpolation [s]
    params.surface_flux_vals = surface_flux_vals;      % Time-varying flux [m/s]
else
    dt_input = 60; Nt_input = round(params.Tmax / dt_input); t_input  = (0:Nt_input - 1) * dt_input;
    params.surface_flux_time = Nt_input; params.surface_flux_vals = 0 * sin(linspace(0, pi, Nt_input)); % Example rainfall or use a constant value
end

% params.surface_flux_vals = -ones(1,length(params.surface_flux_vals)) * 1000 / 1000 / 3600; % 10 mm/h

% === Bottom Flux — (e.g., Recharge, Lateral Leakage, etc.) ==============
% You can assign another hydrologic model or use constant zero
params.bottom_flux_time  = params.surface_flux_time;
params.bottom_flux_vals  = 0 * params.surface_flux_vals;   % No recharge in this example

% === STRUCTURAL DRAINAGE SINKS ==========================================

% === Orifice Flow Parameters =============================================
% General form: Q_orifice = K_orifice * (max(h, 0)) ^ exp_orifice
params.K_orifice     = zeros(1, params.Nz);             % [m^(e)/s] Coefficient per node
params.exp_orifice   = 0.5 * ones(1, params.Nz);         % [-] Exponent (typically 0.5 for orifices)

% Example: Add orifice at node 1 (10 cm diameter)
node_idx = [2];                   % Vector indicating orifice nodes
n_orifices = 0;                 % Vector indicating the number of orifices in a node
Cd       = 0.6;                 % Discharge coefficient [-] (0.5 - 0.6)
D        = [0.1];                % Vector of orifice diameters [m]
Aeff     = pi * D .^2 / 4;      % Effective area [m²]
g        = 9.81;                % Gravity [m/s²]
params.K_orifice(node_idx) = n_orifices .* Cd .* Aeff * sqrt(2 * g);   % Full orifice flow coefficient

% === Spillway Flow at Top Node (for Neumann BC only) =====================
% General form: Q_spillway = c_spillway * (h - h_spill) ^ exp_spillway, when h > h_spill
params.spillway_enabled   = true;              % Enable spillway only for Neumann BC
params.c_spillway         = 0 *1.8 * 1.5;      % Cd * L [m^(e)/s] for weir-type equation
params.h_spill            = 0.05;              % Activation height above soil [m]
params.exp_spillway       = 1.5;               % Exponent for free surface overflow

% === Output Arrays to Track Flow =========================================
Q_orifice_total  = zeros(1, n_save);  % [m³/m²/s] per time step
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s] per time step

% === B.C and Ponding Initial Values ======================================
if params.top_bc_type == "dirichlet"
    ponding_depth = max(params.top_bc_value, 0);
    ponding_prev  = ponding_depth;
    top_val       = params.top_bc_value;
else
    ponding_depth = 0;
    ponding_prev  = 0;
    top_val       = 0;
end

% === 6. SOURCE TERM AND INITIAL CONDITIONS ===============================

% Create source profile [Nz x Nt] – default: zero
params.source_times   = linspace(0, params.Tmax, params.Nt);
params.source_profile = zeros(params.Nz, params.Nt);  % Fill as needed

% Initial condition: hydrostatic profile
p = -1;                       % Uniform suction [m]
h = p * ones(1, params.Nz);   % Default

% Hydrostatic Equilibrium
% p = -0.4;                    % Bottom Pressure
% h = p + [0, cumsum(params.dz(1:end-1))];

if params.bottom_bc_type == "dirichlet"
    h(1) = params.bottom_bc_value;
end
if params.top_bc_type == "dirichlet"
    h(end) = params.top_bc_value;
end

% === 📄 7. Save Simulation Log/Metadata =====================================
logfile = fullfile(base_output_dir, 'Log.txt');
fid = fopen(logfile, 'w');

fprintf(fid, '=====================================================\n');
fprintf(fid, '🔧 Mixed-Form Richards Model — Simulation Log\n');
fprintf(fid, '=====================================================\n\n');

% === Simulation Info
fprintf(fid, '🟢 Simulation Name     : %s\n', sim_name);
fprintf(fid, '📅 Date and Time       : %s\n', datestr(now));
fprintf(fid, '⏱️ Duration            : %.2f hr\n', params.Tmax / 3600);
fprintf(fid, '💾 Save Interval       : %.2f min\n', params.save_interval / 60);
fprintf(fid, '\n');

% === Mesh Info
fprintf(fid, '📐 Mesh Information\n');
fprintf(fid, '    • Number of Nodes      : %d\n', params.Nz);
fprintf(fid, '    • Total Depth          : %.2f m\n', abs(params.z(1)));
fprintf(fid, '    • Top Node Elevation   : %.2f m\n', params.z(end));
fprintf(fid, '    • Min Cell Thickness   : %.4f m\n', min(params.dz));
fprintf(fid, '    • Max Cell Thickness   : %.4f m\n', max(params.dz));
fprintf(fid, '\n');

% === Soil Properties
fprintf(fid, '🌱 Soil Hydraulic Properties (Van Genuchten)\n');
fprintf(fid, '    • α       : %.4f 1/m\n', mean(params.alpha));
fprintf(fid, '    • n       : %.4f [-]\n', mean(params.n));
fprintf(fid, '    • m       : %.4f [-]\n', mean(params.m));
fprintf(fid, '    • θ_s     : %.4f [m³/m³]\n', mean(params.theta_s));
fprintf(fid, '    • θ_r     : %.4f [m³/m³]\n', mean(params.theta_r));
fprintf(fid, '    • K_s     : %.2e [m/s]\n', mean(params.Ks));
fprintf(fid, '    • S_s     : %.2e [1/m]\n', mean(params.S_s));
fprintf(fid, '\n');

% === Boundary Conditions
fprintf(fid, '🔲 Boundary Conditions\n');
fprintf(fid, '    • Top BC Type     : %s\n', params.top_bc_type);
if isfield(params, 'top_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.top_bc_value);
end
fprintf(fid, '    • Bottom BC Type  : %s\n', params.bottom_bc_type);
if isfield(params, 'bottom_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.bottom_bc_value);
end
fprintf(fid, '\n');

% === Drainage Sinks Information ==========================================
fprintf(fid, '💧 Drainage Sinks Information\n');

% Orifice information
if any(params.K_orifice > 0)
    fprintf(fid, '    • Orifices Enabled     : %s\n', "Yes");
else
    fprintf(fid, '    • Orifices Enabled     : %s\n', "No");
end
if any(params.K_orifice > 0)
    fprintf(fid, '    • Orifice Exponent(s)  : %.2f (for each node)\n', mean(params.exp_orifice));
    fprintf(fid, '    • Orifice Coefficients : %.4f [m^(exp_orifice)/s]\n', mean(params.K_orifice(params.K_orifice > 0)));
end

% Spillway information
if params.spillway_enabled
    fprintf(fid, '    • Spillway Enabled     : %s\n', "Yes");
else
    fprintf(fid, '    • Spillway Enabled     : %s\n', "No");
end
if params.spillway_enabled
    fprintf(fid, '    • Spillway Coefficient : %.4f [m^(exp_spillway)/s]\n', params.c_spillway);
    fprintf(fid, '    • Spillway Exponent    : %.2f\n', params.exp_spillway);
    fprintf(fid, '    • Spillway Height      : %.2f [m]\n', params.h_spill);
end

fprintf(fid, '\n');

% === Solver and Time-Stepping
fprintf(fid, '⚙️ Solver Configuration\n');
fprintf(fid, '    • Initial dt            : %.4f s\n', params.dt);
fprintf(fid, '    • Min/Max dt            : [%.4f, %.4f] s\n', params.dt_min, params.dt_max);
fprintf(fid, '    • Max Newton Iterations : %d\n', params.max_iters);
fprintf(fid, '    • Adapt Up/Down         : [%.2f, %.2f]\n', params.adapt_up, params.adapt_down);
fprintf(fid, '    • Convergence Threshold : %.1e\n', params.tol);
fprintf(fid, '\n');

% === Source/Sink Terms
if isfield(params, 'source_profile')
    total_source = sum(params.source_profile(:));
    fprintf(fid, '🛠️  Source Term Included  : Yes\n');
    fprintf(fid, '    • Total Integrated Source [m³/m²]: %.4f\n', total_source * sum(params.dz));
else
    fprintf(fid, '🛠️  Source Term Included  : No\n');
end
fprintf(fid, '\n');

fprintf(fid, '=====================================================\n');
fprintf(fid, '📍 End of Log\n');
fprintf(fid, '=====================================================\n');

fclose(fid);

% === SAVE FILE ===========================================================
save('Examples/Ev6_May13_0110.mat')

%% =========================================================================
% 📘 EXAMPLE 5: Large Bioretention to Receive a small urban catchment
% =========================================================================
clear; clc;

% === 🗂️ Simulation Name and Directory Setup =============================
sim_name = 'Large_Bioretention';

% Define subdirectories using sim_name
base_output_dir   = fullfile('Modeling_Results', sim_name);
figures_dir       = fullfile(base_output_dir, 'Figures');
data_dir         = fullfile(base_output_dir, 'Data');
mesh_dir          = fullfile(figures_dir, 'Mesh');
profiles_dir      = fullfile(figures_dir, 'Profiles');
time_series_dir   = fullfile(figures_dir, 'TimeSeries');
diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
fdc_dir           = fullfile(figures_dir, 'FlowDuration');
wetting_dir       = fullfile(figures_dir, 'WettingFront');

% Create each directory only if it does not already exist
dirs_to_create = {base_output_dir, figures_dir, mesh_dir, profiles_dir, ...
                  time_series_dir, diagnostics_dir, fdc_dir, data_dir, wetting_dir};

for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end

% === DOMAIN SETUP =======================================================
params.Nz = 41;
params.L  = 1.0;                      % Total soil column depth [m]
nonlin_factor = 1.0;                 % Uniform mesh
[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);
params.tol = 1e-8;

% === TIME SETTINGS ======================================================
params.Tmax        = 3600 * 24;       % 48 hours
params.dt          = 15*60;            % [s] = 5 min
params.dt_min      = 60;
params.dt_max      = 15*60;
params.adapt_down  = 0.5;
params.adapt_up    = 1.2;
params.n_up        = 3;
params.n_down      = 10;
params.max_iters   = 100;
params.Nt          = round(params.Tmax / params.dt);

params.save_interval_min = 15*60;        % every 15 min
params.save_interval     = max(params.save_interval_min * 60, params.dt);  % [s]
params.save_steps        = round(0:params.save_interval / params.dt : params.Tmax / params.dt);
n_save                   = length(params.save_steps);
params.plot_interval     = round(0.05 * params.Nt);

% === SOIL PROPERTIES FROM CELIA ET AL. (1990) ============================
params.alpha   = 0.0355 * 100 * ones(1, params.Nz);    % [1/m]
params.n       = 2.0 * ones(1, params.Nz);             % [-]
params.m       = 1 - 1 ./ params.n;
params.theta_r = 0.102 * ones(1, params.Nz);           % [m³/m³]
params.theta_s = 0.368 * ones(1, params.Nz);           % [m³/m³]
params.S_s     = 0 * ones(1, params.Nz);               % [1/m]
params.Ks      = 0.00922 / 100 * ones(1, params.Nz);   % [m/s]

% === BOUNDARY CONDITIONS ================================================
params.top_bc_type     = "neumann";     % Rainfall / Evaporation
params.top_bc_value    = nan;           % [m]
params.bottom_bc_type  = "free";        % Free drainage
params.bottom_bc_value = nan;           % [m]


% Example: load or define your flux input here
% Assume q_flux_top is [Nt_input x 1] in m/s, and t_flux is in seconds
[t_flux, q_flux_top, ~] = Catchment_Outputs(params.dt / 60);  % ← You define this
params.top_flux_time = t_flux;     % Time vector [s]
params.top_flux_vals = q_flux_top; % Flux at top boundary [m/s]


% === SOURCES/SINKS (NONE) ================================================
params.source_times   = linspace(0, params.Tmax, params.Nt);
params.source_profile = zeros(params.Nz, params.Nt);

dt_input        = 60;
Nt_input        = round(params.Tmax / dt_input);
t_input         = (0:Nt_input - 1) * dt_input;
rain_flux_vals  = zeros(Nt_input, 1);
bottom_flux_vals = zeros(Nt_input, 1);

% Example: load or define your flux input here
% Assume q_flux_top is [Nt_input x 1] in m/s, and t_flux is in seconds
[t_flux, q_flux_top, C_top] = Catchment_Outputs(params.dt / 60);  % ← You define this
params.surface_flux_time = t_flux;     % Time vector [s]
params.surface_flux_vals = q_flux_top; % Flux at top boundary [m/s]
params.bottom_flux_time  = t_input;
params.bottom_flux_vals  = bottom_flux_vals;

% === INITIAL CONDITION: DRY PROFILE ======================================
h = -10 * ones(1 , params.Nz);  % dry soil (initial condition)
if params.top_bc_type == "dirichlet"
    h(end) = params.top_bc_value;
end
if params.bottom_bc_type == "dirichlet"
    h(1) = params.bottom_bc_value;
end

% === OUTPUT PREALLOCATION ====================================================
head_out         = nan(params.Nz, n_save);
theta_out        = nan(params.Nz, n_save);
flux_out         = nan(params.Nz + 1, n_save);
ponding_series   = zeros(1, n_save);
seepage_flux     = zeros(1, n_save);
top_flux         = zeros(1, n_save);
save_count       = 0;
mb_error_cumulative = 0;
ponding_depth = 0;

% Plotting time vector
n_plots = 20;
plot_times = linspace(0, params.Tmax, n_plots);
plot_index = 1;

% Pressure Limiter for Evaporation (Feddes)
params.h_lim_upper = -10; % [m]
params.h_lim_down = -0.1; % [m]

% === Set top_val and bottom_val based on BC type ============================
if params.top_bc_type == "dirichlet"
    top_val = params.top_bc_value;
else
    top_val = nan;
end

if params.bottom_bc_type == "dirichlet"
    bottom_val = params.bottom_bc_value;
else
    bottom_val = nan;
end

% === DRAINAGE CONFIGURATION =====================================
params.K_orifice    = zeros(params.Nz, 1);   % [m^(e_o)/s] Orifice coefficients (default: 0 [not activated])
params.exp_orifice    = 0.5 * ones(params.Nz, 1);    % [–] Exponent for each orifice (usually 0.5)

% Example: Place orifice at node node_idx with K = 0.0216 m^(0.5)/s, assuming a 10
% cm orifice
node_idx = [];
params.K_orifice(node_idx) = 0.0216 * pi() * (0.10) ^ 2 / 4 * sqrt(2 * 9.81); % Cd * Aef * sqrt(2 * g)
params.exp_orifice(node_idx) = 0.5;

% Spillway at top (only used if top BC is Neumann)
params.spillway_enabled   = true;
params.c_spillway         = 1.80 * 1.5;   % [m^(e)/s] -> Cd * Lef
params.h_spill            = 0.05;         % [m] (height above soil where spillway activates)
params.exp_spillway       = 1.5;          % [–]

Q_orifice_total = zeros(1, n_save);   % [m³/m²/s]
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s]


Q_orifice_total = zeros(1, n_save);   % [m³/m²/s]
Q_spillway_total = zeros(1, n_save);  % [m³/m²/s]


% === 📄 Save Simulation Log/Metadata =====================================
logfile = fullfile(base_output_dir, 'Log.txt');
fid = fopen(logfile, 'w');

fprintf(fid, '=====================================================\n');
fprintf(fid, '🔧 Mixed-Form Richards Model — Simulation Log\n');
fprintf(fid, '=====================================================\n\n');

% === Simulation Info
fprintf(fid, '🟢 Simulation Name     : %s\n', sim_name);
fprintf(fid, '📅 Date and Time       : %s\n', datestr(now));
fprintf(fid, '⏱️ Duration            : %.2f hr\n', params.Tmax / 3600);
fprintf(fid, '💾 Save Interval       : %.2f min\n', params.save_interval / 60);
fprintf(fid, '\n');

% === Mesh Info
fprintf(fid, '📐 Mesh Information\n');
fprintf(fid, '    • Number of Nodes      : %d\n', params.Nz);
fprintf(fid, '    • Total Depth          : %.2f m\n', abs(params.z(1)));
fprintf(fid, '    • Top Node Elevation   : %.2f m\n', params.z(end));
fprintf(fid, '    • Min Cell Thickness   : %.4f m\n', min(params.dz));
fprintf(fid, '    • Max Cell Thickness   : %.4f m\n', max(params.dz));
fprintf(fid, '\n');

% === Soil Properties
fprintf(fid, '🌱 Soil Hydraulic Properties (Van Genuchten)\n');
fprintf(fid, '    • α       : %.4f 1/m\n', mean(params.alpha));
fprintf(fid, '    • n       : %.4f [-]\n', mean(params.n));
fprintf(fid, '    • m       : %.4f [-]\n', mean(params.m));
fprintf(fid, '    • θ_s     : %.4f [m³/m³]\n', mean(params.theta_s));
fprintf(fid, '    • θ_r     : %.4f [m³/m³]\n', mean(params.theta_r));
fprintf(fid, '    • K_s     : %.2e [m/s]\n', mean(params.Ks));
fprintf(fid, '    • S_s     : %.2e [1/m]\n', mean(params.S_s));
fprintf(fid, '\n');

% === Boundary Conditions
fprintf(fid, '🔲 Boundary Conditions\n');
fprintf(fid, '    • Top BC Type     : %s\n', params.top_bc_type);
if isfield(params, 'top_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.top_bc_value);
end
fprintf(fid, '    • Bottom BC Type  : %s\n', params.bottom_bc_type);
if isfield(params, 'bottom_bc_value')
    fprintf(fid, '      > Value         : %.4f\n', params.bottom_bc_value);
end
fprintf(fid, '\n');

% === Solver and Time-Stepping
fprintf(fid, '⚙️ Solver Configuration\n');
fprintf(fid, '    • Initial dt            : %.4f s\n', params.dt);
fprintf(fid, '    • Min/Max dt            : [%.4f, %.4f] s\n', params.dt_min, params.dt_max);
fprintf(fid, '    • Max Newton Iterations : %d\n', params.max_iters);
fprintf(fid, '    • Adapt Up/Down         : [%.2f, %.2f]\n', params.adapt_up, params.adapt_down);
fprintf(fid, '    • Convergence Threshold : %.1e\n', params.tol);
fprintf(fid, '\n');

% === Source/Sink Terms
if isfield(params, 'source_profile')
    total_source = sum(params.source_profile(:));
    fprintf(fid, '🛠️  Source Term Included  : Yes\n');
    fprintf(fid, '    • Total Integrated Source [m³/m²]: %.4f\n', total_source * sum(params.dz));
else
    fprintf(fid, '🛠️  Source Term Included  : No\n');
end
fprintf(fid, '\n');

% === Structural Drainage (Orifices and Spillways)
fprintf(fid, '🚰 Structural Drainage Sinks\n');
if isfield(params, 'K_orifice') && any(params.K_orifice > 0)
    idxs = find(params.K_orifice > 0);
    for i = 1:length(idxs)
        fprintf(fid, '    • Orifice at Node %d: K = %.2e, e = %.2f\n', ...
            idxs(i), params.K_orifice(idxs(i)), params.exp_orifice(idxs(i)));
    end
else
    fprintf(fid, '    • No orifices enabled.\n');
end

if isfield(params, 'spillway_enabled') && params.spillway_enabled
    fprintf(fid, '    • Spillway at Surface: c = %.2e, h_s = %.3f m, e = %.2f\n', ...
        params.c_spillway, params.h_spill, params.exp_spillway);
else
    fprintf(fid, '    • No spillway enabled.\n');
end
fprintf(fid, '\n');

fprintf(fid, '=====================================================\n');
fprintf(fid, '📍 End of Log\n');
fprintf(fid, '=====================================================\n');

fclose(fid);


% === SAVE FILE ===========================================================
save('Examples/Large_Bioretention.mat')


