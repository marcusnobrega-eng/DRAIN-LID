%% === 🗂️ User-defined Simulation Name and Output Paths ===================
sim_name = 'Large_Bioretention';  % 📝 CHANGE THIS for each run

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

params.Nz = 27;               % Number of vertical nodes [-]
params.L  = 1.0;              % Total soil depth [m]
nonlin_factor = 1;            % Grid refinement factor (1 = uniform)
params.LID_area = 1000;       % 1D column area [m2]

% Generate refined mesh (Hydrus-style, refined near surface)
[params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);

% === 2. TIME DISCRETIZATION ==============================================

params.Tmax = 24*3600;        % Total simulation time [s]
params.dt   = 5*60;           % Initial time step [s]
params.dt_min = 0.01;           % Minimum dt [s]
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
params.h_lim_upper = -0.1; % [m]
params.h_lim_down = -4; % [m]

% === 3. MULTILAYER SOIL HYDRAULIC PROPERTIES =============================

% Define soil layer thicknesses [m], from bottom to top
media_thicknesses = [0.3, 0.2, 0.1, 0.4];
media_interfaces = [-params.L + cumsum(media_thicknesses)];
media_interfaces = [-params.L, media_interfaces];  % Include bottom
n_layers = length(media_thicknesses);

% Van Genuchten + Ks parameters (per layer)
media_props = struct( ...
    'alpha',   [2.0, 2.0, 2.0, 2.0], ...
    'n',       [1.5, 1.5, 1.5, 1.5], ...
    'theta_r', [0.04, 0.04, 0.04, 0.04], ...
    'theta_s', [0.42, 0.42, 0.42, 0.42], ...
    'S_s',     [1e-4, 1e-4, 1e-4, 1e-4], ...
    'Ks',      [5e-4, 5e-4, 5e-4, 5e-4] ...
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

media_props.labels = {'S1', ...
                      'S2', ...
                      'S3', ...
                      'S4'};

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
params.top_bc_value = 0.05;         % Pressure head at top [m] (set NaN for Neumann)

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
n_orifices = 1;                 % Vector indicating the number of orifices in a node
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
p = 0;                        % Uniform suction [m]
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

