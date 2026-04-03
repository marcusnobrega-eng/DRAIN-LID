%% DRAIN-LID database generator
% Builds a hydrologic performance database by running DRAIN-LID
% over a discretized design and hydrologic space.
% Time resolution: 5 minutes.

%% DRAIN-LID database generator

clear all; clc;

% --- NEW: initialize base parameters once ---
base_params = init_base_params();

% 🔹 Simulation Start Time (Used for time-series plotting)
% Format: datetime(YYYY, MM, DD, HH, MM, SS)
start_datetime = datetime(2023, 12, 21, 0, 10, 0);

% 🔹 Empirical Evaporation / Evapotranspiration Coefficients
% These affect the surface water balance.
% Reference values (see manuscript draft for details):
%   - Delta: 0.3 for Permeable Pavement (PP), 1 for others
%   - Gamma: 2 for Permeable Pavement (PP), 0 for others

base_params.Delta = 1;   % General default
base_params.Gamma = 0;   % General default

%% 1. Discrete sets and constants

% Hydrologic inputs (catchment-scale effective precipitation, mm)
% Peff_set_mm = [5, 10, 25, 50, 75, 100, 200, 500];   % [mm] over upstream catchment
%% 1. Discrete sets and constants

% Hydrologic inputs (catchment-scale effective precipitation, mm)
Peff_set_mm = [0.5]*25.4;   % [mm] over upstream catchment

% Upstream catchment areas [m^2]
Aup_set_m2  = [1];   % A_up

% Design variable: contributing-area ratio (A_up / A_TC)
% You can adjust these as you like
area_ratio_set = [0.01:0.01:0.2];   % (ATC / Aup)

% Time to peak of inflow hydrograph (minutes)
% tp_set_min  = [15, 30, 60];        % [min]
% tp_set_min  = [5, 10, 15, 30];        % [min]
tp_set_min = [60]; 
% tp_set_min = 5;


Ld_set_m    = [0.4:0.1:1.2]; % [m] media depth
% Ld_set_m    = [0.3 0.4 0.5 0.6]; % [m] media depth


% Soil hydraulic parameter sets (van Genuchten–Mualem)
% TODO: replace with your actual calibrated table (from HYDRUS / Carsel & Parrish)
Soils = define_soils();   % returns struct array with fields:
                          % name, theta_r, theta_s, alpha_1perm, n, Ks_mm_per_h
Soils = Soils(3);
nSoil      = numel(Soils);
nAup       = numel(Aup_set_m2);
nAreaRatio = numel(area_ratio_set);
nLd        = numel(Ld_set_m);
nPeff      = numel(Peff_set_mm);
nTp        = numel(tp_set_min);

% Ponding constraint
h_max_allow = 1;       % [m] = 100 cm

% Nash hydrograph parameters
m_shape     = 4.69;       % shape factor
b_factor    = 0.75;       % peaking factor (as in your derivation)

% Simulation time settings
T_hours      = 48;   % total simulation duration [h]

% Internal resolution for generating Nash hydrograph
dt_gen_min   = 1;    % [min] (fine grid for Qin and forcing)

% Output / database resolution
dt_save_min  = 5;    % [min] (saved Qin_all, Qout_all, metrics)

% Fine grid for generating inflow and feeding solver
t_minutes_fine = (0:dt_gen_min:(T_hours*60))';  % [min]
t_hours_fine   = t_minutes_fine / 60;           % [h]

% Coarse grid for saved database
t_minutes   = (0:dt_save_min:(T_hours*60))';    % [min]
t_hours     = t_minutes / 60;                   % [h]
nSteps      = numel(t_hours);
%% 2. Preallocate storage

% Total number of combinations (before feasibility filtering)
N_total = nSoil * nAup * nAreaRatio * nLd * nPeff * nTp;

Metrics = struct( ...
    'soilName',      cell(N_total,1), ...
    'Aup_m2',        nan(N_total,1), ...   % upstream area
    'ATC_m2',        nan(N_total,1), ...   % LID area (computed from ratio)
    'area_ratio',    nan(N_total,1), ...   % A_up / A_TC
    'Ld_m',          nan(N_total,1), ...
    'Peff_catch_mm', nan(N_total,1), ...
    'Peff_eq_mm',    nan(N_total,1), ...
    'tp_min',        nan(N_total,1), ...
    'Qp_in_m3s',     nan(N_total,1), ...
    'Qp_out_m3s',    nan(N_total,1), ...
    'eta_p_pct',     nan(N_total,1), ...
    'Delta_tp_min',  nan(N_total,1), ...
    'DetTime_min',   nan(N_total,1), ...
    'DeltaV_m3',     nan(N_total,1), ...
    'max_ponding',   nan(N_total,1));

% Total number of combinations (before feasibility filtering)
N_total = nSoil * nAup * nAreaRatio * nLd * nPeff * nTp;

% Progress counter (attempted scenarios)
combo_idx = 0;

% Time series: store Q_in and Q_out at 5-min resolution
Qin_all  = nan(nSteps, N_total);   % [m^3/s]
Qout_all = nan(nSteps, N_total);   % [m^3/s]
Se_top_all    = nan(nSteps, N_total);   % [-]
Se_mid_all    = nan(nSteps, N_total);   % [-]
Se_bottom_all = nan(nSteps, N_total);   % [-]

k = 0;   % feasible scenario counter

%% 3. Main nested loops

for iSoil = 1:nSoil
    soil = Soils(iSoil);

    for iAup = 1:nAup
        A_up = Aup_set_m2(iAup);    % [m^2] upstream catchment area

        for iAR = 1:nAreaRatio
            area_ratio = area_ratio_set(iAR);   % A_up / A_TC (input)
            ATC        = A_up * area_ratio;     % [m^2] LID area (computed)

            % (Optional) skip if ATC gets too small/large, e.g.
            % if ATC < 1 || ATC > 1000
            %     continue;
            % end

            for iLd = 1:nLd
                Ld = Ld_set_m(iLd);  % [m]

                for iPeff = 1:nPeff
                    Peff_catch_mm = Peff_set_mm(iPeff); % [mm] over A_up

                    % Equivalent effective depth over bioretention footprint [mm]
                    Peff_eq_mm = Peff_catch_mm / area_ratio;  % [mm over A_TC]
                    Peff_eq_m  = Peff_eq_mm / 1000;           % [m]

                    % Runoff volume over bioretention surface [m^3]
                    V_run = Peff_eq_m * ATC;

                    for iTp = 1:nTp
                        tp_min = tp_set_min(iTp);     % [min]
                        tp_hr  = tp_min / 60;         % [h]

                        % ==== 🔄 Progress tracking =======================
                        combo_idx = combo_idx + 1;
                        if mod(combo_idx, 50) == 0 || combo_idx == 1 || combo_idx == N_total
                            progress = 100 * combo_idx / N_total;
                            fprintf('Progress: %6.2f %%  (%d / %d scenarios attempted)\n', ...
                                    progress, combo_idx, N_total);
                        end
                        % =================================================

                        % Generate inflow hydrograph (still using t_hours)
                        Qin = generate_nash_hydrograph( ...
                                  t_hours, V_run, m_shape, tp_hr);

                        % 🔹 Define inflow peak (used later in Metrics)
                        Qp_in = max(Qin);   % [m^3/s]

                        % --- Call DRAIN-LID Richards solver ---
                        [Qout, h_surf, t_save_hours, Se_top, Se_mid, Se_bottom] = run_drain_lid( ...
                            t_hours, Qin, ATC, Ld, soil, base_params);

                        % If solver failed or returned invalid outputs, skip
                        if isempty(Qout) || isempty(h_surf) || all(isnan(Qout))
                            continue;
                        end

                        % Feasibility check: maximum ponding
                        h_max_sim = max(h_surf);

                        if h_max_sim > h_max_allow
                            % undersized for this load; skip storing
                            continue;
                        end

                        % If feasible, increment counter and store results
                        k = k + 1;

                        % Store time series
                        Qin_all(:, k)  = Qin;
                        Qout_all(:, k) = Qout;

                        %%% store Se time series (top, mid, bottom)
                        Se_top_all(:,    k) = Se_top(:);
                        Se_mid_all(:,    k) = Se_mid(:);
                        Se_bottom_all(:, k) = Se_bottom(:);

                        % Compute metrics
                        [eta_p, Delta_tp_min, DetTime_min, Qp_out, DeltaV_m3] = ...
                            compute_metrics(t_minutes, Qin, Qout);

                        % Fill metrics struct
                        Metrics(k).soilName      = soil.name;
                        Metrics(k).Aup_m2        = A_up;
                        Metrics(k).ATC_m2        = ATC;
                        Metrics(k).area_ratio    = area_ratio;
                        Metrics(k).Ld_m          = Ld;
                        Metrics(k).Peff_catch_mm = Peff_catch_mm;
                        Metrics(k).Peff_eq_mm    = Peff_eq_mm;
                        Metrics(k).tp_min        = tp_min;
                        Metrics(k).Qp_in_m3s     = Qp_in;
                        Metrics(k).Qp_out_m3s    = Qp_out;
                        Metrics(k).eta_p_pct     = eta_p;
                        Metrics(k).Delta_tp_min  = Delta_tp_min;
                        Metrics(k).DetTime_min   = DetTime_min;
                        Metrics(k).DeltaV_m3     = DeltaV_m3;
                        Metrics(k).max_ponding   = h_max_sim;

                    end % iTp
                end % iPeff
            end % iLd
        end % iATC
    end % iAup
end % iSoil

%% 4. Trim unused preallocated entries

if k < N_total
    Metrics   = Metrics(1:k);
    Qin_all   = Qin_all(:, 1:k);
    Qout_all  = Qout_all(:, 1:k);
    Se_top_all     = Se_top_all(:, 1:k);      
    Se_mid_all     = Se_mid_all(:, 1:k);      
    Se_bottom_all  = Se_bottom_all(:, 1:k); 
end

% Removing cases with no outflow
idx = max(Qout_all) == 0;      % 1 x k logical (per scenario/column)

Metrics(idx)        = [];
Qin_all(:,  idx)    = [];
Qout_all(:, idx)    = [];
Se_top_all(:,    idx) = [];    %%%% NEW
Se_mid_all(:,    idx) = [];    %%%% NEW
Se_bottom_all(:, idx) = [];    %%%% NEW

%% 5. Save database

save('DRAIN_LID_database_Example2_SandyLoam_Loam.mat', ...
     'Metrics', 'Qin_all', 'Qout_all', ...
     'Se_top_all', 'Se_mid_all', 'Se_bottom_all', ...   %%%% NEW
     't_minutes', 'Soils', ...
     'Peff_set_mm', 'tp_set_min', 'Aup_set_m2', 'area_ratio_set', 'Ld_set_m', ...
     'h_max_allow', 'm_shape', 'b_factor');

disp(['Finished. Number of feasible scenarios: ', num2str(k)]);

function Soils = define_soils()
%DEFINE_SOILS  Returns a struct array with van Genuchten–Mualem parameters
%               for USDA soil texture classes.
%
% Units:
%   theta_r        [-]    residual volumetric water content
%   theta_s        [-]    saturated volumetric water content
%   alpha_1perm    [1/m]  van Genuchten alpha
%   n              [-]    van Genuchten n
%   Ks_mm_per_h    [mm/h] saturated hydraulic conductivity
%   Ss             [1/m]  specific storage (from table; constant here)
%
% Values adapted from Carsel & Parrish (1988), see your Table 5.
% Ks in the paper is in [m/s]; here it is converted to [mm/h] via:
%     Ks_mm_per_h = Ks_m_per_s * 3.6e6

Ss_val = 1e-5;  % [1/m] for all textures in the table

i = 1;
Soils(i).name        = 'Sand';
Soils(i).theta_r     = 0.045;
Soils(i).theta_s     = 0.43;
Soils(i).alpha_1perm = 14.5;        % [1/m]
Soils(i).n           = 2.68;
Soils(i).Ks_mm_per_h = 297.0;       % 8.25e-5 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Loamy Sand';
Soils(i).theta_r     = 0.057;
Soils(i).theta_s     = 0.41;
Soils(i).alpha_1perm = 12.4;
Soils(i).n           = 2.28;
Soils(i).Ks_mm_per_h = 145.8;       % 4.05e-5 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Sandy Loam';
Soils(i).theta_r     = 0.065;
Soils(i).theta_s     = 0.41;
Soils(i).alpha_1perm = 7.5;
Soils(i).n           = 1.89;
Soils(i).Ks_mm_per_h = 44.3;        % 1.23e-5 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Loam';
Soils(i).theta_r     = 0.078;
Soils(i).theta_s     = 0.43;
Soils(i).alpha_1perm = 3.6;
Soils(i).n           = 1.56;
Soils(i).Ks_mm_per_h = 10.4;        % 2.89e-6 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Silt';
Soils(i).theta_r     = 0.034;
Soils(i).theta_s     = 0.46;
Soils(i).alpha_1perm = 1.6;
Soils(i).n           = 1.37;
Soils(i).Ks_mm_per_h = 2.5;         % 6.94e-7 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Silty Loam';
Soils(i).theta_r     = 0.067;
Soils(i).theta_s     = 0.45;
Soils(i).alpha_1perm = 2.0;
Soils(i).n           = 1.41;
Soils(i).Ks_mm_per_h = 4.5;         % 1.25e-6 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Sandy Clay Loam';
Soils(i).theta_r     = 0.100;
Soils(i).theta_s     = 0.39;
Soils(i).alpha_1perm = 5.9;
Soils(i).n           = 1.48;
Soils(i).Ks_mm_per_h = 13.1;        % 3.64e-6 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Clay Loam';
Soils(i).theta_r     = 0.095;
Soils(i).theta_s     = 0.41;
Soils(i).alpha_1perm = 1.9;
Soils(i).n           = 1.31;
Soils(i).Ks_mm_per_h = 2.6;         % 7.22e-7 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Silty Clay Loam';
Soils(i).theta_r     = 0.089;
Soils(i).theta_s     = 0.43;
Soils(i).alpha_1perm = 1.0;
Soils(i).n           = 1.23;
Soils(i).Ks_mm_per_h = 0.70;        % 1.94e-7 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Sandy Clay';
Soils(i).theta_r     = 0.100;
Soils(i).theta_s     = 0.38;
Soils(i).alpha_1perm = 2.7;
Soils(i).n           = 1.23;
Soils(i).Ks_mm_per_h = 1.20;        % 3.33e-7 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Silty Clay';
Soils(i).theta_r     = 0.070;
Soils(i).theta_s     = 0.36;
Soils(i).alpha_1perm = 0.5;
Soils(i).n           = 1.09;
Soils(i).Ks_mm_per_h = 0.20;        % 5.56e-8 m/s
Soils(i).Ss          = Ss_val;

i = i + 1;
Soils(i).name        = 'Clay';
Soils(i).theta_r     = 0.068;
Soils(i).theta_s     = 0.38;
Soils(i).alpha_1perm = 0.8;
Soils(i).n           = 1.09;
Soils(i).Ks_mm_per_h = 2.00;        % 5.56e-7 m/s
Soils(i).Ss          = Ss_val;

end

function base_params = init_base_params()
%INIT_BASE_PARAMS  Base configuration for the mixed-form Richards solver
% used in DRAIN-LID. These settings are scenario-independent and can be
% reused across all simulations in the database.

%% === 🔗 Add Required Paths ==============================================
% Adjust this root path if you move the repository
repo_root = 'C:\Users\marcu\Documents\GitHub\DRAIN-LID';

addpath(fullfile(repo_root, 'Physical_Functions_and_Utilities'));
addpath(fullfile(repo_root, 'Numerical_Solver'));
addpath(fullfile(repo_root, 'Main_Functions'));
addpath(fullfile(repo_root, 'Model_Configurations'));

%% 1. DOMAIN & MESH (generic)
base_params.Nz = 20;        % Number of vertical nodes [-]; uniform grid
base_params.L  = 1.0;       % Placeholder total depth [m]; overwritten per run

% 2. TIME DISCRETIZATION
base_params.Tmax    = 36*3600;  % Total simulation time [s]
base_params.dt      = 0.1;      % Initial time step [s]
base_params.dt_min  = 0.1;    % Minimum dt [s]
base_params.dt_max  = 5*60;     % Maximum dt [s]

% Adaptive timestep control
base_params.adapt_down = 0.5;   % Shrink factor
base_params.adapt_up   = 2.0;   % Growth factor
base_params.n_up       = 5;     % Threshold for fast convergence
base_params.n_down     = 10;    % Threshold for slow convergence

base_params.Nt         = round(base_params.Tmax / base_params.dt);
base_params.max_iters  = 20;    % Newton max iterations
base_params.tol        = 1e-6;  % Newton tolerance [m]

% Output saving interval (5-min resolution)
base_params.save_interval_min = 5;            % [min]
base_params.save_interval     = 5*60;         % [s]

% 3. FEDDES-TYPE EVAPORATION PRESSURE LIMITER (kept generic)
base_params.h_lim_upper = -0.0;  % [m]
base_params.h_lim_down  = -4.0;  % [m]

% 4. BOUNDARY CONDITIONS (types; values will be set per scenario)
% Top boundary: Neumann flux from inflow hydrograph
base_params.top_bc_type     = "neumann";   % 'neumann' enforced
base_params.top_bc_value    = NaN;         % not used for Neumann

% Bottom boundary: free drainage by default
base_params.bottom_bc_type  = "free";      % 'free' (gravity driven)
base_params.bottom_bc_value = 0.0;         % not used for 'free'

% 5. STRUCTURAL DRAINAGE SINKS (default: no orifices; spillway at 15 cm)
base_params.K_orifice   = [];  % will be sized per run (no orifices)
base_params.exp_orifice = [];  % (unused for now)

base_params.spillway_enabled = true;  % overflow beyond 0.15 m ponding
base_params.c_spillway       = 0.0;   % default; set per model if needed
base_params.h_spill          = 0.15;  % [m] max ponding depth
base_params.exp_spillway     = 1.5;   % typical weir exponent

% 6. SOURCE TERM (default: none; can be extended)
base_params.source_times   = [];
base_params.source_profile = [];

% 7. INITIAL CONDITIONS (to be set per soil in run_drain_lid)
% Field capacity suction will be assigned later based on soil type
base_params.initial_suction = NaN;   % placeholder; not used directly

% 8. LID AREA (m²) – will be overwritten each run
base_params.LID_area = 1.0;

end

function [Qout, h_surf, t_save_hours, Se_top, Se_mid, Se_bottom] = run_drain_lid(t_hours, Qin, ATC, Ld, soil, base_params)
%RUN_DRAIN_LID  Configure and run one DRAIN-LID simulation for the database.
%
% Inputs:
%   t_hours     : time vector [h] at which Qin is given (5-min resolution)
%   Qin         : inflow hydrograph [m^3/s] at bioretention inlet
%   ATC         : bioretention surface area [m^2]
%   Ld          : media depth [m]
%   soil        : struct with fields:
%                 theta_r, theta_s, alpha_1perm [1/m], n, Ks_mm_per_h
%   base_params : base configuration from init_base_params()
%
% Outputs:
%   Qout        : outflow discharge [m^3/s] (bottom drainage, + spill if added)
%   h_surf      : surface ponding depth [m]
%   t_save_hours: time vector of outputs [h]

% -------------------------------------------------------------------------
% 1. Clone base parameters and set geometry
% -------------------------------------------------------------------------
params = base_params;

params.LID_area = ATC;    % [m^2]
params.L        = Ld;     % modeled depth [m]
params.Nz       = base_params.Nz;

% Uniform vertical grid from -Ld to 0
params.z  = linspace(-Ld, 0, params.Nz);
params.dz = diff(params.z);
params.dz = [params.dz, params.dz(end)];

% Non-linear mesh from -Ld to 0
nonlin_factor = 1.5;            % Grid refinement factor (1 = uniform)

% Generate refined mesh (Hydrus-style, refined near surface)
% 🧭 Uniformly spaced pseudo-depth coordinate from 0 (top) to 1 (bottom)
s = linspace(0, 1, params.Nz)';

% 🎯 Apply nonlinear stretching (exponentially refine near surface)
s_refined = 1 - (1 - s).^nonlin_factor;

% 📏 Convert to actual depth values from 0 (top) to -L (bottom)
params.z = -Ld + Ld * s_refined;  % z(1) = -L, z(end) = 0

% 📐 Compute cell thicknesses (dz_i = z(i+1) - z(i))
params.dz = diff(params.z)';
params.dz = [params.dz, params.dz(end)];  % Last cell thickness = same as previous

% -------------------------------------------------------------------------
% 2. Soil hydraulic properties (single layer, uniform)
% -------------------------------------------------------------------------
params.theta_r = soil.theta_r     * ones(1, params.Nz);
params.theta_s = soil.theta_s     * ones(1, params.Nz);
params.alpha   = soil.alpha_1perm * ones(1, params.Nz);   % [1/m]
params.n       = soil.n           * ones(1, params.Nz);
params.m       = 1 - 1 ./ params.n;
params.S_s     = 1e-5 * ones(1, params.Nz);               % [1/m]

Ks_m_per_s     = soil.Ks_mm_per_h / 1000 / 3600;          % mm/h → m/s
params.Ks      = Ks_m_per_s * ones(1, params.Nz);

% -------------------------------------------------------------------------
% 3. Time discretization consistent with t_hours
% -------------------------------------------------------------------------
t_seconds   = t_hours(:) * 3600;
params.Tmax = t_seconds(end);

params.save_interval_min = 5;
params.save_interval     = params.save_interval_min * 60;  % [s]

% Initial time step can stay as in base_params (e.g., 1 s)
% params.dt already set in base_params

% -------------------------------------------------------------------------
% 4. Boundary conditions: Neumann top, free bottom
% -------------------------------------------------------------------------
params.top_bc_type    = "neumann";
params.bottom_bc_type = "free";

% Surface flux [m/s] = -Qin / area (negative downward)
params.surface_flux_time = t_seconds;
params.surface_flux_vals = -Qin(:) / ATC;  % [m^3/s] / [m^2] = [m/s]

% Bottom flux: zero (free drainage)
params.bottom_flux_time = params.surface_flux_time;
params.bottom_flux_vals = zeros(size(params.surface_flux_vals));

% -------------------------------------------------------------------------
% 5. Structural drainage sinks — no underdrain, spillway at 0.15 m
% -------------------------------------------------------------------------
params.K_orifice   = zeros(1, params.Nz);
params.exp_orifice = 0.5 * ones(1, params.Nz);  % unused if K_orifice = 0

params.spillway_enabled = true;
params.h_spill          = 1;  % [m] maximum ponding depth
params.exp_spillway     = 1.5;

% Spillway coefficient (can be tuned)
params.c_spillway       = 1.0;   % placeholder

% -------------------------------------------------------------------------
% 6. Initial conditions: 50%-saturation-based suction
% -------------------------------------------------------------------------
% We set the initial condition to a uniform matric head h_init such that
% the effective saturation Se = (theta - theta_r)/(theta_s - theta_r)
% is 0.5 everywhere.
%
% Van Genuchten relations:
%   Se(h) = [1 + (alpha * |h|)^n]^(-m),   m = 1 - 1/n
%
% Inverting for h given Se_target:
%   [1 + (alpha * |h|)^n] = Se^(-1/m)
%   (alpha * |h|)^n       = Se^(-1/m) - 1
%   |h|                   = ( Se^(-1/m) - 1 )^(1/n) / alpha
%   h_init                = -|h|  (matric head is negative in unsaturated zone)

% Target effective saturation
Se_target = 0.5;   % 50% of saturation (degree of saturation)

% Extract soil parameters (scalar for this single-layer case)
alpha   = soil.alpha_1perm;  % [1/m]
n       = soil.n;
m       = 1 - 1/n;

% Safety clamp for Se to avoid singularities at exactly 0 or 1
Se = max(min(Se_target, 1 - 1e-6), 1e-6);

% Compute |h| from inverted van Genuchten curve
abs_h = ((Se.^(-1/m) - 1).^(1/n)) ./ alpha;   % [m]

% Matric head (negative for suction)
h_init = -abs_h;   % [m]

% Store as uniform initial suction (solver will build the full profile)
params.initial_suction = h_init;

% -------------------------------------------------------------------------
% 7. Distributed source term: here set to zero everywhere
% -------------------------------------------------------------------------
% We define a zero source term over the full simulation period so that
% interp1() in the solver always has valid input.
params.source_times   = [0, params.Tmax];                       % [s]
params.source_profile = zeros(params.Nz, numel(params.source_times));  % [Nz x 2]

% -------------------------------------------------------------------------
% 8. Run the mixed-form Richards solver
% -------------------------------------------------------------------------
[head_out, theta_out, flux_out, ponding_series, ...
          outlet_flux, time_series, max_ponding_depth, ...
          Se_top, Se_mid, Se_bottom, success] = Main_Solver_Function(params);   %%%% NEW

% If solver failed, return empty so the outer loop can skip this design
if ~success || isempty(time_series) || isempty(outlet_flux)
    Qout   = [];
    h_surf = [];
    t_save_hours = [];
    Se_top       = [];   
    Se_mid       = [];   
    Se_bottom    = [];   
    return;
end
% -------------------------------------------------------------------------
% 9. Map solver outputs to Qout(t) and h_surf(t)
% -------------------------------------------------------------------------
% time_series returned by Main_Solver is [s]; convert to hours
t_save_hours = time_series(:) / 3600;

% Ponding depth at surface
h_surf = ponding_series(:);   % [m]

% Outflow definition:
%   outlet_flux is bottom vertical flux [m/s] (as set inside Main_Solver)
%   Multiply by ATC to get discharge [m^3/s].
Q_bottom = outlet_flux(:) * ATC;

% If you later store spillway/orifice discharges separately, add them here.
Qout = Q_bottom;   % for now: only bottom drainage

Qout = -Qout;   % changing signal for consider outflow correctly

end

function [Qin, Qp_exact, b_exact] = generate_nash_hydrograph(t_hours, V_run, m_shape, tp_hr)
%GENERATE_NASH_HYDROGRAPH Generates a Nash-type inflow hydrograph
%   directly on the prescribed time grid t_hours (in hours), using the
%   exact Gamma/Nash formulation (no discrete renormalization).
%
%   [Qin, Qp_exact, b_exact] = generate_nash_hydrograph(t_hours, V_run, m_shape, tp_hr)
%
%   INPUTS
%       t_hours : time vector [h], e.g. 0:5/60:36
%       V_run   : total runoff volume reaching the bioretention [m^3]
%       m_shape : Nash shape parameter m [-]
%       tp_hr   : time to peak of the hydrograph [h]
%
%   OUTPUTS
%       Qin      : inflow discharge hydrograph [m^3/s] at t_hours
%       Qp_exact : analytic peak discharge [m^3/s], i.e.
%                  Qp_exact = b_exact * V_run / (tp_hr * 3600)
%       b_exact  : exact peaking factor for given m_shape:
%                  b_exact = (m-1)^m / (Gamma(m)*exp(m-1))
%
%   Construction (continuous-time theory):
%       • Gamma/Nash unit hydrograph u(t) [1/h] with shape m and scale θ:
%             u(t) = t^(m-1) * exp(-t/θ) / (Γ(m) * θ^m),    t >= 0
%         which satisfies ∫_0^∞ u(t) dt = 1.
%       • Choose θ so that t_p = (m - 1) θ  ⇒  θ = t_p / (m - 1).
%       • Then the continuous-time peak of u(t) is:
%             u_peak = b(m) / t_p,
%         where b(m) = (m-1)^m / (Γ(m) * e^(m-1)).
%       • Discharge hydrograph: Q_hr(t) = V_run * u(t)  [m^3/h],
%         so Qp = V_run * u_peak = b(m) * V_run / t_p  [m^3/h].
%       • Convert to m^3/s: divide by 3600.

    % Ensure column vector
    t_hours = t_hours(:);

    % Basic check on time vector
    if numel(t_hours) < 2
        error('t_hours must have at least 2 points.');
    end

    % Gamma scale parameter θ [h] such that t_p = (m - 1) θ
    theta_hr = tp_hr / (m_shape - 1);

    % Nash unit hydrograph u(t) [1/h] (exact gamma pdf, no renorm)
    u = zeros(size(t_hours));
    positive_idx = t_hours >= 0;       % include t = 0; u(0) = 0 for m > 1
    t_pos = t_hours(positive_idx);

    % Gamma PDF: u(t) = t^(m-1) exp(-t/θ) / (Γ(m) θ^m)
    u(positive_idx) = (t_pos.^(m_shape - 1) .* exp(-t_pos ./ theta_hr)) ./ ...
                      (gamma(m_shape) * theta_hr^m_shape);
    % NOTE:
    %   • This u(t) is already normalized in the continuous sense:
    %       ∫_0^∞ u(t) dt = 1  (exact).
    %   • We intentionally do NOT apply any discrete renormalization here.

    % Discharge in m^3/h
    Q_hr = V_run * u;   % continuous integral over 0..∞ is exactly V_run

    % Convert to m^3/s
    Qin = Q_hr / 3600;  % [m^3/s]

    % --- Analytic peaking factor and peak flow (for reference) ----------
    % b(m) = (m-1)^m / (Γ(m) * e^(m-1))
    b_exact = (m_shape - 1)^m_shape ./ (gamma(m_shape) * exp(m_shape - 1));

    % Analytic continuous-time peak discharge [m^3/s]:
    % Qp = (b * V_run) / t_p, with t_p in hours → convert to seconds.
    Qp_exact = (b_exact * V_run) / (tp_hr * 3600);

end


function [eta_p, Delta_tp_min, DetTime_min, Qp_out, DeltaV_m3] = compute_metrics(t_minutes, Qin, Qout)
%COMPUTE_METRICS Compute hydrological performance metrics for DRAIN-LID.
%
% Inputs:
%   t_minutes : time vector [min], same size as Qin/Qout
%   Qin       : inflow discharge [m^3/s]
%   Qout      : outflow discharge [m^3/s]
%
% Outputs:
%   eta_p        : peak flow reduction [%]
%   Delta_tp_min : time-to-peak mitigation [min] (tp_out - tp_in)
%   DetTime_min  : detention time [min]
%   Qp_out       : peak outflow [m^3/s]
%   DeltaV_m3    : runoff volume difference [m^3] = V_in - V_out
%
% Definitions:
%   • Peak flow reduction:
%       η_p = 100 * (Qp_in - Qp_out) / Qp_in
%
%   • Time-to--peak mitigation:
%       Δt_p = t_p^o - t_p^i
%
%   • Detention time:
%       Δt_d = t(Q_out(t) <= τ) - t_p^i,
%       where τ = 0.25% * Qp_in = 0.0025 * Qp_in
%
%   • Runoff difference:
%       ΔV = ∫ Qin dt - ∫ Qout dt   [m^3]

    % --- Ensure column vectors ------------------------------------------
    t    = t_minutes(:);
    Qin  = Qin(:);
    Qout = Qout(:);

    % --- Basic checks ---------------------------------------------------
    if numel(t) ~= numel(Qin) || numel(t) ~= numel(Qout)
        error('compute_metrics:InputSizeMismatch', ...
              't_minutes, Qin, and Qout must have the same length.');
    end

    % Handle NaNs or empty vectors
    if isempty(Qin) || isempty(Qout) || all(isnan(Qin)) || all(isnan(Qout))
        eta_p        = NaN;
        Delta_tp_min = NaN;
        DetTime_min  = NaN;
        Qp_out       = NaN;
        DeltaV_m3    = NaN;
        return;
    end

    % --- 1. Peak flows and times ----------------------------------------
    [Qp_in,  idx_in]  = max(Qin);   % inflow peak
    [Qp_out, idx_out] = max(Qout);  % outflow peak

    tp_in_min  = t(idx_in);         % [min]
    tp_out_min = t(idx_out);        % [min]

    % --- 2. Peak flow reduction η_p [%] ---------------------------------
    if Qp_in > 0
        eta_p = 100 * (Qp_in - Qp_out) / Qp_in;
    else
        eta_p = NaN;  % undefined if no inflow
    end

    % --- 3. Time-to-peak mitigation Δt_p [min] --------------------------
    Delta_tp_min = tp_out_min - tp_in_min;

    % --- 4. Detention time Δt_d [min] -----------------------------------
    % Threshold τ = 0.25% of inflow peak (0.0025 * Qp_in)
    tau = 0.0025 * Qp_in;

    % Consider only times after inflow peak
    idx_after_in = find(t >= tp_in_min);

    if isempty(idx_after_in) || Qp_in <= 0
        DetTime_min = NaN;
    else
        % Find last time where outflow is above threshold
        idx_above = idx_after_in(Qout(idx_after_in) > tau);

        if isempty(idx_above)
            % Outflow never exceeds threshold after inflow peak
            DetTime_min = 0;
        else
            idx_last = idx_above(end);
            t_last   = t(idx_last);
            DetTime_min = t_last - tp_in_min;   % [min]
        end
    end

    % --- 5. Runoff volume difference ΔV [m^3] ---------------------------
    % Convert t from minutes to seconds for integration
    t_sec = t * 60;   % [s]

    % Volumes [m^3] using trapezoidal rule
    V_in  = trapz(t_sec, Qin);
    V_out = trapz(t_sec, Qout);

    if V_in > 0
        DeltaV_m3 = V_in - V_out;   % [m^3]
    else
        DeltaV_m3 = NaN;
    end

end


