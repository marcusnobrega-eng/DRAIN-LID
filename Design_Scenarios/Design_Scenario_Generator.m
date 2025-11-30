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
Peff_set_mm = [10, 25, 50, 75, 100];   % [mm] over upstream catchment

% Contributing-area ratio (upstream : bioretention)
area_ratio  = 10.0;                    % A_up / A_TC

% Time to peak of inflow hydrograph (minutes)
tp_set_min  = [15, 30, 60];            % [min]

% Bioretention geometry
ATC_set_m2  = [5, 10, 20, 40, 80];     % [m^2] footprint
Ld_set_m    = [0.4, 0.6, 0.8, 1.0, 1.2]; % [m] media depth

% Soil hydraulic parameter sets (van Genuchten–Mualem)
% TODO: replace with your actual calibrated table (from HYDRUS / Carsel & Parrish)
Soils = define_soils();   % returns struct array with fields:
                          % name, theta_r, theta_s, alpha_1perm, n, Ks_mm_per_h

nSoil  = numel(Soils);
nATC   = numel(ATC_set_m2);
nLd    = numel(Ld_set_m);
nPeff  = numel(Peff_set_mm);
nTp    = numel(tp_set_min);

% Ponding constraint
h_max_allow = 0.15;       % [m] = 15 cm

% Nash hydrograph parameters
m_shape     = 4.69;       % shape factor
b_factor    = 0.75;       % peaking factor (as in your derivation)

% Simulation time settings
T_hours     = 36;         % total simulation duration [h]
dt_min      = 5;          % output time step for hydrograph [min]
t_minutes   = (0:dt_min:(T_hours*60))';      % [min]
t_hours     = t_minutes / 60;               % [h]
nSteps      = numel(t_hours);

%% 2. Preallocate storage

% Total number of combinations (before feasibility filtering)
N_total = nSoil * nATC * nLd * nPeff * nTp;

% We don't know a priori how many will be feasible, so we can either:
% (a) preallocate at N_total and keep a counter, or
% (b) collect in a growing struct/cell array.
% Here we preallocate and trim at the end for simplicity.

Metrics = struct( ...
    'soilName',      cell(N_total,1), ...
    'ATC_m2',        nan(N_total,1), ...
    'Ld_m',          nan(N_total,1), ...
    'Peff_catch_mm', nan(N_total,1), ...
    'Peff_eq_mm',    nan(N_total,1), ...
    'tp_min',        nan(N_total,1), ...
    'Qp_in_m3s',     nan(N_total,1), ...
    'Qp_out_m3s',    nan(N_total,1), ...
    'eta_p_pct',     nan(N_total,1), ...
    'Delta_tp_min',  nan(N_total,1), ...
    'DetTime_min',   nan(N_total,1), ...
    'DeltaV_m3',     nan(N_total,1));   % 🔹 new field

% Total number of combinations (before feasibility filtering)
N_total = nSoil * nATC * nLd * nPeff * nTp;

% Progress counter (attempted scenarios)
combo_idx = 0;

% Time series: store Q_in and Q_out at 5-min resolution
Qin_all  = nan(nSteps, N_total);   % [m^3/s]
Qout_all = nan(nSteps, N_total);   % [m^3/s]

k = 0;   % feasible scenario counter

%% 3. Main nested loops

for iSoil = 1:nSoil
    soil = Soils(iSoil);

    for iATC = 1:nATC
        ATC = ATC_set_m2(iATC);  % [m^2]

        for iLd = 1:nLd
            Ld = Ld_set_m(iLd);  % [m]

            for iPeff = 1:nPeff
                Peff_catch_mm = Peff_set_mm(iPeff); % [mm] over A_up

                % Equivalent effective depth over bioretention footprint [mm]
                Peff_eq_mm = Peff_catch_mm * area_ratio;  % [mm over A_TC]
                Peff_eq_m  = Peff_eq_mm / 1000;           % [m]

                % Runoff volume over bioretention surface [m^3]
                V_run = Peff_eq_m * ATC;

                for iTp = 1:nTp
                    tp_min = tp_set_min(iTp);     % [min]
                    tp_hr  = tp_min / 60;         % [h]
                
                    % ==== 🔄 Progress tracking ==========================================
                    combo_idx = combo_idx + 1;
                    if mod(combo_idx, 50) == 0 || combo_idx == 1 || combo_idx == N_total
                        progress = 100 * combo_idx / N_total;
                        fprintf('Progress: %6.2f %%  (%d / %d scenarios attempted)\n', ...
                                progress, combo_idx, N_total);
                    end
                    % ====================================================================
                
                    % Generate inflow hydrograph at 5-min resolution using Nash UH
                    Qin = generate_nash_hydrograph( ...
                              t_hours, V_run, m_shape, tp_hr);
                
                    % 🔹 Define inflow peak (used later in Metrics)
                    Qp_in = max(Qin);   % [m^3/s]
                                    
                    % --- Call DRAIN-LID Richards solver (user-defined) ---
                    [Qout, h_surf] = run_drain_lid( ...
                        t_hours, Qin, ATC, Ld, soil, base_params);
                    
                    % If solver failed or returned invalid outputs, skip this design
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

                    % Compute metrics
                    [eta_p, Delta_tp_min, DetTime_min, Qp_out, DeltaV_m3] = ...
                        compute_metrics(t_minutes, Qin, Qout);

                    % Fill metrics struct
                    Metrics(k).soilName      = soil.name;
                    Metrics(k).ATC_m2        = ATC;
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

                end % iTp
            end % iPeff
        end % iLd
    end % iATC
end % iSoil

%% 4. Trim unused preallocated entries

if k < N_total
    Metrics   = Metrics(1:k);
    Qin_all   = Qin_all(:, 1:k);
    Qout_all  = Qout_all(:, 1:k);
end

%% 5. Save database

save('DRAIN_LID_database.mat', ...
     'Metrics', 'Qin_all', 'Qout_all', 't_minutes', 'Soils', ...
     'Peff_set_mm', 'tp_set_min', 'ATC_set_m2', 'Ld_set_m', ...
     'area_ratio', 'h_max_allow', 'm_shape', 'b_factor');

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

function [Qout, h_surf, t_save_hours] = run_drain_lid(t_hours, Qin, ATC, Ld, soil, base_params)
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
params.h_spill          = 0.15;  % [m] maximum ponding depth
params.exp_spillway     = 1.5;

% Spillway coefficient (can be tuned)
params.c_spillway       = 1.0;   % placeholder

% -------------------------------------------------------------------------
% 6. Initial conditions: uniform suction profile
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 6. Initial conditions: field capacity-based suction
% -------------------------------------------------------------------------
% Following standard soil physics and USDA practice:
%   • Coarse-textured soils (Sand, Loamy Sand) → FC ≈ θ(h = -10 kPa) ≈ h = -1.0 m
%   • Medium/fine-textured soils (others)     → FC ≈ θ(h = -33 kPa) ≈ h = -3.3 m
%
% We approximate field capacity by assigning a uniform matric head h_FC
% and letting the van Genuchten curve provide the corresponding θ(h_FC).

coarse_names = {'Sand','Loamy Sand'};

if any(strcmpi(soil.name, coarse_names))
    h_FC = -1.0;   % [m] ≈ -10 kPa for coarse soils
else
    h_FC = -3.3;   % [m] ≈ -33 kPa for medium/fine soils
end

params.initial_suction = h_FC;   % store for bookkeeping (if used elsewhere)

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
[head_out, theta_out, flux_out, ponding_series, outlet_flux, time_series, success] = Main_Solver_Function(params);

% If solver failed, return empty so the outer loop can skip this design
if ~success || isempty(time_series) || isempty(outlet_flux)
    Qout   = [];
    h_surf = [];
    t_save_hours = [];
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

function Qin = generate_nash_hydrograph(t_hours, V_run, m_shape, tp_hr)
%GENERATE_NASH_HYDROGRAPH Generates a Nash-type inflow hydrograph
%   directly on the prescribed time grid t_hours (in hours).
%
%   Qin = generate_nash_hydrograph(t_hours, V_run, m_shape, tp_hr)
%
%   INPUTS
%       t_hours : time vector [h], e.g. 0:5/60:36
%       V_run   : total runoff volume reaching the bioretention [m^3]
%       m_shape : Nash shape parameter m [-]
%       tp_hr   : time to peak of the hydrograph [h]
%
%   OUTPUT
%       Qin     : inflow discharge hydrograph [m^3/s]
%
%   Construction:
%       • Define a Gamma/Nash unit hydrograph u(t) [1/h] with shape m
%         and scale θ such that t_p = (m - 1) θ.
%       • Normalize u(t) discretely so that Σ u(t_i) Δt = 1.
%       • Set Q_hr(t) = V_run * u(t)  [m^3/h].
%       • Convert to Q_s(t) = Q_hr(t) / 3600  [m^3/s].

    % Ensure column vector
    t_hours = t_hours(:);

    % Time step in hours (assumed uniform, e.g. 5 min = 1/12 h)
    if numel(t_hours) < 2
        error('t_hours must have at least 2 points.');
    end
    dt_hr = t_hours(2) - t_hours(1);

    % Gamma scale parameter θ [h] such that t_p = (m - 1) θ
    theta_hr = tp_hr / (m_shape - 1);

    % Nash unit hydrograph u(t) [1/h]
    u = zeros(size(t_hours));
    positive_idx = t_hours > 0;
    t_pos = t_hours(positive_idx);

    % Gamma PDF: u(t) = t^(m-1) exp(-t/θ) / (Γ(m) θ^m)
    u(positive_idx) = (t_pos.^(m_shape - 1) .* exp(-t_pos ./ theta_hr)) ./ ...
                      (gamma(m_shape) * theta_hr^m_shape);

    % Discrete normalization so that Σ u Δt = 1
    area_u = sum(u) * dt_hr;
    if area_u > 0
        u = u / area_u;
    end

    % Discharge in m^3/h
    Q_hr = V_run * u;   % integral over hours ≈ V_run

    % Convert to m^3/s
    Qin = Q_hr / 3600;  % [m^3/s]

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


