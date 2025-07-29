%% =========================================================================
% 🌧️ Catchment_Outputs.m
% -------------------------------------------------------------------------
% Purpose   : Run hydrological and pollutant generation model for a small
%             watershed and prepare top boundary condition (Neumann) for 
%             Richards Equation infiltration model.
% Inputs    :
%   - dt                      : time-step [sec]
%   - LID_area                : Receiving drainage LID area [m2]
% Outputs   :
%   - time_watershed_seconds  : time vector [sec]
%   - qin                     : flow rate [m/s]
%   - Cpol                    : pollutant concentration [mg/L] (if enabled)
% -------------------------------------------------------------------------
% Author    : Marcus Nóbrega, Ph.D.
% Updated   : May 2025
%% =========================================================================

function [time_watershed_seconds, qin, Cpol] = Catchment_Outputs(dt, LID_area, forcing_path)

%% === 1. Watershed Physical Properties ===================================
width_w      = 1;   % Watershed width [m]
length_w     = 1;   % Watershed length [m]
Watershed_Area = width_w * length_w;  % Total area [m²]

%% === 2. SCS Curve Number Model Setup ===================================
CN_per        = 90;     % Curve number (pervious)
h0_w          = 0.01;   % Initial abstraction [m]
n_w           = 0.02;   % Manning's n [-]
Aimp          = 0;      % Impervious fraction [-]
slope_w       = 0.015;  % Slope [m/m]
baseflow_w    = 0;      % Baseflow [m³/s]
ks_w          = 50;     % Saturated K [mm/h]
kr_w          = 5;      % Recovery rate [mm/h]
A_GI          = 0;      % GI area ratio [-]
CN_GI         = 65;     % GI Curve Number [-]
catch_GI_imp  = 0;      % % of impervious areas draining to GI [-]

%% === 3. Read Input Spreadsheet ==========================================
input_table = readtable(forcing_path);
flag_manual_hydrograph = table2array(input_table(1,6));

%% === 4. Run Hydrologic Model ============================================
if ~flag_manual_hydrograph
    tfinal = table2array(input_table(end,1));  % in minutes
    [time_watershed, Qin, ~, ~, time_inflow_hydrograph, inflow_hydrograph_data] = ...
        SCS_Hydrologic_Model(dt, tfinal, width_w, length_w, CN_per, h0_w, ...
        n_w, Aimp, slope_w, baseflow_w, ks_w, kr_w, A_GI, CN_GI, catch_GI_imp);

elseif flag_manual_hydrograph 
    inflow_hydrograph_table = table2array(input_table(2:end, 1:4));
    time_inflow_hydrograph = inflow_hydrograph_table(:,1);  % [min]
    inflow_hydrograph_data = inflow_hydrograph_table(:,4);  % [mm/h]
    Qin = inflow_hydrograph_data / 1000 / 3600 * Watershed_Area;  % m³/s
    time_watershed = time_inflow_hydrograph;
else
    error('Hydrograph flag not enabled or improperly configured.');
end

time_watershed_seconds = time_watershed * 60; % [sec]

%% === 5. Prepare Top Boundary (Flux) for Richards ========================
if flag_manual_hydrograph
    qin = inflow_hydrograph_data / 1000 / 3600;  % Already with the convention signal
else 
    qin = -Qin' / LID_area;           % [m/s]
end

%% === 6. Water Quality Model Parameters ===================================
    % Buildup and Washoff parameters
    C1 = 50;                         % Buildup coefficient in (kg/ha)
    C2 = 0.3;                        % Buildup exponent in (1/days)
    C3 = 0.02;                       % Washoff coefficient in (mm/h)^(-C4)(h^-1)
    C4 = 1.5;                        % Washoff exponent (1/h)
    ADD = 10;                        % Number of antecedent dry days   
    Area_km2 = Watershed_Area/1000/1000; % Converting to km2

    % Run Buildup and Washoff Model
    tint = time_watershed;
    [~, Cpol, ~, ~, ~, ~, ~, ~] = ...
        Buildup_Washoff_Model(C1, C2, C3, C4, ADD, Area_km2, Qin, tint);
    close all
end

