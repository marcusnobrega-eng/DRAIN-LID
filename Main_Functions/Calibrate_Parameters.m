% Example: Assume 3 layers and you want to calibrate alpha, n, and Ks per layer
% Parameter vector: [alpha1, alpha2, alpha3, n1, n2, n3, Ks1, Ks2, Ks3]

% Lower bounds:
%       α (1/m)     n (-)     θr (-)    θs (-)    Ss (1/m)    Ks (m/s)    Ss (1/m) 
lb = [ 3,        3,       3,      ...    % α parameters (1/m)
       1.25,        1.25,       1.25,      ...    % n parameters (-)
       0.005,        0.005,       0.005,      ...    % θr parameters (-)
       0.20,       0.20,      0.20,     ...    % θs parameters (-)
       1e-5,       1e-5,      1e-5,     ...    % Ss parameters (1/m)
       1e-5,       1e-5,      1e-5 ];          % Ks parameters (m/s)

% Upper bounds:
%       α (1/m)     n (-)     θr (-)    θs (-)    Ss (1/m)    Ks (m/s)    Ss (1/m)
ub = [ 15.0,     15.0,    15.0,  ...    % α parameters (1/m)
       3.0,        3.0,       3.0,     ...    % n parameters (-)
       0.005,       0.005,      0.005,    ...    % θr parameters (-)
       0.50,       0.50,      0.50,    ...    % θs parameters (-)
       1e-5,       1e-5,      1e-5,    ...    % Ss parameters (1/m)
       1e-2,       1e-2,      1e-2];         % Ks parameters (m/s)

% Number of parameters
n = length(lb);

% Read observed data
obs_data = readtable('C:\Users\marcu\Documents\GitHub\LID_Tool\Forcing_Event\Observed_Data.xlsx');
obs_data = table2array(obs_data(:,2))';

% Baseline path
baseline_path = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Examples\Monitored_PP_Events_Data.mat';

options = optimoptions('ga', ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'MaxGenerations', 10, ...
    'PopulationSize', 40, ...
    'PlotFcn', {@gaplotbestf}); % This will plot best fitness at each generation

[x_opt, fval] = ga(@(x) calibration_objective(x, obs_data, baseline_path), n, [], [], [], [], lb, ub, [], options);

disp('Optimal Parameters:');
disp(x_opt);


function error = calibration_objective(x, observed_discharge, baseline_path)

    % Run simulation
    % x = [3.38327384068837	4.05346527787615	9.18686475822007	2.00464342959708	1.72616688255844	2.92220754843722	0.00500000000000000	0.00500000000000000	0.00500000000000000	0.235274277694258	0.200576615766029	0.410104712387117	1.00000000000000e-05	1.00000000000000e-05	1.00000000000000e-05	0.00831705588184255	0.00652553064613185	0.00738990250351285];
    try
        [~, sim_output] = run_richards_model(x,baseline_path);
    catch
        % If simulation crashes, penalize the solution
        error = 1e6;
        return;
    end
    
    if isnan(sim_output.flux_out)
        error = 1e6;
        return
    end
    % Extract simulated discharge at the bottom node [m/s]
    sim_flux = sim_output.flux_out(1, :);  % Bottom flux over time

    % Convert simulated flux from [m/s] to [mm/h]
    sim_discharge = (-1) * sim_flux * 1000 * 3600;

    % Start calibration after warm up
    sim_discharge(1:1440/15) = [];
    observed_discharge(1:1440/15) = [];
    % Compute NSE
    if(~isnan(sim_discharge))
        numerator = sum((sim_discharge - observed_discharge).^2);
        denominator = sum((observed_discharge - mean(observed_discharge)).^2);
    
        if denominator == 0
            % Constant observed discharge, penalize solution
            error = 1e6;
        else
            NSE = 1 - (numerator / denominator);
            error = 1 - NSE;  % Minimize (1 - NSE) → Maximize NSE
        end
    else
        NSE = -10000;
        error = 1 - NSE; % Minimize (1 - NSE) → Maximize NSE
    end

end


function [sim_time, sim_output] = run_richards_model(x,baseline_path)
% =========================================================================
% 🌱 Wrapper to Run Richards Model for Calibration Purposes
% =========================================================================
% Inputs:
%   - x: Struct containing all simulation and soil parameters
%
% Outputs:
%   - sim_time: Time vector from the simulation
%   - sim_output: Struct containing key simulation outputs:
%       • ponding_series
%       • head_out
%       • flux_out
%
% Notes:
%   - This function is designed for automatic calibration.
%   - All plotting is disabled to maximize performance.
%   - Assumes solver is modularized as 'Main_Solver_Function.m'.
%
% Author: Marcus Nóbrega, Ph.D.
% =========================================================================

    try
        % 🛠 Run the solver with the provided parameters
        [time_series, ponding_series, head_out, flux_out] = Main_Solver_Function(x,baseline_path);

        % 📦 Package outputs
        sim_time = time_series;
        sim_output.ponding_series = ponding_series;
        sim_output.head_out = head_out;
        sim_output.flux_out = flux_out;

    catch ME
        % 🚨 Handle any solver crash (e.g., divergence) to avoid stopping GA
        warning('Solver failed during run')

        % Return penalty outputs
        sim_time = [];
        sim_output.ponding_series = NaN;
        sim_output.head_out = NaN;
        sim_output.flux_out = NaN;
    end

end
