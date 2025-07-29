%% ============================================================================
%                          📥 INPUT AND SETUP
% ============================================================================
% File paths
close all; clear all; clc;
xlsPath   = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example4_TopNeumann_Sandy\Data\SimulationResults.xlsx';
hydrusOut = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example4_TopNeumann_Sandy\Data\Example4_TopNeumann_Sandy.out';
pdfOutput = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example4_TopNeumann_Sandy\Data\Example_4_Comparison_Figure.pdf';

% Load model parameters
load Examples\Example4_TopNeumann_Sandy.mat

% Load model data
modelData = readtable(xlsPath, 'Sheet', 'Head');
time_s = modelData.Time_s;
model_heads = table2array(modelData(:, startsWith(modelData.Properties.VariableNames, 'h_Node')));  % [nt x nz]

% Load model moisture data
theta_table = readtable(xlsPath, 'Sheet', 'Moisture');
model_theta = table2array(theta_table(:, startsWith(theta_table.Properties.VariableNames, 'theta_Node')));

%% ============================================================================
%                    💧 PARSE HYDRUS .out FILE FOR HEAD & THETA
% ============================================================================
fid = fopen(hydrusOut, 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
lines = lines{1};

hydrus_heads_blocks = {};
hydrus_theta_blocks = {};
hydrus_times = [];

current_heads = [];
current_thetas = [];

for i = 1:length(lines)
    line = strtrim(lines{i});
    
    if startsWith(line, 'Time:')
        if ~isempty(current_heads)
            hydrus_heads_blocks{end+1} = current_heads;
            hydrus_theta_blocks{end+1} = current_thetas;
        end
        current_heads = [];
        current_thetas = [];
        t_val = sscanf(line, 'Time: %f');
        hydrus_times(end+1) = t_val;

    elseif ~isempty(regexp(line, '^\d+', 'once')) && ~contains(line, 'Node')
        tokens = sscanf(line, '%f');
        if length(tokens) >= 4
            current_heads(end+1) = tokens(3);   % Column 3 = Pressure head
            current_thetas(end+1) = tokens(4);  % Column 4 = Moisture
        end
    end
end

% Append last block
if ~isempty(current_heads)
    hydrus_heads_blocks{end+1} = current_heads;
    hydrus_theta_blocks{end+1} = current_thetas;
end

hydrus_heads = cell2mat(hydrus_heads_blocks');         % [nt x nz]
hydrus_theta_matrix = cell2mat(hydrus_theta_blocks');  % [nt x nz]
hydrus_times = hydrus_times(:);

%% ============================================================================
%                    📈 INTERPOLATE TO MODEL TIMES
% ============================================================================
hydrus_interp = interp1(hydrus_times, hydrus_heads(:, 1:size(model_heads,2)), time_s, 'linear', 'extrap');
hydrus_theta_interp = interp1(hydrus_times, hydrus_theta_matrix(:, 1:size(model_heads,2)), time_s, 'linear', 'extrap');

%% ============================================================================
%                    📊 PERFORMANCE METRICS FUNCTIONS
% ============================================================================
nse = @(sim, obs) 1 - sum((sim - obs).^2) / sum((obs - mean(obs)).^2);
rmse = @(sim, obs) sqrt(mean((sim - obs).^2));
pbias = @(sim, obs) 100 * sum(sim - obs) / sum(obs);
kge = @(sim, obs) ...
    (1 - sqrt((corr(sim', obs') - 1).^2 + ...
    (std(sim)/std(obs) - 1).^2 + ...
    (mean(sim)/mean(obs) - 1).^2));

%% ============================================================================
%                  🧮 COMPUTE METRICS FOR PRESSURE HEAD
% ============================================================================
nSteps = length(time_s);
NSE = nan(nSteps, 1); RMSE = nan(nSteps, 1);
PBIAS = nan(nSteps, 1); KGE = nan(nSteps, 1);

for i = 2:(nSteps - 1)
    sim = model_heads(i, :);
    obs = flip(hydrus_interp(i, :));
    if any(isnan(sim)) || any(isnan(obs)) || all(obs == obs(1)), continue; end
    NSE(i) = nse(sim, obs); RMSE(i) = rmse(sim, obs);
    PBIAS(i) = pbias(sim, obs); KGE(i) = kge(sim, obs);
end

%% ============================================================================
%                🧮 COMPUTE METRICS FOR MOISTURE CONTENT
% ============================================================================
NSE_theta = nan(nSteps, 1); RMSE_theta = nan(nSteps, 1);
PBIAS_theta = nan(nSteps, 1); KGE_theta = nan(nSteps, 1);

for i = 2:(nSteps - 1)
    sim_theta = model_theta(i, :);
    obs_theta = flip(hydrus_theta_interp(i, :));
    if any(isnan(sim_theta)) || any(isnan(obs_theta)) || all(obs_theta == obs_theta(1)), continue; end
    NSE_theta(i) = nse(sim_theta, obs_theta); RMSE_theta(i) = rmse(sim_theta, obs_theta);
    PBIAS_theta(i) = pbias(sim_theta, obs_theta); KGE_theta(i) = kge(sim_theta, obs_theta);
end

%% ============================================================================
%                      🎨 PLOTTING SETTINGS AND DEPTH SETUP
% ============================================================================
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();  % Assuming you have this function
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontName','Montserrat');
set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesLineWidth',1.5);

key_fracs = [0.1, 0.25, 0.5, 0.75, 1.0];
key_indices = round(key_fracs * length(time_s));
z = linspace(-params.L, 0, size(model_heads,2));  % Depth

%% ============================================================================
%                         🎯 GENERATE 2x2 PLOT FIGURE
% ============================================================================
fig = figure('Color','w', 'Units','inches', 'Position', [1, 1, 8.27, 11.69]);  % A4 portrait

% --- RMSE Pressure Head
subplot(2,2,1);
plot(time_s/3600, RMSE, '-', 'Color', pallete.red_colors(1,:), 'LineWidth', 2.5);
ylabel('\textbf{RMSE ($h$) [m]}'); xlabel('\textbf{Time [hours]}');
title('\textbf{Performance of $h$}');
grid on; box on; set(gca, 'TickDir', 'out');
ytickformat('%.3f');               % Use 3 decimal places
ax = gca;
ax.YAxis.Exponent = 0;             % Disable automatic scientific notation

% --- RMSE Moisture
subplot(2,2,3);
plot(time_s/3600, RMSE_theta, '-', 'Color', pallete.green_colors(1,:), 'LineWidth', 2.5);
ylabel('\textbf{RMSE ($\theta$) [$\mathrm{cm^3 \cdot cm^{-3}}$]}'); xlabel('\textbf{Time [hours]}');
title('\textbf{Performance of $\theta$}');
grid on; box on; set(gca, 'TickDir', 'out');
ytickformat('%.3f');               % Use 3 decimal places
ax = gca;
ax.YAxis.Exponent = 0;             % Disable automatic scientific notation
% --- Pressure Head Profiles
subplot(2,2,2); hold on;
clear lines
colors = lines(length(key_indices));
for k = 1:length(key_indices)
    idx = key_indices(k);
    plot(flip(hydrus_interp(idx, :)), z, '-', 'LineWidth', 2.5, 'Color', colors(k,:));
    plot(model_heads(idx, :), z, ':', 'LineWidth', 2.5, 'Color', colors(k,:));
end
xlabel('\textbf{Pressure Head [m]}'); ylabel('\textbf{Depth [m]}');
title('\textbf{Pressure Head Profiles}');
grid on; box on; set(gca, 'TickDir', 'out');

% --- Moisture Profiles
subplot(2,2,4); hold on;
for k = 1:length(key_indices)
    idx = key_indices(k);
    plot(flip(hydrus_theta_interp(idx, :)), z, '-', 'LineWidth', 2.5, 'Color', colors(k,:));
    plot(model_theta(idx, :), z, ':', 'LineWidth', 2.5, 'Color', colors(k,:));
end
xlabel('\textbf{Volumetric Moisture $\theta$ [m$^3$/m$^3$]}'); ylabel('\textbf{Depth [m]}');
title('\textbf{Moisture Profiles}');
grid on; box on; set(gca, 'TickDir', 'out');

%% ============================================================================
%                           💾 EXPORT FIGURE TO PDF
% ============================================================================
% Export the figure cleanly
% === Force Montserrat font on all text elements (especially tick labels) ===
set(findall(fig, '-property', 'FontName'), 'FontName', 'Montserrat');
set(findall(fig, '-property', 'FontUnits'), 'FontUnits', 'points');  % Ensure scalable
set(findall(fig, '-property', 'FontSizeMode'), 'FontSizeMode', 'manual');  % Prevent resizing

exportgraphics(fig, pdfOutput, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', 300);

