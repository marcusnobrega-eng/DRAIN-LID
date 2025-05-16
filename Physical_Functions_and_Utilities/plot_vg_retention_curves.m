function plot_vg_retention_curves(layer_properties, params, directory)
% PLOT_VG_RETENTION_CURVES - Plots Van Genuchten-Mualem hydraulic functions for porous media.
%
% This function generates three diagnostic plots for layered unsaturated media:
%   1. Water retention curve:        θ(h)
%   2. Unsaturated conductivity:     K(h)
%   3. Conductivity vs water content: K(θ)
%
% The function is based on the Van Genuchten (1980) model for water retention
% and the Mualem (1976) model for hydraulic conductivity.
%
% INPUTS:
%   layer_properties - Cell array with each row:
%       { 'Layer Name', alpha [1/cm], theta_s, theta_r, n }
%           - 'Layer Name' : string name for the soil layer
%           - alpha        : inverse air-entry pressure [1/cm]
%           - theta_s      : saturated volumetric water content [cm³/cm³]
%           - theta_r      : residual volumetric water content [cm³/cm³]
%           - n            : Van Genuchten shape parameter [-]
%
%   Ks_values - Vector of saturated hydraulic conductivities [cm/s]
%
% EXAMPLE:
%   layers = {
%       'Gravel subbase ASTM No. 8', 0.0045, 0.40, 0.02, 2.20;
%       'Gravel base ASTM #2',       0.05,   0.40, 0.01, 2.30;
%       'Permeable Asphalt (PA)',    0.02,   0.25, 0.05, 2.10
%   };
%   Ks_values = [36, 36, 10];  % [cm/s]
%   plot_vg_retention_curves(layers, Ks_values)

%% --- Define suction pressure head range [cm] ---
% Suction head varies logarithmically from 0.01 to 1000 cm
h = logspace(-2, 3, 400);  % suction head [cm]

%% --- Define Van Genuchten & Mualem functions ---
% θ(h): Van Genuchten water retention equation
%   θ(h) = θ_r + (θ_s - θ_r) / (1 + (αh)^n)^(1 - 1/n)
theta_vg = @(h, theta_s, theta_r, alpha, n) ...
    theta_r + (theta_s - theta_r) ./ ((1 + (alpha .* h).^n).^(1 - 1./n));

% Se(θ): Effective saturation
%   Se = (θ - θ_r) / (θ_s - θ_r)
Se_vg = @(theta, theta_s, theta_r) ...
    (theta - theta_r) ./ (theta_s - theta_r);

% K(Se): Mualem hydraulic conductivity model
%   K(Se) = Ks * Se^0.5 * [1 - (1 - Se^(1 - 1/n))^n]^2
K_mualem = @(Se, Ks, n) ...
    Ks .* (Se).^0.5 .* (1 - (1 - Se.^(1 - 1./n)).^n).^2;

%% --- Set up figure for plots ---
fig = figure('Color','w');
set(fig, 'Units', 'Inches', 'Position', [1, 1, 11.5, 6]);
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Plot 1: θ(h)
nexttile(1); hold on; box on; grid on;
set(gca, 'XScale', 'log');
xlabel('Pressure head $h$ [cm]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\theta$ [cm$^3$/cm$^3$]', 'Interpreter', 'latex', 'FontSize', 14);
title('$\theta(h)$', 'Interpreter', 'latex', 'FontSize', 16);
axis tight

% Plot 2: K(h)
nexttile(2); hold on; box on; grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Pressure head $h$ [cm]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$K(h)$ [cm/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('$K(h)$', 'Interpreter', 'latex', 'FontSize', 16);
axis tight

% Plot 3: K(θ)
nexttile(3); hold on; box on; grid on;
set(gca, 'YScale', 'log');
xlabel('$\theta$ [cm$^3$/cm$^3$]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$K(\theta)$ [cm/s]', 'Interpreter', 'latex', 'FontSize', 14);
title('$K(\theta)$', 'Interpreter', 'latex', 'FontSize', 16);
axis tight

%% --- Define color scheme ---
colors = lines(size(layer_properties.alpha,2));  % consistent colors per layer

%% --- Loop through each soil layer and compute curves ---
for i = 1:size(layer_properties.alpha,2)
    % Extract layer parameters
    label = layer_properties.labels{i};
    alpha = layer_properties.alpha(i);
    theta_s = layer_properties.theta_s(i);
    theta_r = layer_properties.theta_r(i);
    n = layer_properties.n(i);
    Ks = params.Ks(i);

    % Compute volumetric water content from θ(h)
    theta = theta_vg(h, theta_s, theta_r, alpha, n);

    % Compute effective saturation Se(h)
    Se = Se_vg(theta, theta_s, theta_r);
    Se = min(max(Se, 0), 1);  % constrain Se to [0,1]

    % Compute unsaturated hydraulic conductivity K(h)
    K = K_mualem(Se, Ks, n);

    % Plot θ(h)
    nexttile(1);
    plot(-h, theta, '-', 'LineWidth', 2, 'DisplayName', label, 'Color', colors(i,:));
    set(gca,'XDir','reverse')  % suction heads are negative

    % Plot K(h)
    nexttile(2);
    plot(-h, K, '-', 'LineWidth', 2, 'DisplayName', label, 'Color', colors(i,:));
    set(gca,'XDir','reverse')  % suction heads are negative

    % Plot K(θ)
    nexttile(3);
    plot(theta, K, '-', 'LineWidth', 2, 'DisplayName', label, 'Color', colors(i,:));
end

%% --- Legends and final formatting ---
nexttile(1); legend('Location','southwest', 'Interpreter','latex', 'FontSize', 11);
nexttile(2); legend('off');
nexttile(3); legend('off');

set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% --- Save figure to output folder ---
output_folder = directory;
if ~exist(output_folder, 'dir'); mkdir(output_folder); end
exportgraphics(fig, fullfile(output_folder, 'VG_Curves.png'), 'Resolution', 300);

end
