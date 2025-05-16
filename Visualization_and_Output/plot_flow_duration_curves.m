function plot_flow_duration_curves(q_matrix, time_series, z, save_name)
% =========================================================================
% 📈 FLOW DURATION CURVES — Mixed-Form Richards Equation
% =========================================================================
% Description:
%   Plots flow duration curves (FDCs) at selected depths (top, mid, bottom)
%   for both positive (upward) and negative (downward) fluxes.
%   Uses adaptive time vector and SI units.
%
% Inputs:
%   q_matrix     - Matrix of Darcy fluxes [Nz+1 x Nt] in m/s
%   time_series  - Time stamps [1 x Nt] in seconds
%   z            - Depth vector [Nz x 1] in meters (positive upward)
%   save_name    - File base name for saving outputs (e.g., 'flow_fdc')
%
% Output:
%   Saves figure to: 'Figures_And_Exports/Diagnostics/<save_name>.png'
%
% Author: Marcus Nóbrega, Ph.D.
% Updated: May 2025
% =========================================================================

% === 0. Create output folder if it doesn't exist
out_dir = 'Figures_And_Exports/Diagnostics';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% === 1. Plotting settings
fs = 14; lw = 1.8;
colors = [0.00, 0.45, 0.70;   % Blue
          0.85, 0.33, 0.10;   % Red-orange
          0.00, 0.60, 0.50];  % Teal

% === 2. Nodes to analyze (interfaces: 1=bottom, Nz/2=mid, Nz+1=top)
Nz = size(z,1);
interfaces = [1, round((Nz+1)/2), Nz+1];
depth_labels = {'Bottom', 'Midpoint', 'Top'};

% === 3. Setup figure
figure('Color','w','Units','inches','Position',[1 1 11 5.5]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% === 4. Positive Flows (evaporation/upward)
nexttile; hold on; grid on; box on;
title('Upward (Positive) Flux Duration Curve', 'Interpreter','latex');
xlabel('Exceedance Probability [$\%$]', 'Interpreter','latex', 'FontSize', fs);
ylabel('$q$ [m/s]', 'Interpreter','latex', 'FontSize', fs);

for i = 1:length(interfaces)
    q_series = q_matrix(interfaces(i), :);
    q_pos = q_series(q_series > 0);
    q_sorted = sort(q_pos, 'descend');
    n_pos = numel(q_sorted);
    p_exceed = (1:n_pos) / n_pos * 100;

    plot(p_exceed, q_sorted, '-', 'Color', colors(i,:), ...
        'LineWidth', lw, 'DisplayName', depth_labels{i});
end

legend('Interpreter','latex', 'Location','southwest');
set(gca, 'FontSize', fs, 'TickDir','out', 'YScale','log');
xlim([0 100]);

% === 5. Negative Flows (infiltration/downward)
nexttile; hold on; grid on; box on;
title('Downward (Negative) Flux Duration Curve', 'Interpreter','latex');
xlabel('Exceedance Probability [$\%$]', 'Interpreter','latex', 'FontSize', fs);
ylabel('$|q|$ [m/s]', 'Interpreter','latex', 'FontSize', fs);

for i = 1:length(interfaces)
    q_series = q_matrix(interfaces(i), :);
    q_neg = q_series(q_series < 0);
    q_sorted = sort(abs(q_neg), 'descend');
    n_neg = numel(q_sorted);
    p_exceed = (1:n_neg) / n_neg * 100;

    plot(p_exceed, q_sorted, '-', 'Color', colors(i,:), ...
        'LineWidth', lw, 'DisplayName', depth_labels{i});
end

legend('Interpreter','latex', 'Location','southwest');
set(gca, 'FontSize', fs, 'TickDir','out', 'YScale','log');
xlim([0 100]);

% === 6. Export
filename = fullfile(out_dir, [save_name, '.png']);
print(gcf, filename, '-dpng', '-r400');

fprintf('✅ FDC plot saved to: %s\n', filename);
end
