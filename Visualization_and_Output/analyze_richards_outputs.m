function analyze_richards_outputs(head_out, theta_out, flux_out, params, results_dir)
% =========================================================================
% 📈 Post-processing of Richards Model Outputs (Time Series + Profiles)
%
% Inputs:
%   head_out   - [Nz x Nt] matrix of pressure head [m]
%   theta_out  - [Nz x Nt] matrix of volumetric water content [–]
%   flux_out   - [Nz+1 x Nt] matrix of Darcy fluxes [m/s]
%   params     - Struct with z [Nz x 1], dz [Nz x 1], save_times [1 x Nt]
%   results_dir- Path to main results directory (e.g., 'Results_MySim')
%
% Author: Marcus Nóbrega, Ph.D.
% Updated: May 2025
% =========================================================================

% === 1. Setup ============================================================
t_hr = params.save_times / 3600;    % [s] → [hr]
z     = params.z(:);                % [m]
dz    = params.dz(:);               % [m]
Nt    = length(t_hr);
Nz    = length(z);
z_vec = z;                          % For y-axis labeling

% Create output folders
ts_dir = fullfile(results_dir, 'Figures', 'TimeSeries');
pr_dir = fullfile(results_dir, 'Figures', 'Profiles');
if ~exist(ts_dir, 'dir'), mkdir(ts_dir); end
if ~exist(pr_dir, 'dir'), mkdir(pr_dir); end

% === 2. Nodes of Interest ===============================================
top    = Nz;
middle = round(Nz / 2);
bottom = 1;

% === 3. Time Series Plots ===============================================
figure('Color','w','Units','inches','Position',[1 1 9 6]);
tiledlayout(3,1,"TileSpacing","compact","Padding","compact");

% θ(t)
nexttile; hold on; grid on; box on;
plot(t_hr, theta_out(bottom,:), 'LineWidth', 2);
plot(t_hr, theta_out(middle,:), 'LineWidth', 2);
plot(t_hr, theta_out(top,:), 'LineWidth', 2);
ylabel('$\theta~[-]$', 'Interpreter','latex','FontSize',14);
title('Soil Moisture at Selected Nodes','Interpreter','latex');
legend('Bottom','Middle','Top','Location','best','Interpreter','latex');
set(gca,'TickDir','out','FontSize',14);

% h(t)
nexttile; hold on; grid on; box on;
plot(t_hr, head_out(bottom,:), 'LineWidth', 2);
plot(t_hr, head_out(middle,:), 'LineWidth', 2);
plot(t_hr, head_out(top,:), 'LineWidth', 2);
ylabel('$h~[\mathrm{m}]$', 'Interpreter','latex','FontSize',14);
title('Pressure Head at Selected Nodes','Interpreter','latex');
set(gca,'TickDir','out','FontSize',14);

% q(t)
nexttile; hold on; grid on; box on;
plot(t_hr, flux_out(bottom,:),'LineWidth', 2);
plot(t_hr, flux_out(middle,:),'LineWidth', 2);
plot(t_hr, flux_out(top+1,:),'LineWidth', 2);
ylabel('$q~[\mathrm{m/s}]$', 'Interpreter','latex','FontSize',14);
xlabel('Time [hr]', 'Interpreter','latex','FontSize',14);
title('Vertical Flux at Selected Nodes','Interpreter','latex');
set(gca,'TickDir','out','FontSize',14);

% Export
exportgraphics(gcf, fullfile(ts_dir, 'TimeSeries_Theta_Head_Flux.png'), 'Resolution', 300);

% === 4. Depth Profiles at Key Times =====================================
key_idx = round([0.1 0.5 0.9] * Nt);
labels = arrayfun(@(i) sprintf('$t = %.1f$ hr', t_hr(i)), key_idx, 'UniformOutput', false);
colors = lines(length(key_idx));

figure('Color','w','Units','inches','Position',[1 1 10 4]);
tiledlayout(1,3,"TileSpacing","compact","Padding","compact");

% θ(z)
nexttile; hold on; grid on; box on;
for i = 1:length(key_idx)
    plot(theta_out(:,key_idx(i)), z_vec, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('$\theta~[-]$', 'Interpreter','latex','FontSize',14);
ylabel('$z~[\mathrm{m}]$', 'Interpreter','latex','FontSize',14);
title('Water Content Profiles','Interpreter','latex');
legend(labels, 'Interpreter','latex','Location','best');
set(gca,'TickDir','out','FontSize',14);

% h(z)
nexttile; hold on; grid on; box on;
for i = 1:length(key_idx)
    plot(head_out(:,key_idx(i)), z_vec, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('$h~[\mathrm{m}]$', 'Interpreter','latex','FontSize',14);
title('Pressure Head Profiles','Interpreter','latex');
set(gca,'TickDir','out','FontSize',14);

% q(z)
nexttile; hold on; grid on; box on;
for i = 1:length(key_idx)
    plot(flux_out(1:end-1,key_idx(i)), z_vec, '-', 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('$q~[\mathrm{m/s}]$', 'Interpreter','latex','FontSize',14);
title('Flux Profiles','Interpreter','latex');
set(gca,'TickDir','out','FontSize',14);

% Export
exportgraphics(gcf, fullfile(pr_dir, 'Depth_Profiles_KeyTimes.png'), 'Resolution', 300);

end
