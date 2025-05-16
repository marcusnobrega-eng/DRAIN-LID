%% =========================================================================
% 📂 File Location: Figures_And_Exports/Post_Processing.m
%
% 📊 POST-PROCESSING WORKFLOW — Mixed-Form Richards Model
% Description:
%   - Generates high-quality time series and profile figures
%   - Computes flow duration curves and wetting front movement
%   - Saves figures to organized folder with publication-ready style
%
% Author: Marcus Nóbrega, Ph.D.
% Updated: May 2025
%% =========================================================================

% === 1. Pre-checks and setup ============================================
fprintf('\n📊 Starting post-processing...\n');

% === 2. Extract variables ===============================================
z = params.z;            % Depth (m), Nz x 1
dz = params.dz;          % Cell thickness (m), Nz x 1
Nz = params.Nz;          % Number of nodes
Nt = size(head_out, 2);  % Number of time steps

% Time vector in seconds
time_vec = time_series(:)';   % [1 x Nt]
t_days   = time_vec / 86400;  % [days]
t_hours  = time_vec / 3600;   % [hours]

% === 3. Plot time series at key depths ==================================
top_idx = Nz - 1;
mid_idx = round(Nz / 2);
bot_idx = 1;

start_datetime = datetime(2023, 4, 6, 2, 10, 0);  % Example start date
use_datetime = true;  % Set to false to use elapsed time in hours

% === 8. State Animations ================================================
dt_days = 2/24;

animate_profiles(...
    theta_out, head_out, flux_out, ...
    time_series, ...
    params.z, ...
    fullfile(time_series_dir, 'Animations'), ...
    start_datetime, ...
    dt_days);  % ⬅️ Animate every 1 day


plot_time_series_nodes(...
    theta_out, head_out, flux_out, ...
    time_series, z, ...
    [bot_idx, mid_idx, top_idx], ...
    {'Bottom', 'Mid', 'Top'}, ...
    time_series_dir, ...
    Q_orifice_total, Q_spillway_total, ...
    start_datetime, use_datetime);

% === 4. Plot vertical profiles at selected times ========================
if length(t_hours) > 5
    selected_times = [0.05, 0.25, 0.5, 0.75, 0.95];  % Fractions of final time
    t_indices = round(selected_times * Nt);
    t_labels = strcat("t=", string(round(t_hours(t_indices), 2)), " hr");

plot_vertical_profiles(...
    theta_out, head_out, flux_out, ...
    z, t_indices, t_labels, ...
    profiles_dir);
end

% === 5. Compute and plot Flow Duration Curves ===========================
plot_flow_duration_curves(...
    flux_out, time_vec, z, ...
    fdc_dir);

% === 6. Wetting front depth analysis ====================================
compute_wetting_front_depth(...
    theta_out, time_vec, z, ...
    wetting_dir);

% === 7. θ-h retention curves at selected nodes ==========================
plot_theta_h_relationship(...
    theta_out, head_out, ...
    [bot_idx, mid_idx, top_idx], ...
    {'Bottom', 'Mid', 'Top'}, ...
    time_series_dir, t_days);


fprintf('✅ Post-processing completed. If you have any questions, please contact me at marcusnobrega.engcivil@gmail.com. Model developd by Dr. Marcus Nobrega Gomes Jr..\n');

%% Auxiliary Functions
function plot_time_series_nodes(theta, head, flux, time_vec, z, node_indices, node_labels, save_dir, Q_orifice_total, Q_spillway_total, start_datetime, use_datetime)
% =========================================================================
% 📈 PLOT TIME SERIES AT SELECTED NODES (Elapsed Time or Calendar Date)
% Inputs:
%   theta         – Water content matrix [Nz x Nt]
%   head          – Pressure head matrix [Nz x Nt]
%   flux          – Flux matrix [Nz+1 x Nt]
%   time_vec      – Time vector [1 x Nt] in seconds
%   z             – Depth vector [Nz x 1] in meters
%   node_indices  – Vector of selected node indices (e.g., [1, Nz/2, Nz])
%   node_labels   – Cell array with labels for each node
%   save_dir      – Directory to save the figure
%   Q_orifice_total – Flow vector [1 x Nt]
%   Q_spillway_total – Flow vector [1 x Nt]
%   start_datetime – datetime object of simulation start
%   use_datetime  – true for calendar date x-axis, false for elapsed time
% =========================================================================

% === Truncate data to valid range ===
n_saved = length(time_vec);
theta = theta(:,1:n_saved);
head = head(:,1:n_saved);
flux = flux(:,1:n_saved);
Q_orifice_total = Q_orifice_total(:,1:n_saved);
Q_spillway_total = Q_spillway_total(:,1:n_saved);

% === Colors and markers ===
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();
colors = [pallete.blue_colors(1,:); pallete.red_colors(1,:); [192,192,192]/255];
fs = 14; lw = 2;
markers = {'.','.','.'};

% === Time axis handling ===
if use_datetime
    x_vals = start_datetime + seconds(time_vec);
    x_label = 'Date';
else
    x_vals = time_vec / 3600;  % hours
    x_label = 'Time [hr]';
end

% === Start plotting ===
figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 10 7]);
tiledlayout(3,1,'TileSpacing','tight','Padding','compact');

% === θ(t)
nexttile; hold on; grid on; box on;
for i = 1:length(node_indices)
    plot(x_vals, theta(node_indices(i),:), ...
         'LineStyle', 'none', ...
         'Marker', markers{i}, ...
         'MarkerSize', 6, ...
         'Color', colors(i,:));
end
xlim([x_vals(1), x_vals(end)]);
ylabel('$\theta$ [m$^3$/m$^3$]', 'Interpreter', 'latex', 'FontSize', fs);
yyaxis right
    plot(x_vals, 1000*3600*flux(end,:), '-', 'Color', colors(1,:));
ylabel('$q~\mathrm{[mm \cdot h^{-1}]}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'YColor','black');
title('Soil Moisture at Selected Nodes', 'Interpreter', 'latex');
legend([node_labels, 'Top Flux'], 'Interpreter', 'latex', 'Location', 'best');
xlabel(x_label, 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'FontSize', fs, 'TickDir', 'out');

% === h(t)
nexttile; hold on; grid on; box on;
for i = 1:length(node_indices)
    plot(x_vals, head(node_indices(i),:), ...
          'LineStyle', 'none', ...
         'Marker', markers{i}, ...
         'MarkerSize', 6, ...
         'Color', colors(i,:));
end
xlim([x_vals(1), x_vals(end)]);
ylabel('$h$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
yyaxis right
    plot(x_vals, 1000*3600*flux(end,:), '-', 'Color', colors(1,:));
ylabel('$q~\mathrm{[mm \cdot h^{-1}]}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'YColor','black');
title('Pressure Head at Selected Nodes', 'Interpreter', 'latex');
legend([node_labels, 'Top Flux'], 'Interpreter', 'latex', 'Location', 'best');
xlabel(x_label, 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'FontSize', fs, 'TickDir', 'out');

% === q(t)
nexttile; hold on; grid on; box on;
for i = 1:length(node_indices)
    plot(x_vals, flux(node_indices(i)+1,:), ...
          'LineStyle', 'none', ...
         'Marker', markers{i}, ...
         'MarkerSize', 6, ...
         'Color', colors(i,:));

end
xlim([x_vals(1), x_vals(end)]);
plot(x_vals, Q_orifice_total, '-', 'LineWidth', 2, 'DisplayName', 'Orifice Flow');
plot(x_vals, Q_spillway_total, '--', 'LineWidth', 2, 'DisplayName', 'Spillway Flow');
ylabel('$q$ [m/s]', 'Interpreter', 'latex', 'FontSize', fs);
xlabel(x_label, 'Interpreter', 'latex', 'FontSize', fs);
title('Darcy Flux at Selected Nodes', 'Interpreter', 'latex');
legend([[node_labels, 'Top Flux'], 'Orifice','Spillway'], 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', fs, 'TickDir', 'out');

% === Save
exportgraphics(gcf, fullfile(save_dir, 'TimeSeries_SelectedNodes.png'), 'ContentType', 'image', 'Resolution', 300);
end

function plot_vertical_profiles(theta, head, flux, z, t_indices, t_labels, save_dir)
% =========================================================================
% 📉 PLOT VERTICAL PROFILES OF θ, h, q AT SELECTED TIMES
% Inputs:
%   theta      – Water content matrix [Nz x Nt]
%   head       – Pressure head matrix [Nz x Nt]
%   flux       – Flux matrix [Nz+1 x Nt]
%   z          – Depth vector [Nz x 1] (m)
%   t_indices  – Indices of selected time steps (e.g., [10 100 250])
%   t_labels   – Labels for each time (e.g., {'t=1h', 't=5h', ...})
%   save_dir   – Output directory
% =========================================================================

colors = lines(length(t_indices));
fs = 14; lw = 2;

figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 12 5.5]);
tiledlayout(1,3,"TileSpacing","tight","Padding","compact");

% === θ(z)
nexttile; hold on; grid on; box on;
for i = 1:length(t_indices)
    plot(theta(:,t_indices(i)), z, 'LineWidth', lw, 'Color', colors(i,:));
end
xlabel('$\theta$ [m$^3$/m$^3$]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
title('Water Content Profile', 'Interpreter', 'latex', 'FontSize', fs);
legend(t_labels, 'Interpreter','latex', 'Location','best');
set(gca, 'YDir', 'normal', 'FontSize', fs, 'TickDir', 'out');

% === h(z)
nexttile; hold on; grid on; box on;
for i = 1:length(t_indices)
    plot(head(:,t_indices(i)), z, 'LineWidth', lw, 'Color', colors(i,:));
end
xlabel('$h$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
title('Pressure Head Profile', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'YDir', 'normal', 'FontSize', fs, 'TickDir', 'out');

% === q(z)
nexttile; hold on; grid on; box on;
for i = 1:length(t_indices)
    plot(flux(1:end-1,t_indices(i)), z, 'LineWidth', lw, 'Color', colors(i,:));
end
xlabel('$q$ [m/s]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
title('Darcy Flux Profile', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'YDir', 'normal', 'FontSize', fs, 'TickDir', 'out');

% === Save figure
exportgraphics(gcf, fullfile(save_dir, 'VerticalProfiles_SelectedTimes.png'), 'ContentType', 'image', 'Resolution', 300);
end

function plot_flow_duration_curves(q, time_vec, z, save_dir)
% =========================================================================
% 📈 PLOT FLOW DURATION CURVES FOR POSITIVE AND NEGATIVE FLUXES
% Inputs:
%   q         – Flux matrix [Nz+1 x Nt], in m/s
%   time_vec  – Time vector [1 x Nt], in seconds
%   z         – Depth vector [Nz x 1], in meters
%   save_dir  – Output directory for figure
% =========================================================================

% Only internal fluxes (remove top and bottom BC fluxes if needed)
Nz = length(z);
depth_labels = {'Bottom', '1/4', 'Mid', '3/4', 'Top'};
nodes = [1, round(Nz/4), round(Nz/2), round(3*Nz/4), Nz];
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();
colors = [pallete.blue_colors(1,:); pallete.red_colors(1,:); pallete.green_colors(1,:); pallete.blue_colors(2,:) ; pallete.red_colors(2,:)];
% colors = lines(length(nodes));

fs = 14; lw = 2.1;

figure('Color','w','Units','inches','Position',[1 1 11 5.5]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% === POSITIVE FLUXES
nexttile; hold on; box on; grid on;
title('Positive Fluxes (Evaporation/Capilarity/Upward)', 'Interpreter', 'latex', 'FontSize', fs+1);
xlabel('Exceedance Probability [$\%$]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('Flux $q$ [$\mathrm{m/s}$]', 'Interpreter', 'latex', 'FontSize', fs);

plottedHandles = []; % Store handles for legend
plottedLabels = {};  % Store labels for legend

style = {'-','--','-.',':','-'};

for i = 1:length(nodes)
    q_node = q(nodes(i), :);
    pos_flux = q_node(q_node > 0);
    
    if ~isempty(pos_flux)
        [sorted_flux, ~] = sort(pos_flux, 'descend');
        p_exceed = (1:length(sorted_flux)) / length(sorted_flux) * 100;
        h = plot(p_exceed, sorted_flux, 'LineStyle', style{i}, 'LineWidth', lw, 'Color', colors(i,:));
        
        plottedHandles(end+1) = h; 
        plottedLabels{end+1} = depth_labels{i}; 
    end
end

legend(plottedHandles, plottedLabels, 'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontSize', fs, 'TickDir', 'out', 'YScale', 'log', 'XLim', [0 100]);

% === NEGATIVE FLUXES
nexttile; hold on; box on; grid on;
title('Negative Fluxes (Infiltration/Downward)', 'Interpreter', 'latex', 'FontSize', fs+1);
xlabel('Exceedance Probability [$\%$]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$|q|$ [$\mathrm{m/s}$]', 'Interpreter', 'latex', 'FontSize', fs);


plottedHandles = [];
plottedLabels = {};

for i = 1:length(nodes)
    q_node = q(nodes(i), :);
    neg_flux = q_node(q_node < 0);

    if ~isempty(neg_flux)
        [sorted_flux, ~] = sort(neg_flux, 'ascend');
        p_exceed = (1:length(sorted_flux)) / length(sorted_flux) * 100;
        h = plot(p_exceed, abs(sorted_flux), 'LineStyle', style{i}, 'LineWidth', lw, 'Color', colors(i,:));
        
        plottedHandles(end+1) = h; 
        plottedLabels{end+1} = depth_labels{i}; 
    end
end

legend(plottedHandles, plottedLabels, 'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'FontSize', fs, 'TickDir', 'out', 'YScale', 'log', 'XLim', [0 100]);


% === Export
exportgraphics(gcf, fullfile(save_dir, 'FlowDurationCurves.png'), 'Resolution', 300);
end

function compute_wetting_front_depth(theta, time_vec, z, save_dir)
% =========================================================================
% 🌊 COMPUTE AND PLOT WETTING FRONT MOVEMENT OVER TIME
%
% Inputs:
%   theta     – Water content matrix [Nz x Nt], [-]
%   time_vec  – Time vector [1 x Nt], in seconds
%   z         – Depth vector [Nz x 1], in meters (positive down)
%   save_dir  – Directory to save the output figure
%
% Method:
%   - Computes the wetting front depth as the first depth where
%     theta exceeds a dynamic threshold: initial theta + ε
% =========================================================================

n_saved = length(time_vec);
theta = theta(:,1:n_saved);

Nt = n_saved;
Nz = size(theta, 1);
theta_initial = theta(:, 1);
threshold = theta_initial + 0.05;  % 5% increase considered wetting

wetting_depth = nan(1, Nt);
for k = 1:Nt
    idx = find(theta(:,k) > threshold, 1, 'first');
    if ~isempty(idx)
        wetting_depth(k) = z(idx);  % in meters
    else
        wetting_depth(k) = NaN;
    end
end

% Plot wetting front
figure('Color','w','Units','inches','Position',[1 1 6.5 4]);
hold on; grid on; box on;
plot(time_vec / 3600, wetting_depth, '-o', 'LineWidth', 1.6, 'Color', [0.85 0.33 0.1]);
xlabel('Time [hr]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Wetting Front Depth [m]', 'Interpreter', 'latex', 'FontSize', 14);
title('Wetting Front Position Over Time', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 14, 'TickDir', 'out');

% Export
exportgraphics(gcf, fullfile(save_dir, 'WettingFrontMovement.png'), 'Resolution', 300);
end

function plot_theta_h_relationship(theta, head, node_indices, node_labels, save_dir, t_days)
% =========================================================================
% 📈 PLOT θ–h RETENTION CURVES FOR SELECTED NODES
%
% Inputs:
%   theta         – Water content matrix [Nz x Nt], [-]
%   head          – Pressure head matrix [Nz x Nt], [m]
%   node_indices  – Vector with selected node indices (e.g., [1, Nz/2, Nz])
%   node_labels   – Cell array of labels for each selected node
%   save_dir      – Directory to save the figure
%
% Output:
%   - Saves a plot of θ(h) curves for chosen depths
%
% Author: Marcus Nóbrega, Ph.D.
% =========================================================================

colors = lines(length(node_indices));

figure('Color','w','Units','inches','Position',[1 1 7 5]);
hold on; box on; grid on;

% Assume: 
% - head, theta: [Nz x Nt]
% - t_hours: time vector [1 x Nt] in hours
% - node_indices: vector of node indices
% - node_labels: cell array of labels
% - colors: use colormap to color by time

n_saved = length(t_days);
theta = theta(:,1:n_saved);
head = head(:,1:n_saved);

hold on;
cmap = cool(length(t_days));  % Or use 'parula', 'viridis', etc.
Mkt = ['*', 's', 'o'];
for i = 1:length(node_indices)
    idx = node_indices(i);
    
    % Plot scatter: h on x-axis, theta on y-axis, color by time
    scatter(head(idx,:), theta(idx,:), 40, t_days, ...
         Mkt(i), 'filled', 'DisplayName', node_labels{i});
end

% Add colorbar and label
cb = colorbar;
cb.Label.String = 'Time [days]';
cb.Label.Interpreter = 'latex';
cb.FontSize = 12;

colormap(turbo);  % Consistent with scatter color
xlabel('Pressure Head $h$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Water Content $\theta$ [$\mathrm{m^3/m^3}$]', 'Interpreter', 'latex', 'FontSize', 14);
title('$\theta$-$h$ Relationship at Selected Depths', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 14, 'TickDir', 'out', 'Box', 'on');


% Export
exportgraphics(gcf, fullfile(save_dir, 'Theta_H_Retention_Curves.png'), 'Resolution', 300);
end

function animate_profiles(theta, head, flux, time_vec, z, save_path, start_datetime, dt_days)
% =========================================================================
% 🎥 High-quality animation of hydrological vertical profiles
%      with time-evolving top flux hydrograph
%
% Inputs:
%   theta         – Water content matrix [Nz x Nt]
%   head          – Pressure head matrix [Nz x Nt]
%   flux          – Flux matrix [Nz+1 x Nt]
%   time_vec      – Time vector [1 x Nt], in seconds
%   z             – Depth vector [Nz x 1], in meters
%   save_path     – Output folder to save the MP4 animation
%   start_datetime– datetime object: simulation start
%   dt_days       – User-defined time step between frames (in days)
% =========================================================================

% === Setup
[Nz, Nt] = size(theta);
fs = 14; lw = 2;

% === Time conversion
t_datetime = start_datetime + seconds(time_vec);
flux_top = 1000 * 3600 * flux(end, :);  % Convert to mm/h

% === Frame selection: sample by exact dt_days
target_times = start_datetime + days(0:dt_days:days(t_datetime(end) - start_datetime));
frame_indices = arrayfun(@(t) find(t_datetime >= t, 1, 'first'), target_times);
frame_indices = unique([1, frame_indices(~isnan(frame_indices)), Nt]);
Nf = length(frame_indices) - 1;

% === Color palette
[~,~,~,~,~,pallete,~,~,~,~] = coloramps();
cT = pallete.blue_colors(2,:);    % Water content
cH = pallete.red_colors(1,:);     % Pressure head
cQ = pallete.green_colors(1,:);   % Flux

% === Figure and layout
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 13 8]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Hydrological Profile Animation', 'Interpreter', 'latex', 'FontSize', fs+2);

% === Top hydrograph (spans all top columns)
ax_flux_top = nexttile([1 3]); hold on; box on; grid on;
xlabel('Date', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$q_{\rm top}$ [mm/h]', 'Interpreter', 'latex', 'FontSize', fs);
set(ax_flux_top, 'FontSize', fs, 'TickDir', 'out');

top_flux_line = plot(ax_flux_top, t_datetime(1), flux_top(1), '-', ...
    'Color', [0 0 0.8], 'LineWidth', lw);
ax_flux_top.XLim = [t_datetime(1), t_datetime(end)];

% === Bottom panels
ax_theta = nexttile(4); hold on; box on; grid on;
xlabel('$\theta$ [$\mathrm{m^3/m^3}$]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
title('Water Content', 'Interpreter', 'latex', 'FontSize', fs);
set(ax_theta, 'YDir', 'normal', 'FontSize', fs, 'TickDir', 'out');
ylim(ax_theta, [min(z), max(z)]);

ax_head = nexttile(5); hold on; box on; grid on;
xlabel('$h$ [m]', 'Interpreter', 'latex', 'FontSize', fs);
title('Pressure Head', 'Interpreter', 'latex', 'FontSize', fs);
set(ax_head, 'YDir', 'normal', 'FontSize', fs, 'TickDir', 'out');
ylim(ax_head, [min(z), max(z)]);

ax_q = nexttile(6); hold on; box on; grid on;
xlabel('$q$ [m/s]', 'Interpreter', 'latex', 'FontSize', fs);
title('Darcy Flux', 'Interpreter', 'latex', 'FontSize', fs);
set(ax_q, 'YDir', 'normal', 'FontSize', fs, 'TickDir', 'out');
ylim(ax_q, [min(z), max(z)]);

% === Preallocate video frames
F(Nf) = struct('cdata',[],'colormap',[]);

% === Animation loop
for j = 1:Nf
    k = frame_indices(j);
    t_now = t_datetime(k);

    % --- Update hydrograph line
    set(top_flux_line, 'XData', t_datetime(1:k), 'YData', flux_top(1:k));
    if k > 1 && t_now > t_datetime(1)
        ax_flux_top.XLim = [t_datetime(1), t_now];
    end

    % --- Clear and update each profile plot
    cla(ax_theta); cla(ax_head); cla(ax_q);
    plot(ax_theta, theta(:,k), z, '-', 'LineWidth', lw, 'Color', cT);
    plot(ax_head,  head(:,k),  z, '-', 'LineWidth', lw, 'Color', cH);
    plot(ax_q,     flux(1:end-1,k), z, '-', 'LineWidth', lw, 'Color', cQ);

    % --- Global title with correct timestamp
    tl.Title.String = ['\textbf{Hydrological Profile at } ' ...
        datestr(t_now, 'dd-mmm-yyyy HH:MM')];

    drawnow;
    F(j) = getframe(fig);
end

% === Export video
if ~exist(save_path, 'dir')
    mkdir(save_path);
end
vidname = fullfile(save_path, 'HydrologicalProfilesAnimation.mp4');
v = VideoWriter(vidname, 'MPEG-4');
dt_video = 20; % Seconds
v.FrameRate = max(1, floor(dt_video / 60));  % ~1 min total
v.Quality = 100;

fprintf('🎞️ Exporting %d frames at %d fps → %.1f sec\n', ...
    Nf, v.FrameRate, Nf / v.FrameRate);

open(v); writeVideo(v, F); close(v);
fprintf('✅ Animation saved to: %s\n', vidname);
end
