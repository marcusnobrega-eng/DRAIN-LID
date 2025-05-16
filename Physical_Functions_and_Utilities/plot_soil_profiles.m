%% =========================================================================
% 📂 File Location : Visualization_And_Output/plot_soil_profiles.m
%
% 📊 FUNCTION: plot_soil_profiles
% Purpose    : Plots vertical profiles of h, H, θ, and q with high-quality
%              visual style for scientific publication.
%
% Inputs:
%   h             – Current pressure head [Nz x 1]
%   h_old         – Previous step head     [Nz x 1]
%   params        – Struct with z, dz, etc.
%   t             – Current time [s]
%   tstep         – Time step index
%   theta_fun     – Function handle for θ(h)
%   K_fun         – Function handle for K(h)
%   residual_fun  – Function handle returning [F, q] from h
%
% Author     : Marcus Nóbrega, Ph.D.
% Updated    : May 2025
%% =========================================================================

function plot_soil_profiles(h, h_old, params, t, tstep, theta_fun, K_fun, residual_fun)

    %% === 1. Evaluate State Variables =====================================
    theta_now = theta_fun(h, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
    K_now     = K_fun(h, params.Ks, params.theta_r, params.theta_s, params.alpha, params.n, params.m);
    total_head = h + params.z(:)';

    top_bc = 0; bottom_bc = 0; source = zeros(size(h));
    [~, q_now] = residual_fun(h, h_old, theta_now, theta_now, K_now, params, top_bc, bottom_bc, source);

    % Define Nature-style color palette
    nature_colors = [
        0.0, 0.45, 0.70;   % Blue
        0.85, 0.33, 0.10;  % Red-orange
        0.00, 0.60, 0.50;  % Teal
        0.94, 0.89, 0.26;  % Yellow
    ];

    %% === 2. Set Font, Markers, Layout ====================================
    fontSize = 14;
    lineWidth = 2.0;
    markerSize = 6;

    %% === 3. Pressure Head ================================================
    nexttile; hold on; grid on;
    plot(h, params.z, '-', 'Color', nature_colors(1,:), 'LineWidth', lineWidth, 'Marker','*','MarkerSize',markerSize);
    plot(h_old, params.z, '--', 'Color', nature_colors(2,:), 'LineWidth', lineWidth);
    xlabel('$h$ [m]', 'Interpreter', 'latex');
    ylabel('$z$ [m]', 'Interpreter', 'latex');
    title('Pressure Head', 'Interpreter', 'latex');
    legend({'Current','Previous'}, 'Interpreter','latex', 'Location','best');
    set(gca, 'FontSize', fontSize, 'TickDir','out', 'Box','on');

    %% === 4. Total Head ===================================================
    nexttile; hold on; grid on;
    plot(total_head, params.z, '-', 'Color', nature_colors(3,:), 'LineWidth', lineWidth, 'Marker','*','MarkerSize',markerSize);
    xlabel('$h + z$ [m]', 'Interpreter', 'latex');
    ylabel('$z$ [m]', 'Interpreter', 'latex');
    title('Total Head', 'Interpreter', 'latex');
    set(gca, 'FontSize', fontSize, 'TickDir','out', 'Box','on');

    %% === 5. Water Content ================================================
    nexttile; hold on; grid on;
    plot(theta_now, params.z, '-', 'Color', nature_colors(4,:), 'LineWidth', lineWidth, 'Marker','*','MarkerSize',markerSize);
    xlabel('$\theta$ [--]', 'Interpreter', 'latex');
    ylabel('$z$ [m]', 'Interpreter', 'latex');
    title('Water Content', 'Interpreter', 'latex');
    set(gca, 'FontSize', fontSize, 'TickDir','out', 'Box','on');

    %% === 6. Vertical Flux ================================================
    nexttile; hold on; grid on;
    plot(q_now(1:end-1), params.z, '-', 'Color', nature_colors(2,:), 'LineWidth', lineWidth, 'Marker','*','MarkerSize',markerSize);
    xlabel('$q$ [m/s]', 'Interpreter', 'latex');
    ylabel('$z$ [m]', 'Interpreter', 'latex');
    title('Darcy Flux', 'Interpreter', 'latex');
    set(gca, 'FontSize', fontSize, 'TickDir','out', 'Box','on');

    %% === 7. Super Title and Export =======================================
    sgtitle(sprintf('Time = %.2f min  |  Step = %d', t / 60, tstep), ...
            'Interpreter','latex', 'FontSize', fontSize + 2);

    % Export as high-resolution figure
    % print(gcf, sprintf('Visualization_and_Output/VerticalProfiles_T%04d.png', tstep), '-dpng', '-r400');
end 