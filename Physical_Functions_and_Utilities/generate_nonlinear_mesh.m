%% =========================================================================
%  Function : generate_nonlinear_mesh
%  Purpose  : Generates a 1D vertical mesh refined toward the surface
%             using a nonlinear stretching method (Hydrus-style 🌱).
%
%  Author   : Marcus Nóbrega, Ph.D.
%  Updated  : May 2025
%
%  Inputs:
%    Nz            – Number of vertical nodes (e.g., 40)
%    L             – Total profile depth [m] (positive number, e.g., 1.0 m)
%    nonlin_factor – Stretching factor (>1 increases resolution near top)
%
%  Outputs:
%    z  – Node depths from surface to bottom (vector, [m], size = Nz)
%         - Starts at 0 (soil surface), ends at -L (soil bottom)
%    dz – Cell thicknesses [m] (same size as z)
%
%  Notes:
%    - z is negative, representing depth below ground.
%    - The topmost layers become thinner as `nonlin_factor` increases.
% =========================================================================
function [z, dz] = generate_nonlinear_mesh(Nz, L, nonlin_factor, mesh_dir)

    % 🧭 Uniformly spaced pseudo-depth coordinate from 0 (top) to 1 (bottom)
    s = linspace(0, 1, Nz)';

    % Generate Regular Mesh
    z_uniform = linspace(-L, 0, Nz)';

    % 🎯 Apply nonlinear stretching (exponentially refine near surface)
    s_refined = 1 - (1 - s).^nonlin_factor;

    % 📏 Convert to actual depth values from 0 (top) to -L (bottom)
    z = -L + L * s_refined;  % z(1) = -L, z(end) = 0

    % 📐 Compute cell thicknesses (dz_i = z(i+1) - z(i))
    dz = diff(z)';
    dz = [dz, dz(end)];  % Last cell thickness = same as previous
    dz_uniform = diff(z_uniform)';
    dz_uniform = [dz_uniform, dz_uniform(end)];

    % Plot
    figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5.7, 4]);  % ~Half A4 size
    plot(dz_uniform, z_uniform, 'o-', 'LineWidth', 2, 'DisplayName', 'Regular Mesh'); hold on;
    plot(dz, z, 's--', 'LineWidth', 2, 'DisplayName', sprintf('Refined Mesh (n=%.1f)', nonlin_factor));
    
    ylabel('Depth from Surface $z$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Cell Thickness $\Delta z$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
    title('Comparison of Regular and Refined Vertical Meshes', 'Interpreter', 'latex', 'FontSize', 15);
    legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 13, 'TickDir', 'out', 'LineWidth', 1.5);

    % Export as PNG (300 DPI)
    
    print(gcf, fullfile(mesh_dir, 'Spatial_Mesh.png'), '-dpng', '-r300');
end
