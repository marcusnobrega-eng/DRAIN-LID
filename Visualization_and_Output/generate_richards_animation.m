function generate_richards_animation(theta_out, head_out, flux_out, params, results_dir)
% =========================================================================
% 📽️ ANIMATION GENERATOR — Soil Moisture, Pressure Head, and Flux
% Description: Creates three animation files (MP4) showing how 
%              theta, head, and flux evolve over time across the vertical profile.
%
% Inputs:
%   theta_out   - [Nz x Nt] matrix of volumetric water content [–]
%   head_out    - [Nz x Nt] matrix of pressure head [m]
%   flux_out    - [Nz+1 x Nt] matrix of flux [m/s]
%   params      - Struct with fields:
%                 * z (depths, m)
%                 * save_times (1 x Nt) in seconds
%   results_dir - Path to base results folder (e.g., 'Results_MySim')
%
% Author: Marcus Nóbrega, Ph.D.
% Updated: May 2025
% =========================================================================

%% === 1. Setup and Paths ===
output_dir = fullfile(results_dir, 'Figures', 'Animations');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

z     = params.z(:);      % [m]
times = params.save_times / 3600;  % [hr]
Nz    = length(z);
Nt    = length(times);

colors = lines(1);
fs = 14;
lw = 2;

%% === 2. Animation — Soil Moisture (theta) ===
theta_video = VideoWriter(fullfile(output_dir, 'Animation_Theta.mp4'), 'MPEG-4');
theta_video.Quality = 100; theta_video.FrameRate = 12;
open(theta_video);

figure('Color','w','Position',[100 100 500 600]);
for t = 1:Nt
    clf;
    plot(theta_out(:,t), z, 'Color', colors(1,:), 'LineWidth', lw, 'Marker', 'o');
    xlabel('$\theta~[-]$', 'Interpreter','latex','FontSize',fs);
    ylabel('$z~[\mathrm{m}]$', 'Interpreter','latex','FontSize',fs);
    title(sprintf('Volumetric Water Content\nTime = %.2f hr', times(t)), 'Interpreter','latex','FontSize',fs);
    set(gca,'FontSize',fs,'TickDir','out','YDir','normal');
    xlim([min(theta_out(:)), max(theta_out(:))*1.05]);
    ylim([min(z), max(z)]);
    grid on; drawnow;
    writeVideo(theta_video, getframe(gcf));
end
close(theta_video);

%% === 3. Animation — Pressure Head (h) ===
head_video = VideoWriter(fullfile(output_dir, 'Animation_Head.mp4'), 'MPEG-4');
head_video.Quality = 100; head_video.FrameRate = 12;
open(head_video);

figure('Color','w','Position',[100 100 500 600]);
for t = 1:Nt
    clf;
    plot(head_out(:,t), z, 'Color', colors(1,:), 'LineWidth', lw, 'Marker', 's');
    xlabel('$h~[\mathrm{m}]$', 'Interpreter','latex','FontSize',fs);
    ylabel('$z~[\mathrm{m}]$', 'Interpreter','latex','FontSize',fs);
    title(sprintf('Pressure Head\nTime = %.2f hr', times(t)), 'Interpreter','latex','FontSize',fs);
    set(gca,'FontSize',fs,'TickDir','out','YDir','normal');
    xlim([min(head_out(:)), max(head_out(:))*1.05]);
    ylim([min(z), max(z)]);
    grid on; drawnow;
    writeVideo(head_video, getframe(gcf));
end
close(head_video);

%% === 4. Animation — Vertical Flux (q) ===
flux_video = VideoWriter(fullfile(output_dir, 'Animation_Flux.mp4'), 'MPEG-4');
flux_video.Quality = 100; flux_video.FrameRate = 12;
open(flux_video);

figure('Color','w','Position',[100 100 500 600]);
for t = 1:Nt
    clf;
    plot(flux_out(1:end-1,t), z, 'Color', colors(1,:), 'LineWidth', lw, 'Marker', '^');
    xlabel('$q~[\mathrm{m/s}]$', 'Interpreter','latex','FontSize',fs);
    ylabel('$z~[\mathrm{m}]$', 'Interpreter','latex','FontSize',fs);
    title(sprintf('Vertical Flux\nTime = %.2f hr', times(t)), 'Interpreter','latex','FontSize',fs);
    set(gca,'FontSize',fs,'TickDir','out','YDir','normal');
    xlim([min(flux_out(:)), max(flux_out(:))*1.05]);
    ylim([min(z), max(z)]);
    grid on; drawnow;
    writeVideo(flux_video, getframe(gcf));
end
close(flux_video);

end