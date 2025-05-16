%% ---------------------- coloramps.m ----------------------
% Returns various color ramps and color palettes for plotting hydrodynamic data.
%
% Outputs:
% Spectrum       - Full-spectrum color gradient
% depth_ramp     - Smooth water depth colormap
% terrain_ramp   - Terrain elevation colormap
% blue_ramp      - Gamma-corrected blue gradient
% blues_2        - Discrete to smooth blue shades
% pallete        - Struct with blue, red, green color sets
% Depth_RAS      - Maria-style blue gradient for depth
% Terrain_RAS    - Terrain style colormap
% Velocity_RAS   - Velocity-based heatmap
% WSE_RAS        - Water Surface Elevation gradient

function [Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS,WSE_RAS] = coloramps()

%% ========== Full Spectrum Ramp ==========
RGB = [ % Custom full-spectrum colors
    0.1127 0 0.3515
    0.2350 0 0.6663
    0.3536 0 1.0000
    0.4255 0 1.0000
    0.4384 0 1.0000
    0.3888 0 1.0000
    0.2074 0 1.0000
    0      0 1.0000
    0 0.4124 1.0000
    0 0.6210 1.0000
    0 0.7573 0.8921
    0 0.8591 0.6681
    0 0.9642 0.4526
    0 1.0000 0.1603
    0 1.0000 0
    0 1.0000 0
    0 1.0000 0
    0 1.0000 0
    0.4673 1.0000 0
    0.8341 1.0000 0
    1.0000 0.9913 0
    1.0000 0.8680 0
    1.0000 0.7239 0
    1.0000 0.5506 0
    1.0000 0.3346 0
    1.0000 0 0
    1.0000 0 0
    1.0000 0 0
    1.0000 0 0
    0.9033 0 0
    0.7412 0 0
    0.5902 0 0];
Spectrum = interp1(linspace(1, 255, size(RGB,1)), RGB, 1:255);

%% ========== Water Depth Colormap ==========
depth_colors = [
    0.27 0.47 0.68
    0.41 0.68 0.84
    0.95 0.95 0.72
    0.98 0.74 0.43
    0.85 0.33 0.10
    1.00 0.00 0.00
    0.50 0.00 0.25
    0.50 0.00 0.50];
depth_values = linspace(0.0, 4.0, size(depth_colors,1));
depth_ramp = interp1(depth_values, depth_colors, linspace(min(depth_values), max(depth_values), 255), 'pchip');

%% ========== Terrain Elevation Colormap ==========
terrain_colors = [
    [0, 100, 0]
    [152, 251, 152]
    [255, 255, 0]
    [184, 134, 11]
    [139, 69, 19]
    [101, 67, 33]
    [255, 0, 0]
    [139, 0, 0]
]/255;

numColors = 256;
terrain_ramp = [];
for i = 1:length(terrain_colors)-1
    ramp = [linspace(terrain_colors(i,1), terrain_colors(i+1,1), numColors/8)', ...
            linspace(terrain_colors(i,2), terrain_colors(i+1,2), numColors/8)', ...
            linspace(terrain_colors(i,3), terrain_colors(i+1,3), numColors/8)'];
    terrain_ramp = [terrain_ramp; ramp];
end
terrain_ramp = smoothdata(terrain_ramp, 'gaussian', 10);

%% ========== Blue Gradient (Gamma-Corrected) ==========
num_colors = 100;
start_color = [0, 0.2, 0.6];
end_color = [0.8, 0.9, 1];
gamma = 2;

r = linspace(start_color(1), end_color(1), num_colors).^gamma;
g = linspace(start_color(2), end_color(2), num_colors).^gamma;
b = linspace(start_color(3), end_color(3), num_colors).^gamma;
blue_ramp = [r' g' b'];

%% ========== Blues 2 ==========
custom_blues = [
    0.031, 0.365, 0.639
    0.098, 0.569, 0.737
    0.306, 0.675, 0.816
    0.529, 0.765, 0.882
    0.725, 0.851, 0.941
    0.882, 0.937, 0.988
    0.973, 0.988, 1.000];
blues_2 = interp1(linspace(0, 1, size(custom_blues,1)), custom_blues, linspace(0, 1, 256));

%% ========== Color Palettes ==========
% Blue shades
pallete.blue_colors = [
    31, 102, 169
    52, 148, 204
    141, 197, 228
    190, 223, 240
    219, 236, 245
]/255;

% Red shades
pallete.red_colors = [
    159, 0, 0
    196, 103, 102
    216, 165, 166
    240, 210, 210
]/255;

% Green shades
pallete.green_colors = [
    31, 111, 112
    84, 162, 161
    159, 200, 200
    200, 225, 225
]/255;

%% ========== Depth RAS ("Maria" style palette) ==========
temp = [... % 32 shades from cyan to deep blue
    0,255,255; 0,247,251; 0,239,248; 0,230,244; 0,222,240;
    0,214,236; 0,206,233; 0,197,229; 0,189,225; 0,181,221;
    0,173,218; 0,165,214; 0,156,210; 0,148,206; 0,140,203;
    0,132,199; 0,123,195; 0,115,191; 0,107,188; 0,99,184;
    0,90,180; 0,82,176; 0,74,173; 0,66,169; 0,58,165;
    0,49,161; 0,41,158; 0,33,154; 0,25,150; 0,16,146;
    0,8,143; 0,0,139]/255;
Depth_RAS = interp1(linspace(1, 255, size(temp,1)), temp, 1:255);

%% ========== Terrain RAS ==========
temp = [... % From green to gray to white
    102,205,170; 136,216,132; 171,228,93; 206,239,55; 241,250,16;
    222,239,0; 165,210,0; 107,181,0; 49,153,0; 8,129,0;
    66,138,0; 123,146,0; 181,154,0; 239,163,0; 237,139,0;
    210,101,0; 184,64,0; 157,26,0; 141,3,3; 147,12,12;
    152,22,22; 158,31,31; 164,41,41; 158,59,59; 150,78,78;
    141,98,98; 133,117,117; 140,140,140; 169,167,167;
    198,195,195; 227,223,223; 255,250,250]/255;
Terrain_RAS = interp1(linspace(1, 255, size(temp,1)), temp, 1:255);

%% ========== Velocity RAS ==========
temp = [... % Blue to red velocity gradient
    0,0,139; 0,0,167; 0,0,195; 0,0,223; 0,0,252; 30,50,232;
    65,107,205; 100,165,178; 135,223,151; 163,241,119; 191,245,84;
    217,249,49; 244,253,14; 255,242,0; 255,220,0; 255,198,0;
    255,177,0; 255,154,0; 255,131,0; 255,108,0; 255,84,0;
    251,67,0; 240,60,0; 229,53,0; 218,47,0; 206,40,0;
    195,33,0; 184,27,0; 173,20,0; 162,13,0; 150,7,0;
    139,0,0]/255;
Velocity_RAS = interp1(linspace(1, 255, size(temp,1)), temp, 1:255);

%% ========== Water Surface Elevation (WSE) RAS ==========
temp = [... % Green to red to purple
    127,255,0; 103,231,0; 78,206,0; 53,181,0; 28,156,0;
    4,132,0; 41,148,0; 91,173,0; 140,198,0; 189,222,0;
    239,247,0; 255,243,0; 255,226,0; 255,209,0; 255,191,0;
    255,174,0; 255,149,0; 255,117,0; 255,85,0; 255,53,0;
    255,21,0; 255,12,13; 255,50,53; 255,87,92; 255,124,131;
    255,161,171; 255,186,205; 255,149,215; 255,111,225;
    255,74,235; 255,37,245; 255,0,255]/255;
WSE_RAS = interp1(linspace(1, 255, size(temp,1)), temp, 1:255);

end
