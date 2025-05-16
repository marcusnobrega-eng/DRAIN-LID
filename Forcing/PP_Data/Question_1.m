clear; clc; close all;

%% Soil parameters from Table 1
soils = {'Sand','Loamy sand','Sandy loam','Silt loam','Loam','Clay loam','Sandy clay','Silty clay','Clay'};
theta_s = [0.395, 0.410, 0.435, 0.485, 0.451, 0.476, 0.426, 0.492, 0.482];
Ks = [1.76e-2, 1.56e-2, 3.47e-3, 7.20e-4, 6.95e-4, 2.45e-4, 2.17e-4, 1.03e-4, 1.28e-4];
hb = [-12.1, -9.0, -21.8, -78.6, -47.8, -63.0, -15.3, -49.0, -40.5];
b = [4.05, 4.38, 4.90, 5.30, 5.39, 8.52, 10.4, 10.4, 11.4];
hmf = [4.95, 6.13, 11.01, 16.68, 8.89, 20.88, 23.90, 29.22, 31.63];

%% Hm from theta 
colors = lines(length(soils));
figure;

for i = 1:length(soils)
    hb_i = hb(i);       % Air entry pressure (cm)
    b_i = b(i);         % Brooks-Corey exponent
    thetas_i = theta_s(i);
    
    Theta = linspace(0.01, 1, 100);  % Effective saturation (avoid Theta = 0)

    %theta = Theta * thetas_i;
    hm = hb_i * Theta.^(-b_i);  % Vectorized operation

    % Plot
    plot(Theta, hm, 'DisplayName', soils{i}, 'Color', colors(i,:)); hold on;
end

xlabel('\Theta [-]');
ylabel('h_m [cm]');
title('Inverse Brooks-Corey: h_m vs \Theta');
legend show;
set(gca, 'YDir', 'reverse')
set(gca, 'YScale', 'log');  % Log scale to capture full head range
grid on;

%% Question 1

% Create figure
figure('Name','Soil Water Retention & Conductivity','Position',[100 100 1200 500]);

% Plot theta vs hm (semi-log x)
subplot(1,2,1); hold on;
title('\theta vs h_m (Semi-log)');
ylabel('Matric head h_m [cm]'); xlabel('\theta [cm^3/cm^3]');
set(gca, 'YScale', 'log'); grid on;

% Plot K vs theta (log-log)
subplot(1,2,2); hold on;
title('K(\theta) vs \theta (Log-log)');
xlabel('\theta [cm^3/cm^3]'); ylabel('K [cm/s]');
set(gca, 'YScale', 'log', 'YScale', 'log'); grid on;

% Loop through soils
colors = lines(length(soils));
for i = 1:length(soils)
    % Matric head range
    hm = logspace(log10(abs(hb(i))), 3, 1000000); % [cm], log scale
    hm = -hm; % negative above water table
    
    % Compute theta
    theta = theta_s(i) * (abs(hm) ./ abs(hb(i))).^(-1/b(i));
    theta(hm >= hb(i)) = theta_s(i); % saturation
    
    % Compute K
    K = Ks(i) * (theta / theta_s(i)).^(2*b(i)+3);
    
    % Plot theta vs |hm|
    subplot(1, 2, 1);
    plot(theta, abs(hm), 'DisplayName', soils{i}, 'Color', colors(i,:), 'LineWidth', 1); hold on;
    xlabel('\theta [cm³/cm³]');
    ylabel('|h_m| [cm]');
    title('\theta vs |h_m|');
    set(gca, 'YScale', 'log');  % Log scale for pressure head
    legend show; grid on;

    % Plot K vs theta
    subplot(1, 2, 2);
    plot(theta, K, 'DisplayName', soils{i}, 'Color', colors(i,:), 'LineWidth', 1); hold on;
    xlabel('\theta [cm³/cm³]');
    ylabel('K(\theta) [cm/h]');
    title('Hydraulic Conductivity vs \theta');
    set(gca, 'YScale', 'log');  % ✅ Log scale for K
    legend show; grid on;
end

subplot(1,2,1); legend('show','Location','best');
box on;
subplot(1,2,2); legend('show','Location','best');
box on;

%% Question 2

% 2. Constants
depth_total = 100;  % total distance from water table to surface (cm)
n_soils = length(soils);

% 3. Initialize storage capacity vector
Storage_cm = zeros(n_soils,1);  % Storage in cm

% Water table
zwt = -100; % cm

% Initialize
n_soils = length(soils);
Storage_cm = zeros(n_soils,1);

for i = 1:n_soils
    thetaS = theta_s(i);
    hb_i = hb(i);
    b_i = b(i);

    % Find zb
    zb = zwt - hb_i;  % zb is negative!

    % Compute Storage
    term1 = abs(zb);
    term2 = (hb_i^b_i)/(1-b_i) * ( (abs(zb)+100)^(1-b_i) - (100)^(1-b_i) );
    S = thetaS * (term1 - term2);

    Storage_cm(i) = S;
end

% 5. Display results
fprintf('Storage Capacity (cm of water)\n');
for i = 1:n_soils
    fprintf('%-12s : %.2f cm\n', soils{i}, Storage_cm(i));
end

% 6. Optional: Plotting
figure('Color','w');
bar(Storage_cm)
set(gca,'XTickLabel',soils, 'XTickLabelRotation',45)
ylabel('Storage Capacity (cm)')
title('Storage Capacity until Saturation (Corrected)')
grid on;
%% Question 3

% Soil properties (choose contrasting ones from Table 1)
soils = {'Sand', 'Loam', 'Clay'};
Ks = [1.76e-2*3600, 6.95e-4*3600, 1.28e-4*3600];  % cm/h
theta_s = [0.395, 0.451, 0.482];    % cm³/cm³
hb = [-12.1, -47.8, -40.5];         % cm
hf = [4.95, 8.89, 31.63];           % cm

% soils = {'Sand','Loamy sand','Sandy loam','Silt loam','Loam','Clay loam','Sandy clay','Silty clay','Clay'};
% theta_s = [0.395, 0.410, 0.435, 0.485, 0.451, 0.476, 0.426, 0.492, 0.482];
% Ks = 3600 * [1.76e-2, 1.56e-2, 3.47e-3, 7.20e-4, 6.95e-4, 2.45e-4, 2.17e-4, 1.03e-4, 1.28e-4];
% hb = [-12.1, -9.0, -21.8, -78.6, -47.8, -63.0, -15.3, -49.0, -40.5];
% b = [4.05, 4.38, 4.90, 5.30, 5.39, 8.52, 10.4, 10.4, 11.4];
% hf = [4.95, 6.13, 11.01, 16.68, 8.89, 20.88, 23.90, 29.22, 31.63];

% Time domain
t = linspace(0, 20/3600, 20);  % hours
dt = diff(t); dt = [dt dt(end)];

% Green-Ampt cumulative infiltration function F(t)
% compute_ft = @(Ks, hf, theta_s, theta_i, t) arrayfun(@(ti) ...
%      fsolve(@(F) ti - (F/Ks - hf * (theta_s - theta_i)/Ks * log(1 + F / (hf * (theta_s - theta_i)))), ...
%      1, optimset('Display','off')), t);

compute_ft = @(Ks, hf, theta_s, theta_i, t) arrayfun(@(ti) ...
     fsolve(@(F) (F - (Ks*ti + hf * (theta_s - theta_i) * log(1 + F / (hf * (theta_s - theta_i))))), ...
     1, optimset('Display','off')), t);

figure;
for i = 1:length(soils)
    % Initial moisture
    theta_i = 0;
    
    % Compute cumulative infiltration
    F = compute_ft(Ks(i), hf(i), theta_s(i), theta_i, t);
    
    % Compute infiltration rate
    ft = Ks(i) * (1 + (hf(i) * (theta_s(i) - theta_i)) ./ F);

    plot(t, ft / 3600, 'LineWidth', 2, 'DisplayName', soils{i}); hold on;
    % set(gca, 'YScale', 'log');
end


xlabel('Time [h]');
ylabel('Infiltration Capacity f(t) [cm/sec]');
title('Green-Ampt Infiltration Capacity');
legend show; grid on;

% % Explore effect of initial moisture on one soil (Loam)
theta_i_vals = [0.05, 0.15, 0.25];
figure;
for j = 1:length(theta_i_vals)
    theta_i = theta_i_vals(j);
    i = 2; % Loam
    F = compute_ft(Ks(i), hf(i), theta_s(i), theta_i, t);
    ft = Ks(i) * (1 + (hf(i) * (theta_s(i) - theta_i)) ./ F);
    plot(t, ft / 3600, 'LineWidth', 2, 'DisplayName', ['\theta_i = ', num2str(theta_i)]); hold on;
end

xlabel('Time [h]');
ylabel('Infiltration Capacity f(t) [cm/sec]');
title('Effect of Initial Moisture on Loam Infiltration');
legend show; grid on;

