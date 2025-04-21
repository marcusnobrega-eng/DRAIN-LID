% Richards Model - Finite-Difference Explicit Scheme
% Developer: Marcus Nobrega
% Goal - Solve 1-D Richards with Known Boundary Conditions
% Next Steps: Implement a hydrologic model and include boundary conditions
% to allow ponding depth
% Moreover, implement pollutant transport model

% --- Code working only for negative head boundary conditions ---- %

clear all
clc
close all
%% 1.0 - Van Genutchen Model - Soil Parameters

% Van Genuthen Paper Example
% a = 0.005;
% theta_s = 0.5;
% theta_r = 0.1;
% n = 2;
% m = 1 - 1/n;
% Ks = 100/24/60/60; % cm/s

% Web Example
% a = 0.0335; % 1/cm
% theta_s = 0.368;
% theta_r = 0.102;
% n = 2;
% % m = 4;
% m = 1 - 1/n;
% Ks = 0.00922; % cm/s

% Sand Example
% a = 0.087424; % 1/cm
% theta_s = 0.370687;
% theta_r = 0.043019;
% n = 1.57535;
% % m = 0.5;
% m = 1 - 1/n;
% Ks = 50; % mm/h
% Ks = Ks/10/3600; % cm/s

% HYDRUS1D Example
a = 0.036; % 1/cm
theta_s = 0.43;
theta_r = 0.078;
n = 1.56;
% m = 0.5;
m = 1 - 1/n;
Ks = 24.96; % cm/day
Ks = Ks/24/60/60; % cm/s

% Dussailant Paper Example
% a = 0.0355; % 1/cm
% theta_s = 0.368;
% theta_r = 0.102;
% n = 2;
% % m = 0.5;
% m = 1 - 1/n;
% Ks = 0.00922; % cm/s


% ------- Functions - In terms of Pressure --------- %
theta_function = @(h) (abs(sign((sign(h) - 1))).*(theta_r + (theta_s - theta_r)./((1 + (a.*abs(h)).^n).^m)) + (sign(sign(h) + 1).*theta_s)); % Soil Moisture for h < 0
theta_matlab = sym(theta_function); % Creating a Symbolic Function Soil Moisture
syms h
Theta_matlab = abs(sign((sign(h) - 1))).*(theta_matlab - theta_r)./(theta_s - theta_r) + (sign(sign(h) + 1).*1); % Effective Soil Moisture
Kr = Theta_matlab.^(1/2).*(1 - (1 - Theta_matlab.^(1/m)).^m).^2; % Kr = K / Ksat; - Relative K
K_function = Kr.*Ks; % Symbolic Function of Unsaturated Hydraulic Conductivity in terms of Ksat
K_function = matlabFunction(K_function); % Matlab Function of K in terms of h

% Other Way
K_function_head = @(h) (abs(sign((sign(h) - 1))).*((1 - ((a.abs(h)).^(n-2)).*(1 + (a.*abs(h)).^n).^(-m)))./((1 + (a.*abs(h)).^n).^(2*m)) + (sign(sign(h) + 1).*Ks));

% Inverse Functions
theta_matlab = sym(theta_function); % Creating a Symbolic Function
% This head in terms of theta 
% h_theta_function = finverse(theta_matlab); % Calculating the Inverse. Now h is theta. h_theta_function(theta)
% h_theta_function = matlabFunction(h_theta_function); % Calculating the Inverse. Now h is theta. h_theta_function(theta)

% Solving Analytically (by hand)
h_theta_function_2 = @(theta) ((-1)*(1/a*(-1 + 1./((theta - theta_r)./(theta_s - theta_r)).^(1/m)).^(1/n)));

% Van Genutchen Plots
p_plot = linspace(-1e0,-1e6,1e6); % Range of plotted pressure
set(gcf,'units','inches','position',[2,2,8,3])
ax1 = subplot(1,3,1);
plot(theta_function(p_plot),p_plot,'linewidth',2,'Color','red');
hold on
xlabel('$\theta$','Interpreter','latex')
ylabel('Head (cm)','Interpreter','latex')
grid on
xlim([0 0.5])
ax1.TickLength = [0.025 0.015];
set(gca, 'YScale', 'log')
ax1.XTick = [0 0.1 0.2 0.3 0.4 0.5];
ax1.TickLength = [0.03 0.015];
ax2 = subplot(1,3,2);
plot(K_function(p_plot),p_plot,'linewidth',2,'Color','black');
xlabel('$K$ (cm/s)','Interpreter','latex')
ylabel('Head (cm)','Interpreter','latex')
grid on
set(gca, 'YScale', 'log')

ax3 = subplot(1,3,3);
plot(theta_function(p_plot),K_function(p_plot),'linewidth',2,'Color','green');
xlabel('$\theta$ ','Interpreter','latex')
ylabel('$K$ (cm/s)','Interpreter','latex')
grid on
ax1.XTick = [0 0.1 0.2 0.3 0.4 0.5];

exportgraphics(gcf,'Van_Genutchen.pdf','ContentType','vector')


%% 2.0 - Domain Discretization
% Hydrus Example
Soil_Depth = 100; % cm
nz = 100; % Number of nodes
dt = 1; % Seconds
tfinal = 24; % Hours
alpha = 0; % Angle between flow direction and vertical axis


% Dussailant Example
% Soil_Depth = 100; % cm
% nz = 40; % Number of nodes
% dt = 1; % Seconds
% tfinal = 24; % Hours
% alpha = 0; % Angle between flow direction and vertical axis

% Z-axis discretization
dz = Soil_Depth/(nz-1); % cm
soil_depth = linspace(Soil_Depth,0,nz)*cos(deg2rad(alpha));
%% 3.0 - Boundary Conditions
% x is oriented upwards, such that x = 0 is the bottom and x = L is  the
% top

% In this model, only negative boundary conditions of heads are allowed (so
% far)

% Hydrus Example
% % Boundary Condition at x = 0
h0 = -100; % cm

% Boundary Condition at Interior Points => dx < x < L
hint = -100; % cm

% Boundary Condition at Top => x = L (upwards)
hf = -20; % cm

% Dussailant Example
% Boundary Condition at x = 0
% h0 = -1000; % cm
% 
% % Boundary Condition at Interior Points => dx < x < L
% hint = -1000; % cm
% 
% % Boundary Condition at Top => x = L (upwards)
% hf = -75; % cm

% Field Capacity
Sfc= n.^(-0.60*(2 + log10(Ks*86400))); % cm/day

%% 4.0 Pre-Allocation of Arrays
nt = (tfinal*3600)/dt - 1;

% States
theta = zeros(nz,nt);
head = zeros(nz,nt);
K = zeros(nz,nt);
q = zeros(nz,nt);

%% 5.0 - Assigning Boundary Conditions
% 1 => x = 0; nz => x = L

% ---- Heads ---- %
head(1,:) = h0; % Boundary Condition at x = 0
head(end,:) = hf; % Boundary Condition at x = L points
head(2:(end-1),1) = hint; % Boundary Condition at interior points

% ---- K ---- %
K(1,:) = K_function(head(1,:));
K(end,:) = K_function(head(end,:));

% ---- theta ---- %
theta(1,:) = theta_function(h0);
theta(2:(end-1),1) = theta_function(head(2:(end-1),1));
theta(end,:) = theta_function(hf);

% ---- Internal Nodes Vector ---- %
int_nodes = 2:(nz - 1);

% ---- Heads for i, i -1, and i + 1 ---- %
head_i = head(int_nodes,1); % Head at node i for internal nodes
head_prev = head(1:(end-2),1); % Head in internal nodes for i-1
head_post = head(3:(end),1); % Head in internal nodes for i+1

% ------ Unsaturated Hydraulic Conductivity ------ %
K(:,1) = K_function(head(:,1));

% ---- Previous and Post Hydraulic Conductivity ---- %

K_i = K_function(head_i); % Value of K for i
K(int_nodes,1) = K_i; % Value of K for i saved
K_prev = K_function(head_prev); % Value of K for i - 1
K_post = K_function(head_post); % Value of K for i + 1

% ---- Heads for i + 1/2 and i - 1/2 for Internal Nodes ---- %
head_half_post = (head_i + head_post)/2; % (h(i) + h(i+1))/2
head_half_prev = (head_i + head_prev)/2; % (h(i) + h(i-1))/2

% K for i + 1/2 and i - 1/2 for Internal Nodes
K_prev_half =(K_i + K_prev)/2; % Value of K for i - 1/2
K_post_half = (K_i + K_post)/2; % Value of K for i + 1/2

% K half for 1st and last node
K_post_1st = 1/2*(K(2,1) + K(1,1));
K_prev_last = 1/2*(K(end,1) + K(end-1,1));

% ------ Flows -------%
%%% 1st node
q(1,1) = -K_post_1st*((head(2,1) - head(1,1))/(dz) + 1);

%%% Internal Nodes
q(int_nodes,1) = (-K_post_half*dz.*((head_post - head_i)./dz + 1) ...
    -K_prev_half*dz.*((head_i - head_prev)./dz + 1))/(2*dz);

%%% Last node
q(end,1) = -K_prev_last*((head(end,1) - head(end-1,1))/(dz) + 1); % Still need to consider plant uptake

%% 5.0 - Main Loop

% - Richards Equation
% ptheta/pt = 1/px [K(ph/px + cos(alfa))] - S
% where p represents a partial derivative,
% theta is the volumumetric water content,
% t is the time
% x is the spatial coordinate (positive upward, starting from bottom)
% alfa is the angle between the flow direction and the vertical axis (i.e,
% alfa = 0 for vertical flow, 90 for horizontal flow, and between 0 and 90
% for inclined flow
% K is the unsaturated hydraulic conductivity given by
% K(h,x) = Ks(x) Kr(h,x)
% Where Kr is the unsaturated hydraulic conductivity



for k = 1:(nt-1) % 1 -> x = 0, nz => x= L

    % -------- Fully Explicit First Order Method --------- %
    % Water Content
    theta(2:(nz - 1),k+1) = theta(2:(nz - 1),k) + ...
        dt/dz*(K_post_half.*((head_post - head_i)/dz) - K_prev_half.*((head_i - head_prev)/dz))  ...
        +dt*cos(deg2rad(alpha))*((K_post_half- K_prev_half)/(dz));
    % Heads
    head(2:(nz - 1),k+1) = h_theta_function_2(theta(2:(nz - 1),k+1)); % Solving for h with known theta
    
    % 1 => x = 0; nz => x = L

    % ---- Heads for i, i -1, and i + 1 ---- %
    head_i = head(int_nodes,k+1); % Head at node i for internal nodes
    head_prev = head(1:(end-2),k+1); % Head in internal nodes for i-1
    head_post = head(3:(end),k+1); % Head in internal nodes for i+1

    % ------ Unsaturated Hydraulic Conductivity ------ %
    K(:,k+1) = K_function(head(:,k+1));

    % ---- Previous and Post Hydraulic Conductivity ---- %

    K_i = K_function(head_i); % Value of K for i
    K(int_nodes,1) = K_i; % Value of K for i saved
    K_prev = K_function(head_prev); % Value of K for i - 1
    K_post = K_function(head_post); % Value of K for i + 1

    % ---- Heads for i + 1/2 and i - 1/2 for Internal Nodes ---- %
    head_half_post = (head_i + head_post)/2; % (h(i) + h(i+1))/2
    head_half_prev = (head_i + head_prev)/2; % (h(i) + h(i-1))/2

    % K for i + 1/2 and i - 1/2 for Internal Nodes
    K_prev_half =(K_i + K_prev)/2; % Value of K for i - 1/2
    K_post_half = (K_i + K_post)/2; % Value of K for i + 1/2

    % K half for 1st and last node
    K_post_1st = 1/2*(K(2,k+1) + K(1,k+1));
    K_prev_last = 1/2*(K(end,k+1) + K(end-1,k+1));

    % ------ Flows -------%
    %%% 1st node
    q(1,k+1) = -K_post_1st*((head(2,k+1) - head(1,k+1))/(dz) + 1);

    %%% Internal Nodes
    q(int_nodes,k+1) = (-K_post_half*dz.*((head_post - head_i)./dz + 1) ...
        -K_prev_half*dz.*((head_i - head_prev)./dz + 1))/(2*dz);

    %%% Last node
    q(end,k+1) = -K_prev_last*((head(end,k+1) - head(end-1,k+1))/(dz) + 1); % Still need to consider plant uptake

    if sum(imag(((theta(2:(nz - 1),k+1))))) ~= 0  || sum(imag(((head(2:(nz - 1),k+1))))) ~= 0
        k/nt*100
        error('Instability')
    end
end

%% 6.0 - Post Processing
post_processing_richards

