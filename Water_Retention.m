% Water retention curve or the soil moisture characteristics provides the 
% relation between the water content(theta) and water potential(psi). It is 
% used to measure the field capacity, soil quality,soil water storage, in 
% short for soil management. Different models can be used to study the soil
% water retention curve in a specific field. 
% The below program shows water retention curve for different soils 
% (sand, silt, loam silt, clay)using Van Genuchten Model which is widely 
% used for SWRC (soil water retention curve).   
%--------------------------------------------------------------------------
close all;
clear
clc
%--------------------------------------------------------------------------
% Inputs:
% psi - suction pressure ([L] or cm of water)
% Output:
% theta - water retention curve [L^3L^?3] which is same as function
% (Van_Genuchten.m) output-name. 
% The values and the necessary information considered below were referred 
% from the links
% https://en.wikipedia.org/wiki/Water_retention_curve
% https://de.wikipedia.org/wiki/Datei:Wrc.svg
% since they need to be calculated based on the measures taken on field or
% through lab testing.
%--------------------------------------------------------------------------
% For Sand:
%--------------------------------------------------------------------------
psi = linspace(1e0,1e6,1000000);
%-------------------------------
[theta] = Van_Genuchten(0.370687,0.043019, 0.087424, 1.57535 );
plot(theta,psi, 'displayname', 'Sand', 'linewidth', 2)
hold on
%--------------------------------------------------------------------------
% For Silt:
%--------------------------------------------------------------------------
[theta] = Van_Genuchten(0.421256,0, 0.003405, 1.34475 );
plot(theta,psi, 'displayname', 'Silt', 'linewidth', 2)
hold on
%--------------------------------------------------------------------------
% For Loam Silt:
%--------------------------------------------------------------------------
[theta] = Van_Genuchten(0.421217,0, 0.013345, 1.12614 );
plot(theta,psi, 'displayname', 'Loam Silt', 'linewidth', 2)
hold on
%--------------------------------------------------------------------------
% For Clay:
%--------------------------------------------------------------------------
[theta] = Van_Genuchten(0.550541,0, 0.006812, 1.08155 );
plot(theta,psi, 'displayname', 'Clay', 'linewidth', 2)
hold on
%--------------------------------------------------------------------------
title('Water Retention Curve')
xlabel('\theta', 'FontWeight', 'bold')
legend ('show')
set(gca, 'YScale', 'log')
axis([0 0.5 1e0 1e6])
grid on

%---------------------
% Second y-axis:
%---------------------
ay1 = gca;
ay1_pos = ay1.Position;             % position of first axes
ay2 = axes('Position',ay1_pos,...
    'YAxisLocation','right',...
    'Color','none');
ay2.XColor = 'none';
set(ay2, 'YLim', [0 6]);
%--------------------------------------------------------------------------
ylabel(ay1, '- \psi{m}/ hPa', 'FontWeight', 'bold')
ylabel(ay2, 'pF', 'FontWeight', 'bold')
%--------------------------------------------------------------------------