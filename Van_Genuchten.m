% This function performs Van-Genuchten model (direct formula) for soil 
% water retention curve. This function is called in the main program
% (Water_Retention.m) multiple times to plot SWRC for different types of
% soil.
%--------------------------------------------------------------------------
% Inputs:
% psi     - suction pressure ([L] or cm of water)
% theta_s - saturated water content [L^3L^?3]
% theta_r - residual water content [L^3L^?3]
% alpha   - related to the inverse of the air entry suction,
% alpha>0 ([L^?1], or cm^?1)
% n       - measure of the pore-size distribution,n>1 (dimensionless)
% Outputs:
% theta   -  water retention curve [L^3L^?3]
%--------------------------------------------------------------------------
function [ theta ] = Van_Genuchten(theta_s,theta_r, alpha, n )

psi = linspace(1e0,1e6,1000000);
m = 1-(1/n);
theta = theta_r + ((theta_s-theta_r)./(1+(alpha*psi).^n).^m);
end

