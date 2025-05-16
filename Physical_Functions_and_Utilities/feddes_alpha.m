function alpha = feddes_alpha(h, h_lim_upper, h_lim_lower)
% FEDDES_ALPHA calculates evaporation reduction factor (alpha)
% based on surface pressure head (h) following Feddes' approach.
%
% INPUTS:
%   h            - Surface pressure head [cm] (scalar or vector)
%   h_lim_upper  - Upper threshold for full evaporation [cm] (wetter, e.g., -10 cm)
%   h_lim_lower  - Lower threshold for zero evaporation [cm] (drier, e.g., -1500 cm)
%
% OUTPUT:
%   alpha        - Evaporation reduction factor [0-1]

alpha = zeros(size(h));  % Initialize alpha with zeros

% No reduction region (wet condition)
alpha(h >= h_lim_upper) = 1;

% Linear reduction region
linear_zone = (h < h_lim_upper) & (h > h_lim_lower);
alpha(linear_zone) = (h(linear_zone) - h_lim_lower) ./ (h_lim_upper - h_lim_lower);

% Full reduction region (dry condition), alpha remains zero

end
