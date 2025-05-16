function ET0 = hargreaves_ET0(Tavg, Tmax, Tmin, latitude, julian_day)
% hargreaves_ET0 - Calculates potential evapotranspiration using the Hargreaves method
%
% Syntax: ET0 = hargreaves_ET0(Tavg, Tmax, Tmin, latitude, julian_day)
%
% Inputs:
%   Tavg       - Average temperature (°C)
%   Tmax       - Maximum temperature (°C)
%   Tmin       - Minimum temperature (°C)
%   latitude   - Latitude in decimal degrees
%   julian_day - Julian day (1 to 365/366)
%
% Output:
%   ET0        - Potential evapotranspiration (mm/day)

if any(latitude < -90 | latitude > 90)
    error('Latitude must be between -90 and 90 degrees.');
end

% Convert latitude to radians
phi = deg2rad(latitude);

% Constants
Gsc = 1360; % Solar constant [J / m2 / s]

% Calculate solar declination (delta)
delta = 0.409 * sin((2 * pi * julian_day / 365) - 1.39);

% Calculate inverse relative distance Earth-Sun (dr)
dr = 1 + 0.033 * cos((2 * pi * julian_day) / 365);

% Calculate sunset hour angle (omega_s)
omega_s = acos(-tan(phi) * tan(delta));

if any(Tmax < Tmin)
    error('Tmax should be greater than or equal to Tmin.');
end

% Calculate extraterrestrial radiation (Ra) [J / m2 / s]
Ra = (1 / pi) * Gsc * dr .* (omega_s .* sin(phi) .* sin(delta) + cos(phi) .* cos(delta) .* sin(omega_s));

% Converting from J / m2 / s to mm / day
Lambda = 2.45; % MJ / kg (latent heat of vaporization)

Ra = (Ra * 86400) / (Lambda * 1e6); % mm / day

% Hargreaves equation
ET0 = 0.0023 * (Tavg + 17.8) .* ((Tmax - Tmin) .^ 0.5) .* Ra; % mm / day
end
