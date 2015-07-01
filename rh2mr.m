%
% calculate the specific humidity [g/kg]
%
% s=rh2spez(T,p,rh);
%
% T: Temperature [K]
% rh: relative humidity [%]
% p: pressure [hPa]
%
% source: CIMO guide http://www.wmo.int/pages/prog/www/IMOP/publications/CIMO-Guide/CIMO%20Guide%207th%20Edition,%202008/Part%20I/Chapter%204.pdf
%
% haa, 2013-10-16
% haa/rjs, 2014-09-08: documentation fixed, order is T, p, rh 
%
function [qw, qi] = rh2mr(T,p,rh)

% temperature in C
t = T-273.15;

% pressure factor
f = 1.0016 + 3.15e-6 * p - 0.074./p;

% saturation vapour pressure in pure phase over water
ew = 6.112*exp(17.62*t./(243.12 + t));

% saturation vapour pressure of humid air over water
eww = f .* ew;

% saturation vapour pressure in pure phase over ice
ei = 6.112*exp(22.46*t./(272.62 + t));

% saturation vapour pressure of humid air over ice
eii = f .* ei;

% vapor pressure for rh over water
eew = rh .* eww / 100;

% vapor pressure for rh over ice
eei = rh .* eii / 100;

b = 0.62198;

% mixing ratio for rh over water in g/kg
qw = b ./ (p./eew - 1) * 1000;

% mixing ratio for rh over ice in g/kg
qi = b ./ (p./eei - 1) * 1000;