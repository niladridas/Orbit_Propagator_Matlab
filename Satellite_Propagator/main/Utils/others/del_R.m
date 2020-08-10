% Author: Niladri Das (niladridas@tamu.edu)
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% This code is adapted from eta_r_inv.m
% Input: 
% 1. lambda: wavelength in micrometer
% 2. lat: latitude in radian
% 3. H: elevation in Km
% 4. P0: atmospheric pressure in millibars
% 5. T0: temperature in degree Kelvins
% 6. Rh0: relative humidity
% 7. E: elevation
% Output:
% 1. x : delta R

function x = del_R(lambda,lat,H,P0,T0,Rh0,E)
    e0 = (Rh0*6.11/100)*10^(7.5*(T0 - 273.15)/(237.3 + T0 - 273.15)); % eqn 2-5
    f_lambda = 0.9650 + 0.0164/(lambda^2) + 0.000228/(lambda^4);
    f_lat_alt = 1 -0.0026*cos(2*lat) - 0.00031*H;
    K = 1.163 - 0.00968*cos(2*lat) - 0.00104*T0 + 0.00001435*P0;
    A = 0.002357*P0 + 0.000141*e0;
    B = (1.084e-8)*P0*T0*K + (4.734e-8)*(P0^2)*2/(T0*(3-1/K));
    % Define
    C = A+B; D = B/C; M = f_lambda/f_lat_alt;
    x = M*C/(sin(E) + D/(sin(E) +0.01) );
end