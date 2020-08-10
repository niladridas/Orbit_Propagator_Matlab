% Author: Vedang Deshpande
% Email: vedang.deshpande@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
%--------------------------------------------------------------------------
% Version  Date         Author      Mod History
%--------------------------------------------------------------------------
% 1        12 Jan 2018  Vedang D.   Inital version
%--------------------------------------------------------------------------
% Inputs:
%   mjd_utc = mjd time in utc
%   lon0 = longitude of station [rad]
%   lat0 = lattitude of station [rad]
%   alt0 = altitude of station from ellipsoid [m]
%   r_eci = ijk position of satellite in ECI frame [m]
% Outputs:
%   Elev = Elevation angle of satellite, [rad]
%   d_eta_r = d/dE(eta_r) at E
% 
% Calls: This function calls following functions from Initial Orbit
%   Determination (Least Squares Method) module by Meysam Mahooti
%   Position(), LTCMatrix(), IERS(), timediff(), gmst(), R_z(), AzElPa().
%   This also uses global eopdata variable. 
% 
% Called by: Main code to calculate eta_r as function of elevation angle

function Elev = elevation(mjd_utc,lon0,lat0,alt0,r_eci,eopdata,R_Earth ,f_Earth)


Rs = Position(lon0, lat0, alt0,R_Earth ,f_Earth)'; % Topocentric position of station [m]
E = LTCMatrix(lon0,lat0);

[UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eopdata,mjd_utc,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
Mjd_TT = mjd_utc + TT_UTC/86400;
Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

% Topocentric coordinates        
theta = gmst(Mjd_UT1);                  % Earth rotation [rad]
U = R_z(theta);
s = E*(U*r_eci-Rs);                     % Topocentric position [m]

% Observations and partials
[Azim, Elev, dAds, dEds] = AzElPa(s);   % Azimuth, Elevation [rad]
end