%--------------------------------------------------------------------------
%
% NutMatrix: Transformation from mean to true equator and equinox
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   NutMat    Nutation matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function NutMat = NutMatrix(Mjd_TT,Arcs, MJD_J2000)

% Mean obliquity of the ecliptic
ep = MeanObliquity(Mjd_TT, MJD_J2000);

% Nutation in longitude and obliquity
[dpsi deps] = NutAngles(Mjd_TT,Arcs, MJD_J2000);

% Transformation from mean to true equator and equinox
NutMat = R_x(-ep-deps)*R_z(-dpsi)*R_x(ep);

