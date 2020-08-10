%--------------------------------------------------------------------------
%
% gast: Greenwich Apparent Sidereal Time
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
%
% Output:
%   gstime    GAST in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function gstime = gast(Mjd_UT1,Arcs, MJD_J2000)

gstime = mod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1,Arcs, MJD_J2000), 2*pi);

