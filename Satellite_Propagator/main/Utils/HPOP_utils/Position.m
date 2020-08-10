%--------------------------------------------------------------------------
%
% Position: Position vector (r [m]) from geodetic coordinates
%           (Longitude [rad], latitude [rad], altitude [m])
%
% Last modified:   2018/01/04   M. Mahooti
%
%--------------------------------------------------------------------------
function r = Position(lon, lat, h,R_Earth ,f_Earth)

R_equ = R_Earth;
f     = f_Earth;

e2     = f*(2-f);  % Square of eccentricity
CosLat = cos(lat); % (Co)sine of geodetic latitude
SinLat = sin(lat);

% Position vector 
N = R_equ/sqrt(1-e2*SinLat*SinLat);

r(1) = (N+h)*CosLat*cos(lon);
r(2) = (N+h)*CosLat*sin(lon);
r(3) = ((1-e2)*N+h)*SinLat;

