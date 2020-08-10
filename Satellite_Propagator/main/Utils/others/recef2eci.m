%--------------------------------------------------------------------------
%
% ECEF2ECI: Transforms Earth Centered Earth Fixed (ECEF) coordinates to 
%           Earth Centered Inertial (ECI) coordinates
%
% Original Author:   2015/08/12   M. Mahooti
% Modified by:       2017/11/24   Vedang D.
% Modified by:       2018/03/02   Niladri Das
%--------------------------------------------------------------------------
function Y = recef2eci(MJD_UTC, Y0,MJD_J2000,eopdata,Arcs)
 

[UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eopdata,MJD_UTC, 'l' );
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);

MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400; 

% ICRS to ITRS transformation matrix and derivative
P      = PrecMatrix(MJD_J2000,MJD_TT,Arcs, MJD_J2000); % IAU 1976 Precession
N      = NutMatrix(MJD_TT,Arcs, MJD_J2000);            % IAU 1980 Nutation
Theta  = GHAMatrix(MJD_UT1,Arcs, MJD_J2000);           % Earth rotation
Pi     = PoleMatrix(x_pole,y_pole);    % Polar motion

S = zeros(3);
S(1,2) = 1; S(2,1) = -1;               % Derivative of Earth rotation 
Omega = 7292115.8553e-11+4.3e-15*( (MJD_UTC-MJD_J2000)/36525 );
dTheta = Omega*S*Theta;          % matrix [1/s]

U      = Pi*Theta*N*P;                 % ICRS to ITRS transformation
dU     = Pi*dTheta*N*P;                % Derivative [1/s]

% Transformation from WGS to ICRS
r = U'*Y0(1:3,1);
%v = U'*Y0(4:6)' + dU'*Y0(1:3)';
Y = r; %[r;v];
end
