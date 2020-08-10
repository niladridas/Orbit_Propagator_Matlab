% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 25nd Nov 2017


% This code is modified from Accel.m in HPOP


function dY = accel_mee_mod(x,Y,MJD_J2000 , eopdata,  AuxParam,Arcs ,PC,Cnm, Snm)
   
     % constants independent of the satellite
    %SAT_Const
    
    [UT1_UTC,TAI_UTC,x_pole,y_pole,~,~] = IERS(eopdata,AuxParam.Mjd_UTC+x/86400,'l');
    [~, ~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = AuxParam.Mjd_UTC+x/86400+TT_UTC/86400;
    Mjd_UT1 = AuxParam.Mjd_UTC+x/86400+UT1_UTC/86400;

    P = PrecMatrix(MJD_J2000,Mjd_TT,Arcs, MJD_J2000);
    N = NutMatrix(Mjd_TT,Arcs, MJD_J2000);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1,Arcs, MJD_J2000) * T;

    [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
     r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE405(AuxParam.Mjd_UTC+x/86400,PC);
%------------------------------------------------------------------------------------------------
%     Y = coe2eci(mee2coe(Y1'),GM_Earth)';

    % Acceleration due to harmonic gravity field
%     a = AccelHarmonic(AuxParam.Mjd_UTC+x/86400,r_Sun,r_Moon,Y(1:3),E);

   
    
    % Acceleration due to harmonic gravity field
    a = AccelHarmonic_mod(Y(1:3), E, AuxParam.n, AuxParam.m,Cnm, Snm);

    dY = [Y(4:6);a];

end
