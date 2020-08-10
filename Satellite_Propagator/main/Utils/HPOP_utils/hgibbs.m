%--------------------------------------------------------------------------
%
%  hgibbs: implements the herrick-gibbs approximation for orbit
%          determination, and finds the middle velocity vector for the 3
%          given position vectors.
%
%  inputs:
%    r1          - ijk position vector #1         m
%    r2          - ijk position vector #2         m
%    r3          - ijk position vector #3         m
%    jd1        - julian date of 1st sighting    days from 4713 bc
%    jd2        - julian date of 2nd sighting    days from 4713 bc
%    jd3        - julian date of 3rd sighting    days from 4713 bc
%
%  outputs:
%    v2          - ijk velocity vector for r2     m/s
%    theta       - angl between vectors           rad
%    error       - flag indicating success        'ok',...
%
%--------------------------------------------------------------------------
function [v2,theta,theta1,copa,error] = hgibbs(r1,r2,r3,jd1,jd2,jd3)

%SAT_Const
GM_Earth = 3.986004418*10^14;

error =  '          ok';
magr1 = norm(r1);
magr2 = norm(r2);
magr3 = norm(r3);

tolangle= 0.01745329251994;
dt21= (jd2-jd1)*86400;
dt31= (jd3-jd1)*86400;
dt32= (jd3-jd2)*86400;

p = cross(r2,r3);
pn = unit(p);
r1n = unit(r1);
copa = asin(dot(pn,r1n));

if ( abs(dot(r1n,pn)) > 0.017452406 )
    error= 'not coplanar';
end

theta  = angl(r1,r2);
theta1 = angl(r2,r3);

if ( (theta > tolangle) || (theta1 > tolangle) )  
    error= '   angl > 1';
end

term1 = -dt32*( 1/(dt21*dt31) + GM_Earth/(12*magr1*magr1*magr1) );
term2 = (dt32-dt21)*( 1/(dt21*dt32) + GM_Earth/(12*magr2*magr2*magr2) );
term3 =  dt21*( 1/(dt32*dt31) + GM_Earth/(12*magr3*magr3*magr3) );

v2 =  term1*r1 + term2* r2 + term3* r3;