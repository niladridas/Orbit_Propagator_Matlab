%--------------------------------------------------------------------------
%
% Gradient: Computes the gradient of the Earth's harmonic gravity field 
%
% Inputs:
%   Mjd_UT      Modified Julian Date (Universal Time)
%   r           Satellite position vector in the true-of-date system
%   n,m         Gravity model degree and order
% 
% Output:
%   G           Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function G = Gradient ( r, U, n_max, m_max,Cnm, Snm )

d = 1; % Position increment [m]

G = zeros(3);
dr = zeros(3,1);

% Gradient
for i=1:3
    % Set offset in i-th component of the position vector
    dr(:) = 0;
    dr(i) = d;
    % Acceleration difference
    da = AccelHarmonic_mod( r+dr/2,U, n_max, m_max ,Cnm, Snm) - ...
         AccelHarmonic_mod( r-dr/2,U, n_max, m_max,Cnm, Snm );
    % Derivative with respect to i-th axis
    G(:,i) = da/d;
end

