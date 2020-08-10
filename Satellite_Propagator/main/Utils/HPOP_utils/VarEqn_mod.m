%--------------------------------------------------------------------------
%
% VarEqn: Computes the variational equations, i.e. the derivative of the
%         state vector and the state transition matrix
%
% Input:
%   x           Time since epoch AuxParam.Mjd_0 in [s]
%   yPhiS       (6+36+6)-dim vector comprising the state vector (y), the
%               state transition matrix (Phi) and the sensitivity matrix
%               in column wise storage order
%
% Output:
%   yPhip      Derivative of yPhi
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function yPhip = VarEqn_mod(x, yPhi,AuxParam, eopdata,Cnm, Snm)


[UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eopdata,AuxParam.Mjd_UTC_Var+x/86400,'l');
Mjd_UT1 = AuxParam.Mjd_UTC_Var+x/86400+UT1_UTC/86400;

% Transformation matrix
U = R_z(gmst(Mjd_UT1));
 
% State vector components
r = yPhi(1:3);
v = yPhi(4:6);
Phi = zeros(6);

% State transition matrix
for j=1:6
    Phi(:,j) = yPhi(6*j+1:6*j+6);
end

% Acceleration and gradient

a = AccelHarmonic_mod(r, U, AuxParam.n_a, AuxParam.m_a,Cnm, Snm);
G = Gradient(r, U, AuxParam.n_G, AuxParam.m_G,Cnm, Snm);

% Time derivative of state transition matrix
yPhip = zeros(42,1);
dfdy = zeros(6);

for i=1:3
  for j=1:3
    dfdy(i,j) = 0.0;               % dv/dr(i,j)
    dfdy(i+3,j) = G(i,j);          % da/dr(i,j)
    if ( i==j ) 
        dfdy(i,j+3) = 1; 
    else
        dfdy(i,j+3) = 0;           % dv/dv(i,j)
    end
    dfdy(i+3,j+3) = 0;             % da/dv(i,j)
  end
end

Phip = dfdy*Phi;

% Derivative of combined state vector and state transition matrix
for i=1:3
  yPhip(i)   = v(i);               % dr/dt(i)
  yPhip(i+3) = a(i);               % dv/dt(i)
end

for i=1:6
  for j=1:6
    yPhip(6*j+i) = Phip(i,j);      % dPhi/dt(i,j)
  end
end
  
