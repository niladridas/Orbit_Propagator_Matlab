%--------------------------------------------------------------------------
%
% LTCMatrix: Transformation from Greenwich meridian system to local tangent
%            coordinates
%
% Inputs:
%   lon      -Geodetic East longitude [rad]
%   lat      -Geodetic latitude [rad]
%   
%   Output:
%   M        -Rotation matrix from the Earth equator and Greenwich meridian
%             to the local tangent (East-North-Zenith) coordinate system
%--------------------------------------------------------------------------
function M = LTCMatrix(lon, lat)

M = iauRy(-lat,eye(3))*iauRz(lon,eye(3));

for j=1:3
    Aux=M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;
end

