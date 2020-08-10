%--------------------------------------------------------------------------
% Version  Date         Author      Mod History
%--------------------------------------------------------------------------
% 1        22 March 2018  Vedang D.   Inital version
%--------------------------------------------------------------------------

% Input: 
% xECI - column of x ECI coordinates [m]
% yECI - column of y ECI coordinates [m]
% zECI - column of z ECI coordinates [m]
% mjdutc - column of corresponding mjdutc times
% Output: 
% lat - column of lattitude angles [degrees]
% lon - column of longitude angles [degrees]

function [lat,lon] = getLatLon(xECI,yECI,zECI,mjdutc)
    % Initialize
    lon = zeros(numel(mjdutc),1);
    lat = lon;
    for i = 1:numel(mjdutc)
        rECEF = ECItoECEF(mjdutc(i)+2400000.5,[xECI(i),yECI(i),zECI(i)]');
        [lat(i,1),lon(i,1),~] = ecef2geo(rECEF', 0);    
    end
    % convert to degrees
    lat =rad2deg(lat);
    lon =rad2deg(lon);
end