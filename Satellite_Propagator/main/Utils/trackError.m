%--------------------------------------------------------------------------
% Version  Date         Author      Mod History
%--------------------------------------------------------------------------
% 1        25 March 2018  Vedang D.   Inital version
%--------------------------------------------------------------------------

% Input: 
% cpfdata - compiled cpf data [m], [mjd]
% cpf_index - index of reference state in cpfdata
% EST_POS - estimated ECI position vector at same time as cpf_index [m]
% Output: 
% along_track - along track error [m]
% cross_track - cross track error, normal to orbit plane [m]
% inplane_cross_track - cross track error, in the orbit plane [m]
% Note: This function uses following functions from IOD Module by  M.
% Mahooti: hgibbs, unit

% All units are [m], [s].

function [along_track, cross_track_normal, cross_track_inplane] = trackError(cpfdata,cpf_index,EST_POS)

    % SAT_Const
    
    GM_Earth    = 398600.4418e9;                % [m^3/s^2]; WGS-84

    if (cpf_index <1)||(cpf_index>numel(cpfdata.mjdutc))
        error('Invalid input Index for CPD data.')
    elseif (cpf_index ==1)||(cpf_index==numel(cpfdata.mjdutc))
        error('Input index can not be at extremeties of CPF data.')
    end
    
    R1 = [cpfdata.xeci(cpf_index-1),cpfdata.yeci(cpf_index-1),cpfdata.zeci(cpf_index-1)]';
    R2 = [cpfdata.xeci(cpf_index),cpfdata.yeci(cpf_index),cpfdata.zeci(cpf_index)]'; % point of interest
    R3 = [cpfdata.xeci(cpf_index+1),cpfdata.yeci(cpf_index+1),cpfdata.zeci(cpf_index+1)]';
    MJD1 = cpfdata.mjdutc(cpf_index-1);
    MJD2 = cpfdata.mjdutc(cpf_index);
    MJD3 = cpfdata.mjdutc(cpf_index+1);
    
   [V2,theta,theta1,copa,err] = hgibbs(R1,R2,R3,MJD1,MJD2,MJD3);
   
    % Angular momentum vector
    H_vec  =  cross(R2,V2);
    % angular momentum magnitude
    h = norm(H_vec);
    
    % Eccentricity vector equation at R2
    E_vec = cross(V2,H_vec)/GM_Earth - R2/norm(R2);
    % eccentricity magnitude
    e = norm(E_vec);
    
    % Using energy equation to calculate semi major axis
    a = -GM_Earth/(norm(V2)^2 - 2*GM_Earth/norm(R2));
    % semi minor axis
    b = a*sqrt(1-e^2);
    
    % Coordinate transformation - orbital elements unchanged
    % X Y Z - Original Frame
    X = unit(E_vec)';
    Z = unit(H_vec)';
    Y = cross(Z,X);
    % x y z -  Transformed Frame
    % x = [1 0 0]'; y = [0 1 0]'; z = [0 0 1]';
    % Transformation matrix
    dcm  = [X,Y,Z]'; % Takes from original to new frame
    
    % Transformed vectors
    e_vec = dcm*E_vec;
    r2 = dcm*R2;
    est_pos = dcm*EST_POS;
    
    % focus is Earth [0 0 0]'. Find center of ellipse
    center = [0 0 0]' - a*e_vec;
    
    % Find closest point on ellipse to the projection of estimated position
    % cost function
    fun = @(closePt)((closePt(1)-est_pos(1))^2 + (closePt(2)-est_pos(2))^2 );
    closePt0 = r2;
    [closePtSol, fval] = fmincon(fun, closePt0, [], [], [], [], [], [], @(closePt)ellipseEqn(closePt,center,a,b)); %, options);
    
    % Calculate track  errors
    syms xe 
    dydx = -b^2*(xe - center(1))/a^2/sqrt(b^2*(1-(xe-center(1))^2/a^2));
    expr = sqrt(1+dydx^2);
    along_track = abs(vpaintegral(expr,xe,r2(1),closePtSol(1)));
    cross_track_normal = abs(est_pos(3) - r2(3));
    cross_track_inplane = norm(r2(1:2) - est_pos(1:2));
end