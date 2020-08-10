% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 6th Jan 2018

% Function
% Input: CPF observations with time stamp,observation_number
% Output: Full state vector in ECI frame

% This code is has lines copied from Author Meysam Mahooti's implementation

function [Y0,error_pos] = init_state(cpfdata_fun,obs_no,eopdata,AuxParam,mod_acc,Cnm, Snm)
sigma_X = 10;  % [m]
sigma_Y = 10;  % [m]
sigma_Z = 10;  % [m]


nobs = length(cpfdata_fun.mjdutc); % no of observations

n2 = round(nobs/2);
% 3 observations are taken into account
Mjd1 = cpfdata_fun.mjdutc(1);
Mjd2 = cpfdata_fun.mjdutc(n2);
Mjd3 = cpfdata_fun.mjdutc(end);

r1 = [cpfdata_fun.xeci(1), cpfdata_fun.yeci(1), cpfdata_fun.zeci(1)];
r2 = [cpfdata_fun.xeci(n2), cpfdata_fun.yeci(n2), cpfdata_fun.zeci(n2)];
r3 = [cpfdata_fun.xeci(end), cpfdata_fun.yeci(end), cpfdata_fun.zeci(end)];

[v2,~,~,~,~] = hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3);

Y0_apr = [r2, v2];
% Y0_apr_mee = eci2mee(GM_Earth, r2, v2);

Mjd0 = cpfdata_fun.mjdutc(obs_no);%Mjd1 - 1/86400; % Set Mjd0 one second before Mjd1

options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
%
[~,~,~,hour,minu,sec] = invjday(Mjd0+2400000.5);
Tend = hour*60*60 + minu*60 + sec;
[~,~,~,hour,minu,sec] = invjday(Mjd2+2400000.5);
Tstart = hour*60*60 + minu*60 + sec;


[UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eopdata,Mjd2,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
AuxParam.Mjd_TT = Mjd2 + TT_UTC/86400;
[~,yout] = radau(mod_acc,[Tstart  Tend],Y0_apr,options);
Y0 = yout(end,:)';
%
A = zeros(nobs*2,6);
b = zeros(nobs*2,1);
w = zeros(nobs*2,nobs*2);

yPhi = zeros(42,1);
Phi  = zeros(6);
%%
for iterat = 1:5
    for i=1:nobs
        % Time increment and propagation
        Mjd_UTC = cpfdata_fun.mjdutc(i);          % Julian Date
        t       = (Mjd_UTC-Mjd0)*86400;       % Time since epoch [s]
        [UT1_UTC,TAI_UTC,~,~,~,~] = IERS(eopdata,Mjd_UTC,'l');
        [~, ~, ~, TT_UTC, ~] = timediff(UT1_UTC, TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;
%         AuxParam.Mjd_UTC = Mjd_UTC;
%         AuxParam.Mjd_TT = Mjd_TT;
        AuxParam.Mjd_UTC_Var = Mjd_UTC;
        VarEqn = @(x,yPhi) VarEqn_mod(x, yPhi,AuxParam, eopdata,Cnm, Snm);

        for ii=1:6
            yPhi(ii) = Y0(ii);
            for j=1:6  
                if (ii==j) 
                    yPhi(6*j+ii) = 1; 
                else
                    yPhi(6*j+ii) = 0;
                end
            end
        end
        
        yPhi = DEInteg(VarEqn,0,t,1e-13,1e-6,42,yPhi);

        % Extract state transition matrices        
        for j=1:6
            Phi(:,j) = yPhi(6*j+1:6*j+6);
        end
        
        % Change here
        [~,~,~,hour,minu,sec] = invjday(Mjd0+2400000.5);
        Tstart = hour*60*60 + minu*60 + sec;
        [~,~,~,hour,minu,sec] = invjday(Mjd_UTC+2400000.5);
        Tend = hour*60*60 + minu*60 + sec;
        
        if t ~= 0
              [~,yout] = radau(mod_acc,[Tend (Tend+t)],Y0,options);
            Y = yout(end,:)';
        else 
            Y = Y0;
        end

        % Observations and partials        
        dXdY0 = [1, 0, 0,zeros(1,3)]*Phi;
        dYdY0 = [0, 1, 0,zeros(1,3)]*Phi;
        dZdY0 = [0, 0, 1,zeros(1,3)]*Phi;
        
        % i th observed position vector
        ri = [cpfdata_fun.xeci(i), cpfdata_fun.yeci(i), cpfdata_fun.zeci(i)];
        
        % Accumulate least-squares system
        A(3*i-2:3*i,:) = [dXdY0;dYdY0;dZdY0];        
        b(3*i-2:3*i) = [ri(1)-Y(1) ; ri(2)-Y(2) ; ri(3)-Y(3)];        
        w(3*i-2:3*i,3*i-2:3*i) = diag([1/sigma_X^2,1/sigma_Y^2,1/sigma_Z^2]);   
    end
    % Solve least-squares system
    dY0 = inv(A'*w*A)*A'*w*b;
    % Correct epoch state
    Y0 = Y0 + dY0;
end
% Y0 contains the updated full state values in ECI frame for the point of
% interest
r_true = [cpfdata_fun.xeci(obs_no),cpfdata_fun.yeci(obs_no),cpfdata_fun.zeci(obs_no)];
error_pos  = norm(r_true' - Y0(1:3));
end
