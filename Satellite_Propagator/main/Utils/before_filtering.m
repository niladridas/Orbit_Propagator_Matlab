
function [mee_start,full_state,obs_states,Mjd_UTC,del_t] = before_filtering(observation,...
    cpfdata_compiled,AuxParam, GM_Earth,eopdata,MJD_J2000,Arcs,PC,Cnm, Snm) 
    [year,month,day,~,~,~] = invjday(observation.mjdutc(1)+2400000.5);
    Mjd_UTC = Mjday(year, month, day, 0, 0, 0);
    AuxParam.Mjd_UTC = Mjd_UTC;
    % Read CPF data a
    c_i = findclosest(observation.mjdutc(1),cpfdata_compiled.mjdutc);
    del_t = (cpfdata_compiled.mjdutc(c_i) - observation.mjdutc(1))*86400;
    % Some points from cpf data are taken for initial state calculation
    for h = 1:3
        tmp_cpf.date = cpfdata_compiled.date((c_i-h):(c_i+h+2),:);
        tmp_cpf.mjdutc = cpfdata_compiled.mjdutc((c_i-h):(c_i+h+2),:);
        tmp_cpf.xeci = cpfdata_compiled.xeci((c_i-h):(c_i+h+2),:);
        tmp_cpf.yeci = cpfdata_compiled.yeci((c_i-h):(c_i+h+2),:);
        tmp_cpf.zeci = cpfdata_compiled.zeci((c_i-h):(c_i+h+2),:);
        % Least Square Based Orbit Determination Technique is used here
        mod_acc = @(x,Y) accel_mee_mod(x,Y,MJD_J2000 , eopdata,  AuxParam,Arcs,PC,Cnm, Snm);
        [full_state,error_pos] = init_state(tmp_cpf,h+1,eopdata, AuxParam,mod_acc,Cnm, Snm);
        if (error_pos<=0.5)
            break
        end
    end
    mee_start.elements = coe2mee(eci2coe(full_state',GM_Earth));
    mee_start.mjdutc = cpfdata_compiled.mjdutc(c_i);  
    %%
    Nobs = size(observation.mjdutc,1);
    obs_states = zeros(Nobs,2+3); 
    obs.pressure = observation.P; % mbar
    obs.temp = observation.T; % Kelvin
    obs.humid = observation.Rh; % percent
    for k = 1:Nobs
        obs_states(k,1) = observation.mjdutc(k); % mjdutc of the observation
        obs_states(k,2) = observation.range(k); % 
        obs_states(k,3) = obs.pressure(k);
        obs_states(k,4) = obs.temp(k);
        obs_states(k,5) = obs.humid(k);
    end
end