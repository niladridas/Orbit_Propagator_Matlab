function dY = accel_mee(x,Y1,sysfun,eopdata,AuxParam,MJD_J2000,...
    AU,GM_Earth, GM_Sun, GM_Moon, GM_Mercury, GM_Venus, GM_Mars,...
    GM_Jupiter, GM_Saturn, GM_Uranus, GM_Neptune, GM_Pluto, P_Sol,...
    jb2008_path,Arcs,PC,Cnm, Snm,R_Earth, f_Earth,c_light,R_Sun, R_Moon)
  
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
    Y = coe2eci(mee2coe(Y1'),GM_Earth)';
%     if norm(Y(1:3,1))<= R_Earth
%         error("The sattelite has crash landed");
%     end
    AuxParam.n      = 10;
    AuxParam.m      = 10;
    a = AccelHarmonic_mod(Y(1:3), E, AuxParam.n, AuxParam.m,Cnm, Snm);
    if (AuxParam.sun)
        a = a + AccelPointMass(Y(1:3),r_Sun,GM_Sun);
    end
    if (AuxParam.moon)
        a = a + AccelPointMass(Y(1:3),r_Moon,GM_Moon);
    end
    if (AuxParam.planets)
        a = a + AccelPointMass(Y(1:3),r_Mercury,GM_Mercury);
        a = a + AccelPointMass(Y(1:3),r_Venus,GM_Venus);
        a = a + AccelPointMass(Y(1:3),r_Mars,GM_Mars);
        a = a + AccelPointMass(Y(1:3),r_Jupiter,GM_Jupiter);
        a = a + AccelPointMass(Y(1:3),r_Saturn,GM_Saturn);
        a = a + AccelPointMass(Y(1:3),r_Uranus,GM_Uranus);    
        a = a + AccelPointMass(Y(1:3),r_Neptune,GM_Neptune);
        a = a + AccelPointMass(Y(1:3),r_Pluto,GM_Pluto);
    end
    if (AuxParam.sRad)
        a = a + AccelSolrad(Y(1:3),r_Earth,r_Moon,r_Sun,r_SunSSB, ...
            AuxParam.area_solar,AuxParam.mass,AuxParam.Cr,P_Sol,AU,R_Earth, R_Sun, R_Moon);
    end
    if (AuxParam.drag)
        % Atmospheric density
        Omega = 7292115.8553e-11+4.3e-15*( (AuxParam.Mjd_UTC+x/86400-MJD_J2000)/36525 );   
        %%
        % Convert from MJD_UTC time to year, day, hour, mins, secs
        time_MJD  = AuxParam.Mjd_UTC+x/86400;
        [year, mnth, day, hr, min, sec] = invjday(time_MJD+2400000.5);
        % use year, mnth and day to find the day of the year
        dayfyr = dayofyear(year, mnth, day);
        input_time = sprintf(' %i %i %i %i %2.10f',year, dayfyr, hr, min, sec);
        % We still have to add the height , lat and log information
        [XLON, lat, height] = Geodetic(T*Y(1:3),R_Earth, f_Earth);
        [UT1_UTC,~,~,~,~,~] = IERS(eopdata,AuxParam.Mjd_UTC+x/86400,'l');
        MJD_UT1 = AuxParam.Mjd_UTC+x/86400 + UT1_UTC/86400.0;
        GWRAS = gmst(MJD_UT1);
        lon = (mod(GWRAS + XLON, 2*pi))*180/pi; % deg 
        lat = lat*180/pi; % deg 
        ht = height/1000; % Km
        display(ht)
        tmp_loc = sprintf(' %f %f %f',ht, lat, lon);
        command = strcat(jb2008_path,input_time,tmp_loc);
        [status,cmdout] = sysfun(command);
        if status ~= 0
            error('The fortran code did not run')
        end
        cmdout = strtrim(cmdout);
        cmdout = strrep(cmdout,'D','e');
        dens = str2double(cmdout);
        if isnan(dens)
            error('Density is not read correctly')
        end
        a = a + AccelDrag(dens,Y(1:3),Y(4:6),T,AuxParam.area_drag,AuxParam.mass,AuxParam.Cd,Omega);
    end
    if (AuxParam.Relativity)
        a = a + Relativity(Y(1:3),Y(4:6),c_light, GM_Earth);
    end
    eci_val = Y';%coe2eci(mee2coe(Y1),GM_Earth);
    S = dot(a,eci_val(1,1:3)')/norm(eci_val(1,1:3));
    temp1 = cross(eci_val(1,1:3)',eci_val(1,4:6)');
    normal_vec = temp1/norm(temp1);
    N = dot(a,normal_vec);
    temp_2 = cross(normal_vec,eci_val(1,1:3)');
    C_unit = temp_2/norm(temp_2);
    C = dot(a,C_unit);
    S = S + GM_Earth/(norm(eci_val(1:3))^2);
    p = Y1(1);
    f = Y1(2);
    g = Y1(3);
    h = Y1(4);
    k = Y1(5);
    L = Y1(6);
    s = sqrt(1+h^2+k^2);
    w = 1 + f*cos(L) + g*sin(L) ;
    mu = GM_Earth;
    dp = 2*p*C*sqrt(p/mu)/w;
    df = sqrt(p/mu)*(S*sin(L)+(C*((w+1)*cos(L)+f)/w)-(g*N*(h*sin(L)-k*cos(L))/w));
    dg = sqrt(p/mu)*(-S*cos(L)+(C*((w+1)*sin(L)+g)/w)+(f*N*(h*sin(L)-k*cos(L))/w));
    dh = sqrt(p/mu)*s^2*N*cos(L)/(2*w);
    dk = sqrt(p/mu)*s^2*N*sin(L)/(2*w);
    dL = sqrt(mu*p)*(w/p)^2 + sqrt(p/mu)*(N*(h*sin(L)-k*cos(L))/w);
    dY = [dp;df;dg;dh;dk;dL];
end
