%% before_loop 
load('Sat_params.mat');
Aux_fileread;
fid = fopen('InitialStarlette.txt','r');
tline = fgetl(fid);
AuxParam.area_solar = str2double(tline(49:51));
tline = fgetl(fid);
AuxParam.area_drag = str2double(tline(38:44));
tline = fgetl(fid);
AuxParam.mass = str2double(tline(19:23));
tline = fgetl(fid);
AuxParam.Cr = str2double(tline(5:7));
tline = fgetl(fid);
AuxParam.Cd = str2double(tline(5:7));
fclose(fid);
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 1;
AuxParam.n      = 10;
AuxParam.m      = 10;
AuxParam.n_a    = 10;
AuxParam.m_a    = 10;
AuxParam.n_G    = 10;
AuxParam.m_G    = 10;

radius_limit = 1;
options = rdpset('RelTol',1e-13,'AbsTol',1e-16);

choose_jb2008