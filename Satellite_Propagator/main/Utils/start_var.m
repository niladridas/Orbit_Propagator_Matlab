%% before_loop 
load('Sat_params.mat');
Aux_fileread;
load('crddata_complied.mat');
load('cpfdata_compiled.mat');
passes = clustercrd(crddata_compiled);
no_passes = size(passes,2); 
gap_pass = zeros(no_passes-1,1);
for i = (1:no_passes-1)
    gap_pass(i) = (passes(i+1).mjdutc(1)- passes(i).mjdutc(end));
end
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
station_ecefxyz = [-2386279.6636;-4802356.3304;3444883.4298];
[stlat,stlong,stalt] = ecef2geo(station_ecefxyz', 0);
st_latlongalt = [stlat  stlong stalt/1000]; % rad and Km % lat in North 
lambda = 532/1000; % in micrometer