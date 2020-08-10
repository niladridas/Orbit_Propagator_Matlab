% Auxilliary file read
% global eopdata Cnm Snm SOLdata DTCdata APdata PC AuxParam

fid = fopen('eop19620101.txt','r'); % Has a specific format. Before updating check the format
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);
%
load DE405Coeff.mat
PC = DE405Coeff;

Cnm = zeros(71,71);
Snm = zeros(71,71);
fid = fopen('egm96','r');
for n=0:70
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read space weather data
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
line = fgetl(fid);
SOLdata_mod = [];
while ischar(line)
   C =  strsplit(line);
   SOLdata_mod = [SOLdata_mod;[str2double(C(1)),str2double(C(2)),str2double(C(3)),str2double(C(4)),...
       str2double(C(5)),str2double(C(6)),str2double(C(7)),str2double(C(8)),str2double(C(9)),...
       str2double(C(10)),str2double(C(11))]] ;
   line = fgetl(fid);
end
% SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
% SOLdata = SOLdata(1:11,1:end);
SOLdata = SOLdata_mod';
fclose(fid);

% READ GEOMAGNETIC STORM DTC VALUE
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
% DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
% fclose(fid);
line = fgetl(fid);
DTCdata_mod = [];
while ischar(line)
   C =  strsplit(line);
   DTCdata_mod = [DTCdata_mod;[str2double(C(2)),str2double(C(3)),str2double(C(4)),...
       str2double(C(5)),str2double(C(6)),str2double(C(7)),str2double(C(8)),str2double(C(9)),...
       str2double(C(10)),str2double(C(11)),str2double(C(12)),str2double(C(13))...
       ,str2double(C(14)),str2double(C(15)),str2double(C(16)),str2double(C(17)),str2double(C(18))...
       ,str2double(C(19)),str2double(C(20)),str2double(C(21)),str2double(C(22)),str2double(C(23))...
       ,str2double(C(24)),str2double(C(25)),str2double(C(26)),str2double(C(27))]] ;
      line = fgetl(fid);
end
fclose(fid);
DTCdata = DTCdata_mod';
%
% save('DTCdata.mat','DTCdata');

% read space weather data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);
%
%
AuxParam = struct('area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

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