%% Read the observation file 
global AuxParam
crd_source = 'starlette_20170401_MONL.npt';
% Extracted from Dataset No. 1568624
if strcmp(crd_source((end-3):end),'.frd')~=1 && strcmp(crd_source((end-3):end),'.npt')~=1
    error('This file does not have a CRD data file extension');
else
    disp('The file extension looks good but the CRD file needs checking');
end
% CRD checking at this point
fidcrd = fopen(crd_source);
[station,stationepochtimescale,target,cosparid,sicid,...
    noradid,spaceepochtimescale,targettype,datatype,...
    year,month,day,hour,minute,sec,endyear,endmonth,...
    endday,endhour,endminute,endsec, refcorrect, cmcorrect,...
    recampcorrect,stsysdelcorrect,satsysdelcorrect,...
    rngtypeind,dataqual,endofsession,endoffile,...
    rangerecord,normrecord,metdata,usedwavelngth,pointing_data] = crd_parse(fidcrd);
fclose(fidcrd);
%% Read the referece states from CPF file
cpf_source = 'starlette_cpf_20170401.hts';
fidcpf = fopen(cpf_source);
% extract postion vectors
cpfdata = getcpf(fidcpf);
fclose(fidcpf);
%%
if isempty(rangerecord)==true
    full_observation = normrecord;
elseif isempty(normrecord)==true
    full_observation = rangerecord;
else
    error('Neither full rate data nor normal pointing data parsed');
end
Mjd_UTC = Mjday(str2double(year), str2double(month), str2double(day), 0, 0, 0);
AuxParam.Mjd_UTC = Mjd_UTC;
%
Jd_DATE = datestr(datetime(Mjd_UTC+2400000.5,'convertfrom','juliandate'));
%
% Filtering out the same day prediction in ECI frame
% The CPF file often contains prediction of the months before it and after
% it
cpfmoddata.date = [];
cpfmoddata.mjdutc = [];
cpfmoddata.xeci = [];
cpfmoddata.yeci = [];
cpfmoddata.zeci = [];

for i = 1:length(cpfdata.date)
    if strcmp(cpfdata.date(i,:),Jd_DATE)==1
        cpfmoddata.date = [cpfmoddata.date;cpfdata.date(i,:)];
        cpfmoddata.mjdutc = [cpfmoddata.mjdutc;cpfdata.mjdutc(i)];
        cpfmoddata.xeci = [cpfmoddata.xeci;cpfdata.xeci(i)];
        cpfmoddata.yeci = [cpfmoddata.yeci;cpfdata.yeci(i)];
        cpfmoddata.zeci = [cpfmoddata.zeci;cpfdata.zeci(i)];
    end
end
%
save('global_data.mat','eopdata', 'Cnm', 'Snm', 'SOLdata', 'DTCdata',...
    'APdata', 'PC', 'AuxParam','usedwavelngth');
save('crd_data.mat','full_observation','metdata','usedwavelngth','rngtypeind');
save('cpf_data.mat','cpfmoddata');
