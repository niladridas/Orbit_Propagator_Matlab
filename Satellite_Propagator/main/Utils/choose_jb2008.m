if ismac
    sysfun = @(x) system(x);
    jb2008_path_all = dir(which('JB2008.'));
    jb2008_folder = jb2008_path_all.folder;
    jb2008_path = char(strcat({'cd'}, {' '},jb2008_folder,' && ./JB2008'));
elseif isunix
    sysfun = @(x) system(x);
    jb2008_path_all = dir(which('JB2008linux.'));
    jb2008_folder = jb2008_path_all.folder;
    jb2008_path = char(strcat({'cd'}, {' '},jb2008_folder,' && ./JB2008linux'));
elseif ispc
    sysfun = @(x) jsystem(x);
    [jb2008_path_all,~,~] = fileparts(which('JB2008.exe'));
    jb2008_path = char(strcat({'cd'}, {' '},jb2008_path_all,' && JB2008.exe'));
else
    disp('Platform not supported')
end