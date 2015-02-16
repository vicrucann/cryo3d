% Function to read config file for making initial model

function [structfile,sampdeg,coordfile,ctffile,savename,pf,addnoise,SNR_dB,pixfromedge] = read_config_file_init_model(filename)

fileID = fopen(filename);
if fileID < 0
    disp('Error: Cannot open config file for read');
    return
end
C = textscan(fileID,'%s %s','delimiter','=','commentstyle','%');

strinds = strfind(C{2},'%');
cinds = find(~cellfun('isempty',strinds));
for i = 1:length(cinds)
    C{2}{cinds(i)}(strinds{cinds(i)}:end) = [];
end

vars = strtrim(C{1});

ind = find(strcmp(vars,'structfile'));
if isempty(ind)
    ind = find(strcmp(vars,'initprojfile'));
    if isempty(ind)
        disp('Error: No initial structure specified (structfile)');
        return
    else
        structfile = C{2}{ind};
    end
else
    structfile = C{2}{ind};
end

ind = find(strcmp(vars,'sampdeg'));
if isempty(ind)
    sampdeg = 0;
else
    sampdeg = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'coordfile'));
if isempty(ind)
    if sampdeg == 0
        disp('Error: No coordinates file specified (coordfile)');
        return
    end
    coordfile = '';
else
    coordfile = C{2}{ind};
end

ind = find(strcmp(vars,'ctffile'));
if isempty(ind)
    ctffile = '';
else
    ctffile = C{2}{ind};
end

ind = find(strcmp(vars,'savename'));
if isempty(ind)
    s = strfind(structfile,'\');
    if isempty(s)
        s = strfind(structfile,'/');
    end
    savename = structfile(s(end)+1:end-4);
else
    savename = C{2}{ind};
end

ind = find(strcmp(vars,'pf'));
if isempty(ind)
    pf = 0;
else
    pf = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'addnoise'));
if isempty(ind)
    addnoise = 0;
else
    addnoise = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'SNR_dB'));
if isempty(ind)
    SNR_dB = 0;
else
    SNR_dB = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'pixfromedge'));
if isempty(ind)
    pixfromedge = 1;
else
    pixfromedge = str2double(C{2}{ind});
end
