% Function to read config file for fast_multi_ref code

% Created by Nicha C. Dvornek 03/2015
% Last modified 06/2015

function [imfile,imreconfile,ctffile,maxmem,numthreads,dispflag,substep,reconhalf,reconstartind,normprojint,numruns,maxnumiter,rotstart,rotstep,rotend,transmax,transdelta,transwidth,convtol,t,numstruct] = read_config_file_multiref(filename)

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

ind = find(strcmp(vars,'imfile'));
if isempty(ind)
    disp('Error: No image stack specified (imfile)');
    return
else
    imfile = C{2}{ind};
end

ind = find(strcmp(vars,'imreconfile'));
if isempty(ind)
    imreconfile = imfile;
else
    imreconfile = C{2}{ind};
end

ind = find(strcmp(vars,'ctffile'));
if isempty(ind)
    ctffile = '';
else
    ctffile = C{2}{ind};
end

ind = find(strcmp(vars,'maxmem'));
if isempty(ind)
    [~,sys] = memory;
    maxmem = sys.PhysicalMemory.Available/1024/1024 - 1024;
else
    maxmem = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'numthreads'));
if isempty(ind)
    pc = parcluster('local');
    numthreads = pc.NumWorkers;
else
    numthreads = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'dispflag'));
if isempty(ind)
    dispflag = 0;
else
    dispflag = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'substep'));
if isempty(ind)
    substep = 0;
else
    substep = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'reconhalf'));
if isempty(ind)
    reconhalf = 0;
else
    reconhalf = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'reconstartind'));
if isempty(ind)
    reconstartind = 0;
else
    reconstartind = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'normprojint'));
if isempty(ind)
    normprojint = 1;
else
    normprojint = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'numruns'));
if isempty(ind)
    numruns = 2;
else
    numruns = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'maxnumiter'));
if isempty(ind)
    maxnumiter = 15;
else
    maxnumiter = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'rotstart'));
if isempty(ind)
    rotstart = 0;
else
    rotstart = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'rotstep'));
if isempty(ind)
    rotstep = 2;
else
    rotstep = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'rotend'));
if isempty(ind)
    rotend = 359.9;
else
    rotend = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'transmax'));
if isempty(ind)
    transmax = 4;
else
    transmax = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'transdelta'));
if isempty(ind)
    transdelta = 1;
else
    transdelta = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'transwidth'));
if isempty(ind)
    transwidth = 1;
else
    transwidth = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'convtol'));
if isempty(ind)
    convtol = 0.01;
else
    convtol = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'t'));
if isempty(ind)
    t = 0.0001;
else
    t = str2double(C{2}{ind});
end

ind = find(strcmp(vars,'numstruct'));
if isempty(ind)
    numstruct = 2;
else
    numstruct = str2double(C{2}{ind});
end

fclose(fileID);