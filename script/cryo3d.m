% script to run the whole pipeline
% all the parameters are read from user provided general config file

function cryo3d(fname)

if nargin < 1
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % Windows
        fname = 'config_example_win.txt';
    else
        fname = 'config_example_lin.txt';
    end
end

fileID = fopen(fname);
if fileID < 0
    disp('Error: Cannot open config file for read');
    return
end

C = textscan(fileID,'%s %s','delimiter',{'=', ';'},'commentstyle','%');

strinds = strfind(C{2},'%');
cinds = find(~cellfun('isempty',strinds));
for i = 1:length(cinds)
    C{2}{cinds(i)}(strinds{cinds(i)}:end) = [];
end
vars = strtrim(C{1});

cellparams = {'pathout';'pathdata';'caching';'pathcache';'dtheta';'maskfile';'maxmem';'numthreads';'rotstep';'transmax';...
    'structfile';'stackfile';'ctffile';'lpf';'sigm';'ds_ini';'ds_img';...
    'substep';'reconhalf';'reconstartind';'pf';'pixfromedge';'dispflag';...
    'f';'numruns';'maxnumiter';'convtol';'normprojint';...
    'rotstart';'rotend';'transdelta';'transwidth';'alignims';...
    'ipaddrs';'login';'ppath';'remmat';'varmat';'sleeptime';'resfold';'printout';};

len = length(cellparams);
for i=1:len
    s = cellparams{i};
    ind = find(strcmp(vars, s));
    if isempty(ind)
        error('No %s specified', s);
    else
        eval([s '=C{2}{ind}']);      
    end
end

% check for format of path-like variables
pathdata = eval(pathdata);
pathdata = fixslash(pathdata);

pathout = eval(pathout);
pathout = fixslash(pathout);

ppath = eval(ppath);
ppath = fixslash(ppath);

resfold = eval(resfold);
resfold = fixslash(resfold);

pathcache = eval(pathcache);
pathcache = fixslash(pathcache);

%% Preprocessing

structfile = init_volume(pathout, eval(structfile), str2num(lpf), str2num(sigm), str2num(ds_ini));
coordfile = coordinate_axes(pathout, str2num(dtheta));
imfile = preprocess_images(pathout, eval(stackfile), eval(ctffile), str2num(ds_img));

%% Best Match

% Generate config file
configfile = generate_config(pathout, imfile, ctffile, str2num(substep), str2num(reconhalf), str2num(reconstartind), ...
    structfile, coordfile, eval(maskfile), str2num(pf), str2num(pixfromedge), str2num(maxmem), str2num(numthreads), str2num(dispflag),...
    str2num(f), str2num(numruns), str2num(maxnumiter), str2num(convtol), str2num(normprojint), str2num(rotstart), str2num(rotstep), ...
    str2num(rotend), str2num(transmax), str2num(transdelta), str2num(transwidth), str2num(alignims),...
    eval(ipaddrs), eval(login), ppath, eval(varmat), str2num(sleeptime), resfold, str2num(printout));

best_match(pathout, configfile, str2num(caching), pathcache);
end

function foldname = fixslash(foldname)
slash = foldname(end);
if (~isequal(slash, '\') && ~isequal(slash, '/'))
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % Windows
        foldname = [foldname '\'];
    else % Linux
        foldname = [foldname '/'];
    end
end
end
