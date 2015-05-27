% script to run the whole pipeline
% all the parameters are read from user provided general config file

function cryo3d(fname)

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
    'rotstart';'rotend';'transdelta';'transwidth';'alignims';};

len = length(cellparams);
for i=1:len
    s = cellparams{i};
    ind = find(strcmp(vars, s));
    if isempty(ind)
        error('No %s specified', s);
    else
        %if (ischar(C{2}{ind}))
        %    eval([s '=C{2}{ind}(2:end-1)']);
        %else
            eval([s '=C{2}{ind}']);
        %end        
    end
end

% check for format of pathdata variable
pathdata = pathdata(2:end-1);
slash = pathdata(end);
if (~isequal(slash, '\') && ~isequal(slash, '/'))
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % Windows
        pathdata = [pathdata '\'];
    else % Linux
        pathdata = [pathdata '/'];
    end
end

% check for format of pathout variable
pathout = pathout(2:end-1);
slash = pathout(end);
if (~isequal(slash, '\') && ~isequal(slash, '/'))
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % Windows
        pathout = [pathout '\'];
    else % Linux
        pathout = [pathout '/'];
    end
end

% check for format of pathcache variable
pathcache = pathcache(2:end-1);
slash = pathcache(end);
if (~isequal(slash, '\') && ~isequal(slash, '/'))
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % Windows
        pathcache = [pathout '\'];
    else % Linux
        pathcache = [pathout '/'];
    end
end

%% Preprocessing

structfile = [pathdata structfile(2:end-1)];
stackfile = [pathdata stackfile(2:end-1)];
if (length(ctffile)>2)
    ctffile = [pathdata ctffile(2:end-1)];
else
    ctffile = [];
end
if length(maskfile)>2
    maskfile = [pathdata maskfile(2:end-1)];
else
    maskfile = [];
end
structfile = init_volume(pathout, structfile, str2num(lpf), str2num(sigm), str2num(ds_ini));
coordfile = coordinate_axes(pathout, str2num(dtheta));
imfile = preprocess_images(pathout, stackfile, ctffile, str2num(ds_img));

%% Best Match

% Generate config file
configfile = generate_config(pathout, imfile, ctffile, str2num(substep), str2num(reconhalf), str2num(reconstartind), ...
    structfile, coordfile, maskfile, str2num(pf), str2num(pixfromedge), str2num(maxmem), str2num(numthreads), str2num(dispflag),...
    str2num(f), str2num(numruns), str2num(maxnumiter), str2num(convtol), str2num(normprojint), str2num(rotstart), str2num(rotstep), str2num(rotend), str2num(transmax), str2num(transdelta), str2num(transwidth), ...
    str2num(alignims));

best_match(pathout, configfile, str2num(caching), pathcache);
