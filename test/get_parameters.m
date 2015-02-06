function params_cell_out = get_parameters(varargin)
% function Cell = get_parameters(varargin) takes number of paired
% parameters as input, e.g. {'var_name_1', 'default_val_1', ...} and a filename
% where used parameters are read from. Returns a Cell {param1, param2, ...}
% of the parameter values.
% Created by Victoria Rudakova, 01/2015

n = size(varargin, 2);
assert(mod(n,2) == 1 && n >= 3);
fname = char(varargin(1,n));
if ~exist(fname)
    error('Data filename is incorrect or does not exist');
end
fid = fopen(fname);
data = textscan(fid, '%s', 'Delimiter', {'\n', '\t'}); % data{1}{i}
len = size(data{1},1);
params_cell_out = cell(2,len/2+(n-1)/2);
idx = 1;
for i = 1 : (n-1)/2
    pname = char(varargin(1,i*2-1));
    pval = char(varargin(1,i*2));
    for j = 1 : len/2
        pn = char(data{1}{j*2-1});
        pv = char(data{1}{j*2});
        if strcmp(pname, pn)
            pval = pv;
            break;
        end
    end
    params_cell_out{1,idx} = pname;
    params_cell_out{2,idx} = pval;
    idx = idx + 1;
end

for j = 1:len/2
    found = 0;
    pn = char(data{1}{j*2-1});
    pv = char(data{1}{j*2});
    for k = 1:len/2
        if strcmp(char(params_cell_out{1,k}), pn)
            found = 1;
            break;
        end
    end
    if (~found)
        params_cell_out{1,idx} = pn;
        params_cell_out{2,idx} = pv;
        idx = idx + 1;
    end
end

params_cell_out = params_cell_out(:,1:idx-1);
fclose(fid);

