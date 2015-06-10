% Function to read data corresponding to given label names from .star file.
% Only works for one data block/file.
% Usage:
%       data = read_star_data_for_labels(starfile,labels,numrows)
% Inputs:
%       starfile = Name of .star file
%       labels = Cell array of label names
%       numrows = Number of rows of data to read
% Outputs:
%       data = Cell array of data values corresponding to labels

% Created by Nicha C. Dvornek, 06/2015

function data = read_star_data_for_labels(starfile,labels,numrows)

if nargin < 3
    numrows = -1;
end

% Read header to get label names and number of header lines
[all_labels,numhlines] = read_star_header(starfile);

% Check that given labels are in the star file
cols = zeros(length(labels));
for i = 1:length(labels)
    tf = strcmp(labels{i},all_labels);
    inds = find(tf == 1);
    if isempty(inds)
        disp('Error: Label name does not exist in .star file');
        return;
    elseif length(inds) > 1
        disp('Warning: Label name occurs multiple times in .star file. Using first occurance');
    end
    cols(i) = inds(1);
end

% Check if can open file
fid = fopen(starfile);
if fid < 0
    disp('Error: Cannot open .star file for read');
    return
end

% Read in first line of data to get appropriate data format
firstline = textscan(fid,repmat('%s ',[1,length(all_labels)]),1,'HeaderLines',numhlines);
datatypes = cell(length(firstline),1);
for i = 1:length(firstline)
    temp = str2double(firstline{i});
    if isnan(temp)
        datatypes{i} = '%s ';
    else
        datatypes{i} = '%f ';
    end
end

% Read in the columns corresponding to label names
frewind(fid);
fmt = [];
ind = 0;
for i = 1:length(labels)
    fmt = [fmt repmat('%*s ',[1,cols(i)-ind-1]) datatypes{cols(i)}];
    ind = cols(i);
end
fmt = [fmt '%*[^\n]'];
if numrows > 0
    data = textscan(fid,fmt,numrows,'HeaderLines',numhlines);
else
    data = textscan(fid,fmt,'HeaderLines',numhlines);
end

fclose(fid);