% Function to read .star file data. Only works for one data block/file.
% Usage:
%       [labels, data] = read_star_file(starfile,numrows)
% Inputs:
%       starfile = Name of .star file
%       numrows = (Optional) Number of rows of data to read
% Outputs:
%       labels = Cell array of label names
%       data = Cell array of data values (numbers or text)

% Created by Nicha C. Dvornek, 06/2015

function [labels, data] = read_star_file(starfile,numrows)

if nargin < 2
    numrows = -1;
end

% Read header to get label names and number of header lines
[labels,numhlines] = read_star_header(starfile);

% Check if can open file
fid = fopen(starfile);
if fid < 0
    disp('Error: Cannot open .star file for read');
    return
end

% Read in first line of data to get appropriate data format
firstline = textscan(fid,repmat('%s ',[1,length(labels)]),1,'HeaderLines',numhlines);
fmt = [];
for i = 1:length(firstline)
    temp = str2double(firstline{i});
    if isnan(temp)
        fmt = [fmt '%s '];
    else
        fmt = [fmt '%f '];
    end
end

% Read in the data
frewind(fid);
if numrows > 0
    data = textscan(fid,fmt,numrows,'HeaderLines',numhlines);
else
    data = textscan(fid,fmt,'HeaderLines',numhlines);
end

fclose(fid);