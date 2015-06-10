% Function to read .star file data. Only works for one data block/file
% Usage:
%       [labels, data] = read_star_file(starfile)
% Inputs:
%       starfile = Name of .star file
% Outputs:
%       labels = Cell array of label names
%       data = Cell array of data values (numbers or text)

% Created by Nicha C. Dvornek, 06/2015

function [labels, data] = read_star_file(starfile)

[labels,numhlines] = read_star_header(starfile);

fid = fopen(starfile);
if fid < 0
    disp('Error: Cannot open .star file for read');
    return
end

% Read in all the data as text
data = textscan(fid,repmat('%s ',[1,length(labels)]),'HeaderLines',numhlines);

% Fix data type for numeric data values
for i = 1:length(labels)
    temp = str2double(data{i});
    if ~isnan(temp)
        data{i} = temp;
    end
end

fclose(fid);