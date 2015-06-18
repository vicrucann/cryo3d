% Function to read .star file header. Only works for one data block/file
% Usage:
%       [labels, numhlines] = read_star_header(starfile)
% Inputs:
%       starfile = Name of .star file
% Outputs:
%       labels = Cell array of label names
%       numhlines = number of header lines

% Created by Nicha C. Dvornek, 06/2015

function [labels,numhlines] = read_star_header(starfile)

% Check if can open file
fid = fopen(starfile);
if fid < 0
    disp('Error: Cannot open .star file for read');
    return
end

% Read in lines until get to first label
tline = fgetl(fid);
numhlines = 1;
% ~strcmp(tline(1),'_')
while numel(tline) == 0 || ~strcmp(tline(1),'_')
    tline = fgetl(fid);
    numhlines = numhlines + 1;
end

% Read in the label names
temp = textscan(tline,'%s');
labelind = 1;
labels{labelind} = temp{1}{1,1};
while strcmp(tline(1),'_')
    tline = fgetl(fid);
    if strcmp(tline(1),'_')
        temp = textscan(tline,'%s');
        labelind = labelind + 1;
        labels{labelind} = temp{1}{1,1};
        numhlines = numhlines + 1;
    end
end

fclose(fid);