% script to split / convert one *.mrcs file to a new *.mrcs stack file(-s)

clear; clc;
original = 'G:\20150205_sdp\stackfile_selected_12.mrcs';
vox_size = 1.32;

addpath(fullfile(cd, '../src/mrc'));

all = ReadMRC(original);
s = size(all,3);
% generate indices for splitting
odd = (1:2:s);
even = (2:2:s);
split1 = all(:,:,odd);
split2 = all(:,:,even);

writeMRC(split1, vox_size, 'G:\20150205_sdp\stackfile_selected_12_odd.mrcs');
writeMRC(split2, vox_size, 'G:\20150205_sdp\stackfile_selected_12_even.mrcs');