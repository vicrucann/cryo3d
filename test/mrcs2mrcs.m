% script to split / convert one *.mrcs file to a new *.mrcs stack file(-s)

clear; clc;

addpath(fullfile(cd, '../src/mrc'));
vox_size = 1.32;

s1 = ReadMRC('G:\20150205_sdp\stackfile_selected_1.mrcs');
s2 = ReadMRC('G:\20150205_sdp\stackfile_selected_2.mrcs');
s12 = cat(3, s1, s2);
len = size(s12,3);
perms = randperm(len);

fst = (1:len/2);
snd = (len/2+1:len);
split1 = s12(:,:,fst);
split2 = s12(:,:,snd);

writeMRC(split1, vox_size, 'G:\20150205_sdp\stackfile_selected_12_fst.mrcs');
writeMRC(split2, vox_size, 'G:\20150205_sdp\stackfile_selected_12_snd.mrcs');

%original = 'G:\20150205_sdp\stackfile_selected_12.mrcs';

%all = ReadMRC(original);
%s = size(all,3);
% generate indices for splitting
%odd = (1:2:s);
%even = (2:2:s);
%split1 = all(:,:,odd);
%split2 = all(:,:,even);

%writeMRC(split1, vox_size, 'G:\20150205_sdp\stackfile_selected_12_odd.mrcs');
%writeMRC(split2, vox_size, 'G:\20150205_sdp\stackfile_selected_12_even.mrcs');