% Function to make an initial 3D volume from a reference volume
% Parameters, example:
% structfile = 'C:\Users\vicrucann\Home\server\sample-db\test_init_model.mrc';

function passed = init_volume(structfile, lpf, sigma)

passed = 0;
if (nargin < 3)
    sigma = 1;
end
if (nargin < 2)
    lpf = 60;
end

%% Access the necessary functions
addpath(fullfile(cd, '../src/preprocessing'));
addpath(fullfile(cd, '../src/best_match'));
addpath(fullfile(cd, '../src/mrc'));

%% Load
disp('Load reference structure');
if strcmp(structfile(end-2:end),'mrc')
    structure = double(ReadMRC(structfile));
    [~,h] = ReadMRC(structfile,1,-1);
    voxelsize = h.pixA;
else
    disp('Error: File extension must be .mrc');
    return
end

%% Low pass filter
disp(['Low Pass Filtering at ' num2str(lpf) ' Angstrom']);
if lpf > 0
    structure = lpf_at_res_sigma(structure,voxelsize,lpf,sigma);
    % OR do we want to use a sharp filter?
%     fc = 1/(lpf/voxelsize);
%     structure = SharpFilt(structure,fc,0);
end

%% Save volume
disp('Save initial volume');
savefile = [structfile(1:end-4) '_lpf_' num2str(lpf) 'A.mrc'];
writeMRC(structure,voxelsize,savefile);
passed = 1;