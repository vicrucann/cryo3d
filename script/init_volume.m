% Function to make an initial 3D volume from a reference volume
% Parameters, example:
% structfile = 'C:\Users\vicrucann\Home\server\sample-db\test_init_model.mrc';

function passed = init_volume(structfile, lpf, sigma, ds)

passed = 0;
if (nargin < 4)
    ds = 1;
end
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

if ds > 1
    imgSx = floor(size(structure,1)/ds);
    if mod(imgSx,2) ~= 0
        imgSx = imgSx + 1;
    end
    factor = size(structure,1) / imgSx;
    structure = DownsampleGeneral(structure,imgSx,1);
    voxelsize = voxelsize*factor;
end

%% Save volume
disp('Save initial volume');
savefile = [structfile(1:end-4) '_lpf' num2str(lpf) 'A_ds' num2str(ds) '.mrc'];
writeMRC(structure,voxelsize,savefile);
passed = 1;