% Function to make an initial 3D volume from a reference volume
% Required parameters, example:
% pathout = 'G:\workspace\';
% structfile = 'G:\db-sample\test_init_model.mrc';
% lpf = 60; Cutoff for low-pass filter in Angstrom
% Optional parameters, example:
% sigma = 1; Controls width of Gaussian edge for low-pass filter
% ds = 2; Factor by which to downsample the volume

function outfile = init_volume(pathout, structfile, lpf, sigma, ds)

if (nargin < 5)
    ds = 1;
end
if (nargin < 4)
    sigma = 1;
end
if (nargin < 3)
    disp('ERROR: Not enough input parameters');
    return;
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
    structure = DownsampleGeneral(structure,imgSx);
    voxelsize = voxelsize*factor;
end

%% Save volume
disp('Save initial volume');
f_path = strfind(structfile, '\');
if (isempty(f_path))
    f_path = strfind(structfile, '/');
end
savefile = structfile(max(f_path)+1:end);
savefile = [savefile(1:end-4) '_lpf' num2str(lpf) 'A'];
if ds > 1
    savefile = [savefile '_ds' num2str(ds)];
end
writeMRC(structure,voxelsize,[pathout savefile '.mrc']);
outfile = [pathout savefile];