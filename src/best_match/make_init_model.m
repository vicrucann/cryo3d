% Function to be used in fast_best_match.m to get all the initial model
% related things from just an initial structure and config parameters

% Created by Nicha C. Dvornek, 01/2015


function savename = make_init_model(configfile, pathout)

[structfile,sampdeg,coordfile,ctffile,ctfvar,savename,...
    structvoxelres,ds,lpf,pf,sigma,addnoise,SNR_dB,...
    pixfromedge] = read_config_file_init_model(configfile);

savename = [pathout savename];

%% Load
disp('Load structure');
if strcmp(structfile(end-2:end),'mat')
    load(structfile,'structure_final');
    if exist('structure_final','var')
        structure = double(structure_final);
    else
        s = load(structfile);
        if isfield(s,'vol')
            structure = double(s.vol);
        elseif isfield(s,'x')
            structure = double(s.x);
        else
            varname = fieldnames(s);
            structure = double(s.(varname{1}));
        end
    end
elseif strcmp(structfile(end-2:end),'mrc')
    structure = double(ReadMRC(structfile));
    [~,h] = ReadMRC(structfile,1,-1);
    structvoxelres = h.pixA;
else
    disp('Error: File extension must be .mat or .mrc');
    return
end
if isempty(ctffile)
    ctfs = ones(size(structure,1)/ds,size(structure,2)/ds);
else
    if ~isempty(ctfvar)
        s = load(ctffile,ctfvar);
        ctfs = double(s.(ctfvar));
    else
        s = load(ctffile);
        varname = fieldnames(s);
        structure = double(s.(varname{1}));
    end
end

if sampdeg > 0
    coord_axes = create_approx_uniform_axes_whole_sphere(sampdeg,0,1);
else
    load(coordfile,'coord_axes');
end

%% Low pass and Resize
disp('Low Pass Filtering and Resizing');
if lpf > 0
    structure = lpf_at_res_sigma(structure,structvoxelres,lpf,sigma);
end
n = size(ctfs,1);
nm = size(structure,1);
structure = Crop(DownsampleGeneral(structure,round(nm/ds)),n,0,0);
prc = prctile(structure(structure~=0),90);

%% Project
disp('Projecting');
numdir = size(coord_axes,2);
if pf
    ctfs = abs(ctfs);
end
numctfs = size(ctfs,3);
imsize = size(structure,1);

%if (~isempty (gcp('nocreate')) ) % matlab 2014, may not be needed
%    delete(gcp('nocreate'));
%end
%poolobj = parpool;
%if matlabpool('size') > 0 % Matlab 2013
%    matlabpool close
%end
%matlabpool open local
mask=get_mask_struct_ncd(size(structure),pixfromedge);
structure=structure.*mask.bin;
[proj,data_axes]=project_in_all_directions_w_ctf_par(structure,mask,coord_axes,ctfs);
%delete(poolobj); % matlab 2014, may not be needed
%matlabpool close % matlab 2013

% fourier -> spatial
numproj = size(proj,3);
numpix = size(proj,1) * size(proj,2);
proj_struct = zeros(size(proj));
data = zeros(numproj,numpix);
maskim = mask.bin(:,:,ceil(end/2));
for j = 1:numproj
    if addnoise
        temp = make_noisy_im(real(ifft2(ifftshift(proj(:,:,j)))),SNR_dB).*maskim;
    else
        temp = real(ifft2(ifftshift(proj(:,:,j)))).*maskim;
    end
    proj_struct(:,:,j) = temp;
    data(j,:) = temp(:);
end

disp('Saving');
if isempty(ctffile)
    savename = [savename '_noctf'];
else
    savename = [savename '_' num2str(numctfs) 'ctf'];
end
if ds ~= 1
    savename = [savename '_ds' num2str(ds)];
end
savename = [savename '_' num2str(size(coord_axes,2)) 'd'];
if lpf == 0
    savename = [savename '_nolpf'];
else
    savename = [savename '_lpf' num2str(lpf)];
end
if addnoise
    savename = [savename 'n'];
end
if pf
    savename = [savename '_pf'];
end
savename = [savename '.mat'];
if exist(savename,'file')
    save(savename,'-append','maskim','ctfs','data_axes','coord_axes','mask','structure','proj_struct');
else
    save(savename,'-v7.3','maskim','ctfs','data_axes','coord_axes','mask','structure','proj_struct');
end
