% Function to be used in fast_best_match.m to get all the initial model
% related things from just an initial structure and config parameters

% Created by Nicha C. Dvornek, 01/2015


function savename = make_init_model(configfile, pathout)

[structfile,maskfile,sampdeg,coordfile,ctffile,savename,pf,...
    addnoise,SNR_dB,pixfromedge] = read_config_file_init_model(configfile);

savename = [pathout '/' savename];

%% Load
disp('Load structure, ctfs, and projection coordinates');
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
else
    disp('Error: File extension must be .mat or .mrc');
    return
end

if isempty(maskfile)
    structmask = ones(size(structure));
else
    structmask = double(ReadMRC(maskfile));
end

if isempty(ctffile)
    ctfs = ones(size(structure,1),size(structure,2));
else
    load(ctffile,'ctfs');
    ctfs = double(ctfs);
end

if sampdeg > 0
    coord_axes = create_approx_uniform_axes_whole_sphere(sampdeg,0,1);
else
    load(coordfile,'coord_axes');
end

%% Project
disp('Projecting');
if pf
    ctfs = abs(ctfs);
end

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
proj_struct = zeros(size(proj));
maskim = mask.bin(:,:,ceil(end/2));
for j = 1:numproj
    if addnoise
        temp = make_noisy_im(real(ifft2(ifftshift(proj(:,:,j)))),SNR_dB).*maskim;
    else
        temp = real(ifft2(ifftshift(proj(:,:,j)))).*maskim;
    end
    proj_struct(:,:,j) = temp;
end

disp('Saving');
if isempty(ctffile)
    savename = [savename '_noctf'];
else
    savename = [savename '_' num2str(size(ctfs,3)) 'ctf'];
end
savename = [savename '_' num2str(size(coord_axes,2)) 'd'];
if addnoise
    savename = [savename '_' num2str(SNRdB) 'dB'];
end
if pf
    savename = [savename '_pf'];
end
savename = [savename '.mat'];
save(savename,'-v7.3','maskim','ctfs','data_axes','coord_axes','mask','structure','proj_struct','structmask');
