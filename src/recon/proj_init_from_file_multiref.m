% Function to create new projection initialization file from results of a
% previous fast fast_multi_ref algorithm run
% Usage:   
%           proj_init_from_file_multiref(structurefile,initfile,writefile)
% Inputs:
%           structurefile = File name of results from previous run
%           initfile = File name of projection initialization used in
%           previous run
%           writefile = File name for new projection initialization
%           usestructure = Flag (0/1) to use structure variable
%           imvoxelres = (Optional, for lpf) Voxel size of volume 
%           lpf = (Optional) Low-pass filter cutoff in Angstrom
% Outputs:
%           (none)

% Created by Nicha C. Dvornek, 03/2015
% Last modified 06/2015


function proj_init_from_file_multiref(structurefile,initfile,writefile,usestructure,imvoxelres,lpf)

if ~strcmp(structurefile(end-3:end),'.mat')
    structurefile = [structurefile '.mat'];
end

if usestructure == 0
    load(structurefile, 'structure_final');
    structure = structure_final;
    if nargin == 6
        for s = 1:size(structure,4)
            structure(:,:,:,s) = lpf_at_res(structure_final(:,:,:,s),imvoxelres,lpf);
        end
    end
else
    load(structurefile,'structure');
end
load(initfile, 'maskim','mask','coord_axes','data_axes','ctfs','structmask')

if (~exist('structmask','var'))
    structmask = ones(size(structure_final));
end
structure = structure .* structmask;

numprojstruct = size(coord_axes,2);
for s = 1:size(structure,4)
    disp(['Projecting structure ' num2str(s)]);
    fproj_struct(:,:,(s-1)*numprojstruct+1:s*numprojstruct) = project_in_all_directions_w_ctf_par(structure(:,:,:,s),mask,coord_axes,ctfs);
end
numproj = size(fproj_struct,3);
proj_struct = zeros(size(fproj_struct));
for j = 1:numproj 
    proj_struct(:,:,j) = real(ifft2(ifftshift(fproj_struct(:,:,j)))).*maskim;
end
clear fproj_struct

disp('Saving')
save(writefile,'-v7.3','maskim','ctfs','data_axes','coord_axes','mask','structure','proj_struct','structmask');
