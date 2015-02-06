% Function to create new projection initialization file from results of a
% previous fast em-algorithm run
% Usage:   
%           proj_init_from_file(structurefile,initfile,writefile)
% Inputs:
%           structurefile = File name of results from previous run
%           initfile = File name of projection initialization used in
%           previous run
%           writefile = File name for new projection initialization
% Outputs:
%           (none)

% Created by Nicha C. Dvornek, 08/2013


function proj_init_from_file(structurefile,initfile,writefile,s,imvoxelres,lpf)

if s == 0
    load(structurefile, 'structure_final');
    structure = structure_final;
    if nargin == 6
        structure = lpf_at_res(structure_final,imvoxelres,lpf);
    end
else
    load(structurefile,'structure');
end
load(initfile, 'maskim','mask','coord_axes','data_axes','ctfs','structmask')

if (~exist('structmask','var'))
    structmask = ones(size(structure_final));
end
structure = structure .* structmask;

disp('Projecting');

fproj_struct = project_in_all_directions_w_ctf_par(structure,mask,coord_axes,ctfs);
numproj = size(fproj_struct,3);
proj_struct = zeros(size(fproj_struct));
data = zeros(numproj,numel(maskim));
for j = 1:numproj
    temp = real(ifft2(ifftshift(fproj_struct(:,:,j)))).*maskim;
    proj_struct(:,:,j) = temp;
    data(j,:) = temp(:);
end
clear fproj_struct temp

disp('PCA')
[coeffproj,scoreproj,latentproj] = princomp(data);

disp('Saving')
save(writefile,'-v7.3','maskim','coeffproj','scoreproj','latentproj',...
    'ctfs','data_axes','coord_axes','mask','structure','proj_struct',...
    'structmask');