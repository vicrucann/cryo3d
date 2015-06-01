% Modification of H. Tagare's projection code to work with parfor.
% Requires that a matlabpool has already been opened - function will NOT
% automatically start a matlabpool

function [proj, data_axes] = project_in_all_directions_w_ctf_par(x,mask,coord_axes,ctfs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This routine projects the 3D data x in along normals
%  x: 3D data. Has to be N X N X N
%  mask: 0-1 valued matrix for masking data, N X N X N 
%  ctf_params: M X 2 matrix. Each row is a set of ctf params
%  normals: K X 3 matrix. Each row is a normal vector
%  disp_flag: when 1 the routine displays the projections as images
%  fig_num: Figure number for displaying images. Can be null if
%                   disp_flag==0
%
%    The routine requires environmental variables CG_UTIL 
%   Code history:
%       Written by H. Tagare  30/06/10.
%       Modified by H. Tagare 9/14/10 to include normal and ctf arrays
%       Modified by Nicha C. Dvornek 04/2014 to for parallel processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set the utilities directory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Get the utilities directory
% util_dir=getenv('CG_UTIL');
% if isempty(util_dir)
%     disp('project_in_all_directions: CG_UTIL environment variable not set');
%     return
% else
%     path(util_dir,path);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin computation by initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize projection array
n_normals=size(coord_axes,2);
if size(ctfs,3) > n_normals
    n_ctfs=size(ctfs,3)/n_normals;
else
    n_ctfs = size(ctfs,3);
end
proj=complex(zeros(size(x,1),size(x,2),n_normals*n_ctfs), ...
             zeros(size(x,1),size(x,2),n_normals*n_ctfs));
               
%if matlabpool('size') > 0 % Matlab2013
if (~isempty(gcp('nocreate'))) % Matlab2014, to be compatible with future matlab versions
    % Project using parfor
    disp('Using matlab pool');
    numproj = n_normals*n_ctfs;
    data_axes = zeros(size(coord_axes,1),numproj);
    for c = 1:n_ctfs
        data_axes(:,c:n_ctfs:end) = coord_axes;
    end
    mask_r = mask.r;
    if n_ctfs == size(ctfs,3)
        parfor i = 1:numproj
            c = mod(i-1,n_ctfs) + 1;
            proj(:,:,i) = project_volume_along_n_with_ctf_par(x,mask_r,data_axes(:,i),ctfs(:,:,c));
        end
    else
        parfor i = 1:numproj
            proj(:,:,i) = project_volume_along_n_with_ctf_par(x,mask_r,data_axes(:,i),ctfs(:,:,i));
        end
    end
else
    %Project along the normals serially
    if n_ctfs == size(ctfs,3)
        img_i=1;
        for i=1:n_normals,
            for j=1:n_ctfs,
                proj(:,:,img_i)=project_volume_along_n_with_ctf_nd(x,mask,coord_axes(:,i),ctfs(:,:,j));
                img_i=img_i+1;
            end
        end
    else
        img_i=1;
        for i=1:n_normals,
            for j=1:n_ctfs,
                proj(:,:,img_i)=project_volume_along_n_with_ctf_nd(x,mask,coord_axes(:,i),ctfs(:,:,img_i));
                img_i=img_i+1;
            end
        end
    end
    
end




