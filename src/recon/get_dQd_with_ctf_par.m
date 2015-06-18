% Function to be used in conjugate gradient reconstruction

% Original code by Hemant Tagare
% Modified for parallel processing by Nicha C. Dvornek, 04/2014

function dqd=get_dQd_with_ctf_par(d,mask_r,data_axes,l_norm,l_smooth,ctfs)

n_proj=size(data_axes,2);
dqd=0;

%if matlabpool('size') > 0
if (~isempty(gcp('nocreate')))
    if n_proj == size(ctfs,3)
        parfor i=1:n_proj
            img=project_volume_along_n_with_ctf_par(d,mask_r,data_axes(:,i),ctfs(:,:,i));
            img=real(ifft2(ifftshift(img)));
            dqd=dqd+img(:)'*img(:);
        end
    else
        n_ctfs = size(ctfs,3);
        parfor i=1:n_proj
            c = mod(i-1,n_ctfs) + 1;
            img=project_volume_along_n_with_ctf_par(d,mask_r,data_axes(:,i),ctfs(:,:,c));
            img=real(ifft2(ifftshift(img)));
            dqd=dqd+img(:)'*img(:);
        end
    end
else
    if n_proj == size(ctfs,3)
        for i=1:n_proj
            img=project_volume_along_n_with_ctf_par(d,mask_r,data_axes(:,i),ctfs(:,:,i));
            img=real(ifft2(ifftshift(img)));
            dqd=dqd+img(:)'*img(:);
        end
    else
        n_ctfs = size(ctfs,3);
        for i=1:n_proj
            c = mod(i-1,n_ctfs) + 1;
            img=project_volume_along_n_with_ctf_par(d,mask_r,data_axes(:,i),ctfs(:,:,c));
            img=real(ifft2(ifftshift(img)));
            dqd=dqd+img(:)'*img(:);
        end
    end
end

if l_norm > 0
    dqd = dqd+l_norm*(d(:)'*d(:));
end

if l_smooth > 0
    l_d=laplacian_3d(d);
    dqd=dqd+l_smooth*(l_d(:)'*l_d(:));
end
