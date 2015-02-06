function g=get_grad_with_ctf_par(x,mask_r,proj,data_axes,l_norm,l_smooth,ctfs)

n_proj=size(proj,3);
g=zeros(size(x));

%if matlabpool('size') > 0
if (~isempty(gcp('nocreate')))
    if n_proj == size(ctfs,3)
        parfor i=1:n_proj,
                g = g + get_grad_component(x,mask_r,proj(:,:,i),...
                                           data_axes(:,i),...
                                           ctfs(:,:,i));
        end
    else
        n_ctfs = size(ctfs,3);
        parfor i=1:n_proj,
                c = mod(i-1,n_ctfs) + 1;
                g = g + get_grad_component(x,mask_r,proj(:,:,i),...
                                           data_axes(:,i),...
                                           ctfs(:,:,c));
        end
    end
else
    if n_proj == size(ctfs,3)
        for i=1:n_proj,
            g=g+get_grad_component(x,mask_r,proj(:,:,i),...
                                          data_axes(:,i),...
                                          ctfs(:,:,i));
        end
    else
        n_ctfs = size(ctfs,3);
        for i=1:n_proj,
            c = mod(i-1,n_ctfs) + 1;
            g=g+get_grad_component(x,mask_r,proj(:,:,i),...
                                          data_axes(:,i),...
                                          ctfs(:,:,c));
        end
    end
end

if l_norm > 0
    g=g+l_norm*x;
end
if l_smooth > 0
    g=g+l_smooth*laplacian_3d(x);
end


function g_comp=get_grad_component(x,mask_r,img,cur_axes,ctf)

img_r=project_volume_along_n_with_ctf_par(x,mask_r,cur_axes,ctf);
g_comp=back_project_image_along_n_with_ctf_par(img_r-img,mask_r,cur_axes,ctf);
