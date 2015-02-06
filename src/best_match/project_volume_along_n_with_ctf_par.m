function img=project_volume_along_n_with_ctf_par(x,mask_r,n,ctf)

% %Get the insertion directions
% n=n/sqrt(sum(n.^2));
% [tform,tforminv,R,TDIMS_A,TDIMS_B,TSIZE_B,TMAP_B,F]=get_tform_param_from_n(size(x),n);
% %Rotate the array and sum
% tmp = tformarray(x, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
% tmp=tmp.*mask;
% img=squeeze(sum(tmp,3));


%Get the insertion directions
%T=get_tform_coord_from_n(n);
img=mex_forward_project(x,n,mask_r);
img=fftshift(fft2(img)).*ctf;
%img=real(ifft2(ifftshift(img)));
%img=apply_ctf(img,ctf_param(1),ctf_param(2));