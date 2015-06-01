function x=back_project_image_along_n_with_ctf_par(img,mask_r,n,ctf)

%apply CTF
%img=fftshift(fft2(img));
img=real(ifft2(ifftshift(img.*ctf)));
x=mex_back_project(img,n,mask_r);
