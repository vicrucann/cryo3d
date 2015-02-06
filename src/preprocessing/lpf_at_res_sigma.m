% Function to low-pass filter volume at given Gaussian width
% Usage:   
%           outvol = lpf_at_res_sigma(vol,voxel_size,resolution,sigma)
% Inputs:
%           vol = Volume to low pass
%           voxel_size = Size of voxel (in Angstroms)
%           resolution = Resolution at which to low-pass filter the volume
%           sigma = Standard deviation of Gaussian used to smooth data
% Outputs:
%           outvol = Low-pass filtered volume

% Original code by Alp, modified by Nicha C. Dvornek, 08/2013

function outvol = lpf_at_res_sigma(vol,voxel_size,resolution,sigma)

n = size(vol,1);

if(ndims(vol)==2)
  [x y] =    ndgrid(-n/2:n/2-1);  % zero at element n/2+1.
  R = sqrt(x.^2 + y.^2);
elseif (ndims(vol==3))
    [x y z] = ndgrid(-n/2:n/2-1);  % zero at element n/2+1.
    R       = sqrt(x.^2 + y.^2 + z.^2);
end

Fs = (1/voxel_size)/1e-10;
num = size(1:n/2,2);
f_ind = 1./(Fs/2*linspace(0,1,num));
f_ind = f_ind./1e-10;

[~, closestIndex] = min(abs(f_ind - resolution));

R = (R<closestIndex-1);

if(ndims(vol)==2)
    G = fspecial('gaussian',7,sigma);
    R = filter2(G,R);
elseif(ndims(vol)==3)
    R = smooth3(R,'gaussian',[7 7 7],sigma);
end

outvol = real(ifftn(ifftshift(R.*fftshift(fftn(vol)))));
