% script to test post-processing part of cryo3d

% given: aligned_ims, im_projinds, mu, mu_projinds, defocus
% pixA

%load('defocus.mat');
%load('fbm_aligned_ims.mat');
%im_projinds = projinds;

addpath('\cryoEm');
addpath('\cryoEm\EMBase');
addpath('cryoEm\EMIO');
addpath('cryoEm\EMSpec');
addpath('cryoEm\GriddingLib');

lambda = EWavelength(300);
pixA = 1.32; % is it the same as vos_size?

% simulate mu
mx = max(im_projinds);
mu = zeros(size(aligned_ims,1), size(aligned_ims,2), mx);
mu_projinds = zeros(mx,1);
for j=1:mx
    k=find(im_projinds==j);
    k1=k(1);
    mu(:,:,j)=aligned_ims(:,:,k1);
    mu_projinds(j) = k1;
end

I_fft_ctf = ctf_filtered_img(aligned_ims, defocus, pixA, lambda);
map = projection_map(im_projinds, mu_projinds);
mu_fft = get_ffts(mu);
norms = get_norms(I_fft_ctf, mu_fft, map);