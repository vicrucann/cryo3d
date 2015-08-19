% Script to test post-processing part of cryo3d
%% Additional libraries

addpath('cryoEm');
addpath('cryoEm\EMBase');
addpath('cryoEm\EMIO');
addpath('cryoEm\EMSpec');
addpath('cryoEm\GriddingLib');
addpath('other');
addpath('data');

%% The input variables

% aligned_ims, mu
% im_projinds, mu_projinds
% defocus
% pixA, lambda

load('defocus.mat');
%load('fbm_aligned_ims.mat');
load('eigStruct.mat');
load('coord_axes_6.mat');
%im_projinds = projinds;

lambda = EWavelength(300);
pixA = 1.32; % obtained by using ReadMRC

%% Calculate mu

%=======
% simulate mu
% mx = max(im_projinds);
% mu = zeros(size(aligned_ims,1), size(aligned_ims,2), mx);
% mu_projinds = zeros(mx,1);
% for j=1:mx
%     k=find(im_projinds==j);
%     k1=k(1);
%     mu(:,:,j)=aligned_ims(:,:,k1);
%     mu_projinds(j) = k1;
% end
%=======

eig1 = eigStruct.eigs3d(:,:,:,1);
eig2 = eigStruct.eigs3d(:,:,:,2);
eigVal1 = norm(eig1(:));
eigVal2 = norm(eig2());
eigVec1 = eigStruct.eigs3d(:,:,:,1)./eigVal1;
eigVec2 = eigStruct.eigs3d(:,:,:,2)./eigVal2;
n1 = normrnd(0, eigVal1);
n2 = normrnd(0, eigVal2);
M = size(eigStruct.meanVol, 3);
rMax = floor(M/2-2);
num_projdir = size(eigStruct.eigsProj, 3);

mu_mean = mex_forward_project(eigStruct.meanVol, M, coord_axes, num_projdir, rMax);
mu_vec1 = n1 * mex_forward_project(eigVec1, size(eigVec1), coord_axes, size(eigStruct.eigsProj, 3), rMax);
mu_vec2 = n2 * mex_forward_project(eigVec2, size(eigVec2), coord_axes, size(eigStruct.eigsProj, 3), rMax);

mu = mu_mean + mu_vec1 + mu_vec2;

%% Get the norms

I_fft_ctf = ctf_filtered_img(aligned_ims, defocus, pixA, lambda);
map = projection_map(im_projinds, mu_projinds);
mu_fft = get_ffts(mu);
norms = get_norms(I_fft_ctf, mu_fft, map);