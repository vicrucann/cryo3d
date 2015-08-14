% Script to test post-processing part of cryo3d
%% Additional libraries

addpath('\cryoEm');
addpath('\cryoEm\EMBase');
addpath('cryoEm\EMIO');
addpath('cryoEm\EMSpec');
addpath('cryoEm\GriddingLib');

%% The input variables

% aligned_ims, mu
% im_projinds, mu_projinds
% defocus
% pixA, lambda

%load('defocus.mat');
%load('fbm_aligned_ims.mat');
%load('eigStruct.mat');
%load('coord_axes_12.mat');
%im_projinds = projinds;

lambda = EWavelength(300);
pixA = 1.32; % obtained by using ReadMRC

%% Calculate mu

% % simulate mu
% mx = max(im_projinds);
% mu = zeros(size(aligned_ims,1), size(aligned_ims,2), mx);
% mu_projinds = zeros(mx,1);
% for j=1:mx
%     k=find(im_projinds==j);
%     k1=k(1);
%     mu(:,:,j)=aligned_ims(:,:,k1);
%     mu_projinds(j) = k1;
% end

eigVal1 = norm(eigStruct.eigs3d(:,:,:,1));
eigVal2 = norm(eigStruct.eigs3d(:,:,:,2));
eigVec1 = eigStruct.eigs3d(:,:,:,1)./eigVal1;
eigVec2 = eigStruct.eigs3d(:,:,:,2)./eigVal2;
n1=normrnd(0, eigVal1);
n2=normrnd(0, eigVal2);
M=size(eigStruct.meanVol, 3);
rMax=floor(M/2-2);
%mu = zeros(M,M,size(eigStruct.eigsProj, 3));
%for k=1:size(eigStruct.eigsProj, 3)
 mu = mex_forward_project(eigStruct.meanVol, size(eigStruct.meanVol), coordAxes, size(eigStruct.eigsProj, 3), rMax)...
        + n1 * mex_forward_project(eigVec1, size(eigVec1), coordAxes, size(eigStruct.eigsProj, 3), rMax)...
        + n2 * mex_forward_project(eigVec2, size(eigVec2), coordAxes, size(eigStruct.eigsProj, 3), rMax);
%end

%% Get the norms

I_fft_ctf = ctf_filtered_img(aligned_ims, defocus, pixA, lambda);
map = projection_map(im_projinds, mu_projinds);
mu_fft = get_ffts(mu);
norms = get_norms(I_fft_ctf, mu_fft, map);