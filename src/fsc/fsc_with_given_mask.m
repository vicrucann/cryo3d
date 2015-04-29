% Function to to calculate FSC between two volumes
% Usage:   
%           [fsc,f_ind,mask3d] = fsc_with_mask(m1,m2,voxel_size,tick_spacing,pixfromedge,maskwidth)
% Inputs:
%           m1 = Volume 1
%           m2 = Volume 2
%           voxel_size = Physical size of voxel in Angstrom
%           tick_spacing = Tick spacing along x-axis
%           mask3d = The mask in Fourier space
% Outputs:
%           fsc = Fourier Shell Correlation
%           f_ind = Frequency axis
%           mask3d = Mask used to calculate FSC

% Modified from Alp's code by Nicha C. Dvornek, 07/2014

function [fsc,f_ind] = fsc_with_given_mask(m1,m2,voxel_size,tick_spacing,mask3d)


n       = size(m1,1);
[x y z] = ndgrid(-n/2:n/2-1);  % zero at element n/2+1.
R       = sqrt(x.^2 + y.^2 + z.^2);
eps     = 1e-4;

f1  = fftshift(fftn(m1.*mask3d));
f2  = fftshift(fftn(m2.*mask3d));

fsc = zeros(1,n/2-1);

for i=1:n/2-1
    ring   = find(R<0.5+i+eps & R >= i-0.5+eps);
    num    = real(sum(f1(ring).*conj(f2(ring))));
    den    = norm(f1(ring))*norm(f2(ring));
    fsc(i) = num/den;
end;

Fs = 1/voxel_size;

figure
plot(fsc)
ylim([0 1])
set(gca,'YTick',[0,0.143,0.5,1])
set(gca,'XTick',1:tick_spacing:n/2-1)
xlim([1 n/2])
f_ind = Fs/2*linspace(0,1,n/2+1);
f_ind = f_ind(2:end-1);
f_ind_s = 1./interp1(f_ind,1:tick_spacing:length(f_ind));
num = length(f_ind_s);
strOneOver = repmat('1/',num,1);
endBrace = repmat('',num,1);
f_ind_s = horzcat(strOneOver,strjust(num2str(f_ind_s','%2.1f'),'left'),endBrace);
set(gca,'XTickLabel',f_ind_s)
ylabel('FSC')
xlabel('Resolution [A^{-1}]')
grid on