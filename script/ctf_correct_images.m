% Script to ctf correct an image stack

% Code originally from Fred Sigworth
% Modified by Nicha C. Dvornek, 06/2015

% Parameter inputs, example:
%stackfile = 'C:\Users\Nicha\Documents\Data\Frank_Data\Rotated70swithEFGparticle.mrcs';
%paramfile = 'C:\Users\Nicha\Documents\Data\Frank_Data\Rotated70swithEFGparticle.star';
%B = 30; B-factor in units of A^2
%maxnim = 2000; (Optional) Number of images to read in for noise/signal spectra estimation. If stack is large, spectra estimation is slow. Default maxnim = half the stack size

function passed = ctf_correct_images(stackfile, paramfile, B, maxnim)

passed = 0;

%% Read in parameters
data = read_star_data_for_labels(paramfile,{'_rlnDetectorPixelSize','_rlnMagnification','_rlnVoltage','_rlnAmplitudeContrast'},1);
pixA = data{1} / data{2} * 10000; % pixel size in angstroms: rlnDetectorPixelSize / _rlnMagnification * 10000
volt = data{3}; % _rlnVoltage
ampContrast = data{4}; % _rlnAmplitudeContrast

% Get defocus values
disp('Read parameter file for defocus values');
params = read_star_data_for_labels(paramfile,{'_rlnDefocusU','_rlnDefocusV'});
d = mean(cell2mat(params),2) ./ 10000; % Average defocus values in um

%% Estimate noise and signal spectra from the stack
disp('Read data');
[~, s] = ReadMRC(stackfile,1,-1);
nim = s.nz;
% in case whole stack is very large, only read in maxnim images to reduce
% processing time for spectra estimation
if nargin < 4
    maxnim = floor(nim/2); % Default set to half the stack
end
if maxnim > nim
    disp('Warning: maxnim > number of images in the stack. Using all images for spectra estimation');
    maxnim = nim;
end
skip = floor(nim/maxnim);
stk = zeros(s.nx,s.ny,maxnim,'single');
for i = 1:maxnim
    stk(:,:,i) = ReadMRC(stackfile,i*skip,1);
end

n=size(stk,1);  % Had better be even, or my functions will have problems

% Estimate noise and signal spectra
disp('Estimate noise spectra');
msk=fuzzymask(n,2,0.45*n,.05*n);  % soft mask to separate rim from center
% The following call to RimSpectrum will take a long time to compute if the stack is large.
noiseSpec1d=RimSpectrum(stk,1-msk);  % get the average 1d spectrum using the rim of the imageswide on the sides)

disp('Estimate signal spectra');
noiseSpec2d=ifftshift(ToRect(noiseSpec1d,n));  % Polar to rect: in this case, makes it circularly symmetric
sigSpec1d=RadialPowerSpectrum(stk);
sigSpec1d(round(0.4*n):end)=mean(sigSpec1d(round(0.3*n):round(0.4*n)));  % force it to be constant up to Nyquist
sigSpec2d=ifftshift(ToRect(sigSpec1d,1));

k=noiseSpec2d./sigSpec2d;  % Or whatever, your Wiener 'constant'
   
%% Do Wiener filtering
disp('Do Wiener filtering');
clear stk;
filtImgs=zeros(n,n,nim,'single');
lambda=EWavelength(volt);
for i=1:nim
    c=ifftshift(CTF(n,pixA,lambda,d(i),0,B,ampContrast));  % This gives zero frequency at (1,1)
    img = ReadMRC(stackfile,i,1);
    filtImgs(:,:,i)=real(ifftn(fftn(img).*c./(k+c.^2)));
end;

writeMRC(filtImgs,s.pixA,[stackfile(1:strfind(stackfile,'.')-1) '_ctf_corrected.mrcs']);
passed = 1;