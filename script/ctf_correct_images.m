% Script to ctf correct an image stack

% Code originally from Fred Sigworth
% Modified by Nicha C. Dvornek, 02/2015


%% Parameter inputs
stackfile = 'C:\Users\Nicha\Documents\Data\Frank_Data\Rotated70swithEFGparticle.mrcs';
paramfile = 'C:\Users\Nicha\Documents\Data\Frank_Data\Rotated70swithEFGparticle.star';

% WOULD BE NICE IF WE COULD READ THIS IN FROM THE .STAR FILE
pixA=6.35/60765.550200*10000;  % pixel size in angstroms: rlnDetectorPixelSize / _rlnMagnification * 10000
volt = 300; % _rlnVoltage
ampContrast=.1628;  % _rlnAmplitudeContrast

B=30;  % not sure what to use for B-factor, maybe 30 in units of A^2

%% Read in parameters
% Get defocus values
fid = fopen(paramfile);
% THIS ONLY WORKS FOR THE GIVEN .STAR FILE
params = textscan(fid,'%f %f %f %f %f %f %f %f %f %s %f %f %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f','headerlines',31);
fclose(fid);
d = cell2mat(params(2:3));
clear params
d = mean(d,2) ./ 10000; % in um

%% Estimate noise and signal spectra from the stack
disp('Read data');
[h, s] = ReadMRC(stackfile,1,-1);
nim = s.nz;
% in case whole stack is very large, only read in maxnim images to reduce
% processing time for spectra estimation
maxnim = 2000;
if nim > maxnim
    skip = floor(nim/maxnim);
end
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

writeMRC(filtImgs,s.pixA,[stackfile(1:strfind(imfile,'.')-1) '_ctf_corrected.mrcs']);
