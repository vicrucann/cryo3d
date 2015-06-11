% Function for pre-processing the particle images, including downsampling,
% phase flipping, prewhitening, and normalizing the intensities to have
% 0 mean, 1 std dev

% Nicha C. Dvornek 02/2015

% Input parameters, example:
% 1. pathout - user's working directory, all the intermediates and results
% will be saved there, pathout = 'G:\workspace\';
% 2. stackfile - *.mrc file, stackfile = 'G:\db-frank\stack_ds4.mrc';
% 3. ctffile = *.mat file with clustering info, ctffile = 'G:\db-frank\stack_ds4_5ctfs.mat';
% 4. downsample = factor by which to downsample images, downsample = 2;
% 5. normflag = flag (0/1) to normalize image intensities to have 0 mean, 1 std dev, normflag = 1;
% 6. pfflag = flag (0/1) to phase flip the images, pfflag = 1;
% 7. pwflag = flag (0/1) to prewhiten images (input 7 and 8 required if pwflag = 1), pwflag = 1;
% 8. npsfile = normally NPS.txt file, npsfile = 'G:\db-frank\NPS.txt';
% 9. Apix - pixel size of micrograph for NPS data in Angstroms, Apix = 1.045;

function outfile = preprocess_images(pathout, stackfile, ctffile, downsample, normflag, pfflag, pwflag, npsfile, Apix)

addpath(fullfile(cd, '../src/preprocessing'));
addpath(fullfile(cd, '../src/mrc'));

% Data and Params
if (nargin < 7)
    pwflag = 0; % Flag for whether or not to prewhiten the images
end
if (nargin < 6)
    pfflag = 0; % Flag for whether or not to phase flip the images
end
if (nargin < 5)
    normflag = 1; % Flag for whether to normalize image intensities
end
if (nargin < 4)
    downsample = 1; % Factor by which to downsample
end
if (nargin < 3)
    disp('ERROR: Not enough input parameters');
    return;
end
 

%% Load things

disp('Loading CTF info, NPS info, and image stack');

% Load NPS info
% THIS IS HOW I READ IT IN FOR THE EXAMPLE FILE - NOT SURE IF THERE IS A
% STANDARD FILE TYPE THAT IS TO BE READ IN
if pwflag
    f = importdata(npsfile);
    dqe_freq = f.data(:,1)/(2*Apix);
    nps = f.data(:,3);
    clear f
end

% Load image stack
[noisyims, h] = ReadMRC(stackfile);
Apixstack = h.pixA;
imgSx = size(noisyims,1);

% Make each image zero mean
numim = size(noisyims,3);
for i = 1:numim
    noisyims(:,:,i) = noisyims(:,:,i) - mean(reshape(noisyims(:,:,i),imgSx^2,1));
end

% Load CTF things
if isempty(ctffile)
    ctfinds = ones(numim,1);
    K = 1;
else
    load(ctffile,'ctfinds','ctfParams');
    K = size(ctfParams,1);
end

% Downsample each image
if downsample > 1
    disp(['Downsample images by ' num2str(downsample)]);
    imgSx = floor(imgSx/downsample);
    if mod(imgSx,2) ~= 0
        imgSx = imgSx + 1;
    end
    Apixstack = Apixstack * size(noisyims,1) / imgSx;
    tempnoisyims = zeros(imgSx,imgSx,size(noisyims,3),'single');
    for i = 1:size(noisyims,3)
        temp = DownsampleGeneral(noisyims(:,:,i),imgSx);
        tempnoisyims(:,:,i) = temp - mean(temp(:));
    end
    noisyims = tempnoisyims;
    clear tempnoisyims temp
end


%% IMAGE PROCESSING

% For each CTF cluster
for i = 1:K
    
    fprintf('CTF CLUSTER: %d\n',i);

    % Get all images from stack
    singleSet = noisyims(:,:,ctfinds == i);
    numinset = size(singleSet,3);
    montageSx = ceil(sqrt(size(singleSet,3)));
    montageHandle = montage(reshape(singleSet,[imgSx, imgSx, 1, numinset]),'Size',[montageSx montageSx],'DisplayRange',[]);
    clear singleSet
    singleMontage = getimage(montageHandle);
    close
    
    % Do phase flipping if flag is set
    if pfflag 
        fprintf('  Phase Flipping Images...');
        
        % Create CTF image
        ctfImg = CTF(montageSx*imgSx,Apixstack,ctfParams{i,1});

        % Make into phase flipping image
        ctfImg(ctfImg<0) = -1; 
        ctfImg(ctfImg>0) = 1;
        ctfImg           = ctfImg.*-1;

        % Phase flip 
        singleMontage = real(ifft2(ifftshift(ctfImg.*fftshift(fft2(singleMontage)))));
        clear ctfImg
        fprintf('DONE.\n');       
    end
    
    % Do prewhitening if flag is set
    if pwflag       
        % Pre-whitening each "micrograph" bc whole set is too much memory
        fprintf('  Create prewhitening filter....\n');
        sm_size = size(singleMontage,1);
        f_ind = get_freq_axis(Apixstack,sm_size);
        pfilt = interp1(dqe_freq,1./sqrt(nps),1./f_ind);
        [x, y] = ndgrid(single(-sm_size/2:sm_size/2-1)); % zero at element n/2+1.
        R = sqrt(x.^2 + y.^2);
        clear x y
        pfilt2d = interp1(0:sm_size/2-1,pfilt,R(:));
        clear R
        pfilt2d(isnan(pfilt2d)) = pfilt(end);
        pfilt2d = reshape(pfilt2d,[sm_size,sm_size]);
        fprintf('  Prewhiten the images...\n');
        singleMontage = real(ifft2(ifftshift(pfilt2d.*fftshift(fft2(singleMontage)))));
        clear pfilt2d      
    end
    
    % Rearrange montage image back into stack and normalize intensities
    if normflag
        fprintf('  Normalize image intensities...\n');
    end
    tempset = zeros(imgSx,imgSx,numinset,'single');
    idx = 1;
    for j = 1 : imgSx : size(singleMontage,1)
        for k = 1 : imgSx : size(singleMontage,1)
            temp = singleMontage(j:j+(imgSx-1),k:k+(imgSx-1));
            if normflag
                tempset(:,:,idx) = (temp - mean(temp(:))) ./ std(temp(:));
            else
                tempset(:,:,idx) = temp;
            end
            idx = idx + 1;
        end
    end
    noisyims(:,:,ctfinds == i) = tempset(:,:,1:numinset);
    clear singleMontage tempset

    fprintf('DONE.\n');
    
end


%% SAVE THE PROCESSED IMAGES

f_type = strfind(stackfile,'.mrc');
f_path = strfind(stackfile, '\');
if (isempty(f_path))
    f_path = strfind(stackfile, '/');
end
if isempty(f_path)
    f_path = 0;
end
savefile = stackfile(max(f_path)+1:f_type-1); %stackfile(1:f-1);
if downsample > 1
    savefile = [savefile '_ds' num2str(downsample)];
end
if pfflag
    savefile = [savefile '_pf'];
end
if pwflag
    savefile = [savefile '_pw'];
end
if normflag
    savefile = [savefile '_norm'];
end
savefile = [savefile '.mrcs'];
writeMRC(noisyims,Apixstack,[pathout savefile]);
outfile = [pathout savefile];