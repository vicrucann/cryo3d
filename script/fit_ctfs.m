% function to cluster and fit ctfs to the data
% minimal input parameters, example:
% pathout = 'G:\workspace\';
% paramfile = 'G:\db-frank\Rotated70swithEFGparticle.star';
% stackfile = 'G:\db-frank\stack_ds4.mrc';
% optional parameters, example:
% num_clusters = 5; Number of defocus classes
% ds = 2; Factor by which to downsample
% fitflag = 1; Flag (0/1) to run ctfit2 to get ctf parameters
% saveflag = 1; Flag (0/1) to save the ctf parameters and images in .mat

function ctffile = fit_ctfs(pathout, paramfile, stackfile, num_clusters, ds, fitflag, saveflag)

if (nargin < 7)
    saveflag = 1;
end
if (nargin < 6)
    fitflag = 0;
end
if (nargin < 5)
    ds = 1;
end
if (nargin < 4)
    num_clusters = 5;
end
if (nargin < 3)
    disp('ERROR: Not enough input parameters');
    return;
end

%%
addpath(fullfile(cd, '../src/preprocessing'));
addpath(fullfile(cd, '../src/mrc'));

%% Data Params (read from paramfile in .star format)
data = read_star_data_for_labels(paramfile,{'_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration'},1);
volt = data{1}; % _rlnVoltage, in kV
ampcont = data{2}; % _rlnAmplitudeContrast 
Cs = data{3}; % _rlnSphericalAberration, in mm 

%% Preprocessing Params
K = num_clusters; % num ctf clusters

%% CTF CLUSTERING

% Load defocus values (_rlnDefocusU and _rlnDefocusV) 
disp('Read parameter file for defocus values');
params = read_star_data_for_labels(paramfile,{'_rlnDefocusU','_rlnDefocusV'});
ctfdata = cell2mat(params) / 10000; % Defocus values in .star file are in angstrom; need values in um

% Cluster CTFs into K sets
disp('K-means clustering of defocus values');
[ctfinds, C, ~, D] = kmeans(ctfdata, K,'Distance','cityblock',...
                                'EmptyAction','singleton',...
                                'Replicates',10);

% Visualize clustering
f = figure;
hold on
colorsALP = [1 0 0; 0 0 1; 0 1 0;...
             1 0 1; 0 1 1; 1 1 0];
for i = 1:K
    plot(ctfdata(ctfinds==i,1),ctfdata(ctfinds==i,2),'.','color',colorsALP(mod(i-1,6)+1,:))
end
plot(C(:,1),C(:,2),'ks','LineWidth',2);
set(gca,'XMinorTick','on');
set(gca,'yMinorTick','on');
grid on;
xlabel('rlnDefocusU (um)');
ylabel('rlnDefocusV (um)');

% Visually select images close to clusters.
cluDist = 2.0e9;

% Calculate how many images you get
% (I'm actually not sure what this is about....from Alp's code)
total = 0;
for i = 1:K
    total = total + sum(ctfinds==i&D(:,i)<cluDist);
end
total
disp('DONE.\n');

%% GET CTF PARAMETERS

% Some initial things
ctfParams = cell(K,1);
[h, s] = ReadMRC(stackfile,1,-1);
Apix = s.pixA;
imgSx = s.nx;

if fitflag  % Fit CTF parameters to the images
    
    disp('Fit CTF parameters to the images');
    
    % Load image stack
    disp('Load the image stack');
    noisyims_raw = ReadMRC(stackfile);

    % Make each image zero mean
    numim = size(noisyims_raw,3);
    for i = 1:numim
        noisyims_raw(:,:,i) = noisyims_raw(:,:,i) - mean(reshape(noisyims_raw(:,:,i),imgSx^2,1));
    end

    % Split images into K CTF clusters
    images = cell(K,1);
    for i = 1:K
        images{i,1} = noisyims_raw(:,:,ctfinds==i&D(:,i)<cluDist);
    end
    clear noisyims_raw

    % Run Fred/Niko's CTFfit code on each cluster to be sure of CTF parameters

    % Set search parameters
    Pa.lambda   = EWavelength(volt);
    Pa.defocus  = 0:0.5:10;
    Pa.deltadef = -2:.2:2;
    Pa.theta    = 0;
    Pa.alpha    = ampcont;
    Pa.Cs       = Cs;
    Pa.B        = 180;

    % For each CTF cluster
    for i = 1:K
        fprintf('CTF CLUSTER: %d\n',i);

        % Get max 4k images from stack... more makes CTFfit run really slowly
        singleSet = images{i,1};
        randomIdx = randperm(size(singleSet,3));
        numForCTFfit = min(size(singleSet,3),4000);
        montageSx = ceil(sqrt(numForCTFfit));
        figure; montageHandle = montage(reshape(singleSet(:,:,randomIdx(1:numForCTFfit)),[imgSx, imgSx, 1, numForCTFfit]),'Size',[montageSx montageSx],'DisplayRange',[]);
        singleMontage = getimage(montageHandle);
        close
        clear singleSet

        % Run CTFfit
        fprintf('  Running ctfit2...');
        [P, ~] = ctfit2(singleMontage,Pa,Apix,8,55,0);
        clear singleMontage
        fprintf('DONE.\n');
        P.B   = 0;

        % Save CTFfit outputs
        ctfParams{i,1} = P;
    end
    
else        % Use K-means cluster centers as defocus value, ignore astigmatism
    
    disp('Use K-means results to set ctf parameters');
    
    P.lambda   = EWavelength(volt);
    P.alpha    = ampcont;
    P.Cs       = Cs;
    P.B        = 0;
    defocus_means = mean(C,2);
    
    for i = 1:K
        P.defocus = defocus_means(i);
        ctfParams{i,1} = P;
    end
    
end

%% WRITE FILES OUT

if saveflag
    if ds > 1
        imgSx = floor(imgSx/ds);
        if mod(imgSx,2) ~= 0
            imgSx = imgSx + 1;
        end
        Apix = Apix* s.nx / imgSx;
    end
    ctfs = zeros(imgSx,imgSx,K,'single');
    for i = 1:K
        ctfs(:,:,i) = CTF(imgSx,Apix,ctfParams{i,1});
    end
    %savename = [stackfile(1:strfind(stackfile,'.')-1) '_' num2str(K) 'ctfs' ];
    f_type = strfind(stackfile,'.mrc');
    f_path = strfind(stackfile, '\');
    if (isempty(f_path))
        f_path = strfind(stackfile, '/');
    end
    savename = [stackfile(max(f_path)+1:f_type-1) '_' num2str(K) 'ctfs'];
    if ds > 1
        savename = [savename '_ds' num2str(ds)];
    end
    save([pathout savename],'ctfs','ctfParams','ctfinds');
end
ctffile = [pathout savename '.mat'];
