% Script to run SubspaceEM algorithm

% Created by Nicha C. Dvornek, 09/2013
% Last modified 06/2015


function passed = subspaceEM(pathout, configfile)

clear; reset(gpuDevice);

% Add paths
addpath(fullfile(cd, '../src/recon'));
addpath(fullfile(cd, '../src/mrc'));

totaltime = tic;

% Set up - user parameters
[imfile,imreconfile,ctffile,initprojfile,maxmem,numthreads,dispflag,substep,reconhalf,reconstartind,normprojint,numruns,maxnumiter,rotstart,rotstep,rotend,transmax,transdelta,transwidth,convtol,t] = read_config_file_subspaceEM(configfile);

% Reconstruction parameters
l_norm = 0.1;
l_smooth = -10;
iter_lim = 10;
stop_lim = 0.03;

% Start pool
%if numthreads > 0 % matlab2013
%   if matlabpool('size') > 0
%        matlabpool close;
%    end
%    matlabpool('local',numthreads);
%end
if (~isempty (gcp('nocreate')) ) % matlab 2014, may not be needed
    delete(gcp('nocreate'));
end
parpool(numthreads);

%% SubspaceEM outer loop
for run = 1:numruns
       
    runtime = tic;
    disp(' ');
    disp(['Run ' num2str(run)]);
        
    if run > 1  
        reset(gpuDevice);
        writefile = [pathout '/' savename '_reinit.mat'];  
        proj_init_from_file([pathout '/' savename '.mat'],initprojfile,writefile,0);
        initprojfile = writefile;        
    end
  
    % Load things
    if run == 1 
        disp('Loading images and setting up image subspace'); tic;

        % Check if pca of images already exists
        pcafile = [imfile(1:strfind(imfile,'.')-1) '_pca.mat'];
        if ~exist(pcafile,'file');
            % Do PCA
            disp('PCA of images');
            data = single(ReadMRC(imfile));
            data = reshape(data,[size(data,1)*size(data,2),size(data,3)])';  
            tic; [coeffim, scoreim, latentim] = pca(data); toc;
            clear data
            save(pcafile,'-v7.3','coeffim','scoreim','latentim');
        end
        
        % Determine subspace size of images
        if ~exist('latentim','var');
            load(pcafile,'latentim');
        end
        latent_der = latentim(2:end) - latentim(1:end-1);
        latent_der2 = latent_der(2:end) - latent_der(1:end-1);
        latent_der2_avg = conv([1 1 1 ]./3,abs(latent_der2));
        numimcoeffs = find(abs(latent_der2_avg) < t,1) + 3;
        if isempty(numimcoeffs)
            disp('Warning: threshold too low - using full PCA basis');
            numimcoeffs = length(latentim);
        end
        clear latentim;
        disp(['Number of image basis elements: ' num2str(numimcoeffs)]);
        
        % Set up image subspace and coeffs
        if ~exist('coeffim','var')
            load(pcafile,'coeffim');
        end
        imbasis = coeffim(:,1:numimcoeffs-1); % Each col is numpixel-length basis elem
        clear coeffim
        if ~exist('scoreim','var')
            load(pcafile,'scoreim');
        end
        imcoeffs = [ones(size(scoreim,1),1), scoreim(:,1:numimcoeffs-1)]; % Each row is set of coeffs for an image
        clear scoreim;
        
        % Load images/ctf indices and put in mean vector into subspace
        noisyims = single(ReadMRC(imfile));
        if isempty(ctffile)
            ctfinds = ones(size(noisyims,3),1);
        else
            load(ctffile,'ctfinds');
        end
        meanim = mean(noisyims,3);
        imbasis = [meanim(:), imbasis];
        clear meanim;
        
        %%%% FOR TESTING %%%%%%%%%%%%%%%%%%%%
        if substep > 0
            disp(['DIVIDING TOTAL DATA BY ' num2str(substep)]);
            noisyims = noisyims(:,:,1:substep:end);
            ctfinds = ctfinds(1:substep:end);
            imcoeffs = imcoeffs(1:substep:end,:);
        end

        % Reconstruct half?
        if reconhalf
            disp(['Reconstruct half ' num2str(reconstartind) ' of the data']);
            noisyims = noisyims(:,:,reconstartind:2:end);
            ctfinds = ctfinds(reconstartind:2:end);
            imcoeffs = imcoeffs(reconstartind:2:end,:);
        end

        toc;
    end
    
    disp('Loading initial structure and projections'); tic;
    if run > 1
        clear proj_struct
        load(initprojfile,'coeffproj','scoreproj','latentproj','maskim','mask','coord_axes','data_axes','ctfs','structure','proj_struct','structmask');
    else    
        % Get initial model parameters      
        initprojfile = make_init_model(configfile,pathout);
        load(initprojfile,'maskim','mask','coord_axes','data_axes','ctfs','structure','proj_struct','structmask');
        
        % Rescale projection intensities to match images
        if normprojint
            disp('Rescaling projection intensities and initialize basis');
            proj_struct = rescale_proj_ints(proj_struct,noisyims);
        end
        
        disp('PCA of initial projections');
        data = reshape(proj_struct,[size(proj_struct,1)*size(proj_struct,2),size(proj_struct,3)])';
        pcatic = tic; [coeffproj, scoreproj, latentproj] = princomp(data); toc(pcatic);
        clear data;
    end
    
    % Set up projection subspace and coeffs
    disp('Setting up initial projection subspace');
    latent_der = latentproj(2:end) - latentproj(1:end-1);
    latent_der2 = latent_der(2:end) - latent_der(1:end-1);
    latent_der2_avg = conv([1 1 1 ]./3,abs(latent_der2));
    numprojcoeffs = find(abs(latent_der2_avg) < t,1) + 3;
    if (latentproj(numprojcoeffs-2) == 0)
        numprojcoeffs = find(latentproj == 0,1);
    end
    clear latentproj latent_der latent_der2 latent_der2_avg
    disp(['Number of projection basis elements: ' num2str(numprojcoeffs)]);
    meanproj = mean(proj_struct,3);
    numproj = size(proj_struct,3);
    projbasis = [meanproj(:), coeffproj(:,1:numprojcoeffs-1)];
    clear coeffproj meanproj;
    projcoeffs = [ones(numproj,1), scoreproj(:,1:numprojcoeffs-1)];
    clear scoreproj;    
    
    % Stuff for timing
    avgitertime = 0;
    avgssdtime = 0;
    
    % Initialize things
    disp('Initializing things'); tic;
    numim = size(imcoeffs,1);
    rots = rotstart:rotstep:rotend;
    numrot = length(rots);
    [x, y] = meshgrid(-transmax:transdelta:transmax,-transmax:transdelta:transmax);
    trans = [x(:), y(:)];
    clear x y
    numtrans = size(trans,1);
    if run == 1
        sigmastart = std(noisyims(:));
    end
    sigma1 = sigmastart;
    clear noisyims;
    sigma2 = sigma1*100;
    sigmaconst = 0;
    sigmathresh1 = 0.05;
    sigmathresh2 = 0.05;
    maskimcol = maskim(:);
    numpix = size(imbasis,1);
    numpixsqrt = sqrt(numpix);
    nummaskpix = sum(maskim(:));
    numctf = double(max(ctfinds));
    alphas = ones(numproj,1)/numproj;
    onesprojcoeff = ones(numprojcoeffs,1,'int8');
    onesproj = ones(numproj,1,'int8');
    proj_est = reshape(projbasis*projcoeffs'.*maskimcol(:,onesproj),[numpixsqrt,numpixsqrt,numproj]);
    eps = 1e-20;
    
    % Show initial approximated projections
    if dispflag
        figure; subplot(1,3,1); imshow(proj_est(:,:,1),[]);
        subplot(1,3,2); imshow(proj_est(:,:,floor(numproj/2)),[]);
        subplot(1,3,3); imshow(proj_est(:,:,end),[]);
        prc = prctile(structure(structure~=0),90);
        figure; plot_surface(structure,prc);
        pause(0.05);
    end
    
    % Set up initial search domains for translations
    disp('Set up translation search domains'); transtime = tic;
    searchtrans = get_trans_domains(1:numim,ceil(numtrans/2)*ones(numim,1),ones(numim,1),trans,floor(transmax/2),floor(transmax/2),numim);
    toc(transtime);
    
    % Compute image norms for translated mask
    if run == 1
        disp('Calc image norms'); calcnorm = tic;
        imnorms = comp_im_norms(imbasis,imcoeffs,maskim,trans,maxmem);
        toc(calcnorm);
    end
       
    %% Main loop
    for n = 1:maxnumiter
        
        itertime = tic;
        disp(' ');
        disp(['*****EM iter ' num2str(n) '*****']);
        
        % E-Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('E-Step');
        
        % Calculate projection norms squared
        disp('Calc norms of approximated projections'); pause(0.05); ssdtime = tic; tic;
        temp = reshape(imrotate(proj_est,45,'bilinear','crop'), [numpix numproj]);
        projnorms = dot(temp,temp)';
        clear temp;
        toc;
        
        % Compute inner products
        disp('Calc inner products'); pause(0.05); tic;
        ips = comp_inner_prods(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans);
        toc;
        
        % Calculate the SSDs
        disp('Calc SSDs'); pause(0.05); tic;
        [projinds,rotinds,iminds,lpcs_vals,SSDs,transinds] = comp_SSDs_fast(projnorms,projcoeffs,imcoeffs,ips,sigma1,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem);
        toc;
        clear ips;
        ssdtime = toc(ssdtime);
        avgssdtime = avgssdtime + ssdtime;
        ssdtimes(n) = ssdtime;
        
        % Update latent probabilities and related things for M-Step
        disp('Update E-Step probs'); pause(0.05); tic;
        lpcs_vals = alphas(projinds) .* lpcs_vals;
        tempnorm = zeros(numim,1);
        for k = 1:length(lpcs_vals)
            tempnorm(iminds(k)) = tempnorm(iminds(k)) + lpcs_vals(k);
        end
        lpcs_vals = lpcs_vals ./ (tempnorm(iminds) + eps);
        clear tempnorm;
        sumimrt = zeros(numproj,1);
        for k = 1:length(lpcs_vals)
            sumimrt(projinds(k)) = sumimrt(projinds(k)) + lpcs_vals(k);
        end
        toc;
        
        % M-Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('M-Step'); pause(0.05);
        
        proj_last = proj_est;
        sigma1_2 = sigma1^2;
        sigma2_2 = sigma2^2;
        
        % Update mixing coefficients
        disp('Update mixing coeffs and noise standard deviations'); tic;
        alphas = sumimrt ./ numim;
        
        % Update noise variances
        sigma1 = sqrt(1/nummaskpix/numim*(lpcs_vals(:)'*SSDs(:))) + sigmaconst;
        if sigma1 < sigmathresh1
            sigma1 = sigmathresh1 + sigmaconst;
        end
        if n > 1
            sigma2 = sqrt(1/nummaskpix/numproj*sum(sum(sum((proj_struct - proj_est).^2)))) + sigmaconst;
            if sigma2 < sigmathresh2
                sigma2 = sigmathresh2 + sigmaconst;
            end
        end
        disp(['sigma1 = ' num2str(sigma1)]);
        disp(['sigma2 = ' num2str(sigma2)]);
        toc;
        
        % Calculate weighted average of rotated images using E-step weights (ALL ON GPU)
        disp('Calc weighted average of rotated images'); pause(0.05); tic;
        avgrotim = avg_rot_ims(imbasis,imcoeffs,projinds,rotinds,iminds,lpcs_vals,rots,numpix,numim,numproj,numrot,numimcoeffs,numpixsqrt,trans,transinds);
        toc;
        
        % Update projection basis matrix
        disp('Update projection basis matrix'); pause(0.05); tic;
        projbasis = double((sigma2_2*avgrotim+sigma1_2*reshape(proj_struct,[numpix,numproj]))*(projcoeffs*inv(sigma2_2*projcoeffs'*(sumimrt(:,onesprojcoeff).*projcoeffs)+sigma1_2*projcoeffs'*projcoeffs)));
        toc;
        
         % Update projection image cofficients
        disp('Update projection image coefficients'); pause(0.05); tic;
        projcoeffs = double((inv(projbasis'*projbasis)*(projbasis'*(sigma2_2*avgrotim + sigma1_2*reshape(proj_struct,[numpix,numproj])))./(sigma2_2*sumimrt(:,onesprojcoeff)' + sigma1_2))');
        clear avgrotim proj_struct;
        toc;
        
        % Update structure
        disp('Update structure'); pause(0.05); tic;
        proj_est = reshape(projbasis*projcoeffs',[numpixsqrt,numpixsqrt,numproj]);
        fproj_est = zeros(size(proj_est));
        for j = 1:numproj
            fproj_est(:,:,j) = fftshift(fft2(proj_est(:,:,j)));
        end
        if n > 1
            structure = reconstruct_by_cg_w_ctf_par(fproj_est,data_axes,ctfs,mask,l_norm,l_smooth,iter_lim,stop_lim,structure);
        else
            structure = reconstruct_by_cg_w_ctf_par(fproj_est,data_axes,ctfs,mask,l_norm,l_smooth,iter_lim,stop_lim);
        end
        structure = structure.*structmask;
        clear fproj_est
        toc;
            
        % Update projections of structure
        disp('Update structure projections'); pause(0.05); tic;
        fproj_struct = project_in_all_directions_w_ctf_par(structure,mask,coord_axes,ctfs);
        proj_struct = zeros(size(proj_est));
        for j = 1:numproj
            proj_struct(:,:,j) = real(ifft2(ifftshift(fproj_struct(:,:,j)))).*maskim;
        end
        clear fproj_struct
        toc;
        proj_est = proj_est .* maskim(:,:,onesproj);
        projbasis = projbasis .* maskimcol(:,onesprojcoeff);
        
        % Update translation search
        if numtrans > 1
            disp('Update translation search domain'); tic;
            searchtrans = get_trans_domains(iminds,transinds,lpcs_vals,trans,transwidth,transdelta,numim);
            toc;
        end
        
        % Form the projections and check convergence
        disp('Check convergence'); tic;
        [done,err,pind] = check_convergence(proj_last,proj_est,convtol,numproj);
        if pind > 0
            disp(['Error: ' num2str(err) ' for proj ' num2str(pind)]);
        else
            disp('Converged!');
        end
        toc;
        
        % Timing stuff
        itertime = toc(itertime);
        itertimes(n) = itertime;
        disp(['Iteration time: ' num2str(itertime) ' seconds']);
        avgitertime = avgitertime + itertime;
        
        % Plot results
        if dispflag
            figure; subplot(1,3,1); imshow(proj_est(:,:,1),[]);
            subplot(1,3,2); imshow(proj_est(:,:,floor(numproj/2)),[]);
            subplot(1,3,3); imshow(proj_est(:,:,end),[]);
            figure; subplot(1,3,1); imshow(proj_struct(:,:,1),[]);
            subplot(1,3,2); imshow(proj_struct(:,:,floor(numproj/2)),[]);
            subplot(1,3,3); imshow(proj_struct(:,:,end),[]);
            if (n == 1)
                prc = prctile(structure(structure~=0),90);
            end
            figure; plot_surface(structure,prc);
            pause(0.05)
        end
        
        % Save iteration results
        if mod(n-1,1) == 0
            savename = ['run_' num2str(run) '_iter_' num2str(n)];
            save([pathout '/' savename],'-v7.3','structure','projbasis','projcoeffs','lpcs_vals','projinds','rotinds','iminds','transinds','searchtrans','err');
        end
        if done == 1
            break;
        end
        
    end
    
    % Calculate estimate using final E-step params and original noisy images
    disp(' '); disp('Calc final estimate!'); disp('Load original particle images'); tic;
    lpcs_vals = lpcs_vals ./ (sumimrt(projinds)+eps);
    noisyims = single(ReadMRC(imreconfile));
    
    %%%% FOR TESTING %%%%%%%%%%%%%%%%%%%%
    if substep > 0
        disp(['DIVIDING TOTAL DATA BY ' num2str(substep)]);
        noisyims = noisyims(:,:,1:substep:end);
        numim = size(noisyims,3);
    end
    
    if reconhalf
        noisyims = noisyims(:,:,reconstartind:2:end);
        numim = size(noisyims,3);
    end
    toc; 
    
    disp('Update final projection estimates'); tic;
    noisyims_g = gpuArray(noisyims);
    clear noisyims
    proj_est = update_templates2(noisyims_g,iminds,projinds,rotinds,transinds,lpcs_vals,rots,trans,numpixsqrt,numim,numproj);
    lpcs_vals = lpcs_vals .* sumimrt(projinds);
    toc;
    
    % Reconstruct
    disp('Final reconstruction'); tic;
    fproj_est = zeros(size(proj_est));
    for j = 1:numproj
        fproj_est(:,:,j) = fftshift(fft2(proj_est(:,:,j)));
    end
    numfp = size(fproj_est,3);
    takeoutinds = [];
    zeroim = zeros(numpixsqrt,numpixsqrt);
    for j = 1:numfp
        if isequal(fproj_est(:,:,j),zeroim)
            takeoutinds = [takeoutinds; j];
        end
    end
    keepinds = 1:numfp;
    keepinds(takeoutinds) = [];
    structure_final = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,keepinds),data_axes(:,keepinds),ctfs(:,:,mod(keepinds-1,numctf)+1),mask,l_norm,l_smooth,iter_lim,stop_lim);
    toc;
    
    % Display final results and stats summary
    if dispflag
        figure; subplot(1,3,1); imshow(proj_est(:,:,1),[]);
        subplot(1,3,2); imshow(proj_est(:,:,floor(numproj/2)),[]);
        subplot(1,3,3); imshow(proj_est(:,:,end),[]);
        prc = prctile(structure_final(structure_final~=0),90);
        figure; plot_surface(structure_final,prc);
    end
    avgitertime = avgitertime / n;
    avgssdtime = avgssdtime / n;
    runtime = toc(runtime);
    disp(['Average time per iteration: ' num2str(avgitertime) ' seconds']);
    disp(['Average time per SSD calc: ' num2str(avgssdtime) ' seconds']);
    disp(['Total wall time: ' num2str(runtime) ' seconds']);
    
    % Save
    dotinds = strfind(imfile,'.');
    slashinds = strfind(imfile,'/');
    if isempty(slashinds)
        slashinds = strfind(imfile,'\');
    end
    if isempty(slashinds)
        slashinds = 0;
    end
    savename = imfile(slashinds(end)+1:dotinds(end)-1);
    savename = [savename '_subspaceEM_' num2str(numimcoeffs) '_' num2str(numprojcoeffs) '_' num2str(rotstep) 'd_' num2str(transmax) num2str(transdelta) num2str(transwidth) 't_' num2str(run) 'x'];
    if substep > 0
        savename = [savename '_sub' num2str(substep)];
    end
    if reconhalf
        savename = [savename '_h' num2str(reconstartind)];
    end
    save([pathout '/' savename ,'.mat'],'-v7.3','structure_final','structure','proj_struct','proj_est','alphas','lpcs_vals','projinds','rotinds','rots','iminds','SSDs','totaltime','itertimes','avgitertime','ssdtimes','avgssdtime','n','sigma1','sigma2','projbasis','projcoeffs','searchtrans','transinds','trans','numprojcoeffs','numimcoeffs','t','imfile','initprojfile','imreconfile','convtol','coord_axes','structmask');
    
    % Some clean up
    clear structure_final structure proj_est proj_last
    
end

% Recon with largest mask possible for FSC calcs
disp('Reconstruct one more time with largest mask possible');
mask = get_mask_struct_ncd([numpixsqrt numpixsqrt numpixsqrt],1); % reconstruct with largest mask possible
recon = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,keepinds),data_axes(:,keepinds),ctfs(:,:,mod(keepinds-1,numctf)+1),mask,l_norm,l_smooth,iter_lim,stop_lim);
delete(gcp('nocreate'));
save([pathout '/' savename ,'.mat'],'-append','recon');
[~,h] = ReadMRC(imreconfile,1,-1);
writeMRC(recon,h.pixA,[pathout '/' 'subspaceEM_recon.mrc']);

if ~isequal(structmask,ones(size(structmask)))
    writeMRC(recon.*structmask,h.pixA,[pathout '/' 'subspaceEM_recon_masked.mrc']);
end

disp(['Total wall time for subspaceEM.m: ' num2str(toc(totaltime)/60/60) ' hrs']);

passed = 1;
