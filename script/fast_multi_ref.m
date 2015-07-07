% Script for running multiparticle reconstructions using subspace
% approximations

% Created by Nicha C. Dvornek, 12/2014
% Last modified 06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function passed = fast_multi_ref(pathout, configfile)

clear; reset(gpuDevice);

% Add paths
addpath(fullfile(cd, '../src/recon'));
addpath(fullfile(cd, '../src/mrc'));

% Stuff for timing
totaltime = tic;

% Set up - user parameters
[imfile,imreconfile,ctffile,maxmem,numthreads,dispflag,substep,reconhalf,reconstartind,normprojint,numruns,maxnumiter,rotstart,rotstep,rotend,transmax,transdelta,transwidth,convtol,t,numstruct] = read_config_file_multiref(configfile);

% Reconstruction parameters
l_norm = 0.1;
l_smooth = -10;
iter_lim = 10;
stop_lim = 0.03;

% Start pool
% if numthreads > 0  % matlab2013
%    if matlabpool('size') > 0
%         matlabpool close;
%     end
%     matlabpool('local',numthreads);
% end
if (~isempty (gcp('nocreate')) ) % matlab 2014, may not be needed
    delete(gcp('nocreate'));
end
parpool(numthreads);

%% Outer loop
for run = 1:numruns
    
    totaltime = tic;
    disp(' ');
    disp(['Run ' num2str(run)]);
    
    reset(gpuDevice);
    
    % Load image-related things
    if run == 1 
        disp('Loading images and setting up image subspace'); tic;
        
        % Check if pca of images already exists
        pcafile = [imfile(1:strfind(imfile,'.')-1) '_pca.mat'];
        if ~exist(pcafile,'file');
            % Do PCA
            data = single(ReadMRC(imfile));
            data = reshape(data,[size(data,1)*size(data,2),size(data,3)])';  
            disp('PCA of images');
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
        [noisyims,h] = ReadMRC(imfile);
        noisyims = single(noisyims);
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
    
    disp('Set up initial structure and projections'); tic;
    if run == 1
        % Get initial model parameters      
        initprojfile = make_init_model(configfile,pathout);
    else
        writefile = [savename '_reinit.mat'];  
        proj_init_from_file_multiref(savename,initprojfile,writefile,0);
        initprojfile = writefile;        
    end
    
    disp('Loading initial structure and projections');
    clear proj_struct
    load(initprojfile,'maskim','mask','coord_axes','data_axes','ctfs','structure','proj_struct','structmask'); 
    
    if run == 1 && normprojint   
        % Rescale projection intensities to match images
        disp('Rescaling projection intensities');
        proj_struct = rescale_proj_ints(proj_struct,noisyims); 
    end

    disp('PCA of initial projections');
    data = reshape(proj_struct,[size(proj_struct,1)*size(proj_struct,2),size(proj_struct,3)])';
    pcatic = tic; [coeffproj, scoreproj, latentproj] = princomp(data); toc(pcatic);
    clear data;
    
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
    sigma2 = sigma1*100;
    sigmaconst = 0;
    sigmathresh1 = 0.05;
    sigmathresh2 = 0.05;
    maskimcol = maskim(:);
    numpix = size(imbasis,1);
    numpixsqrt = sqrt(numpix);
    nummaskpix = sum(maskim(:));
    numctf = double(max(ctfinds));
    iminds = 1:numim;
    iminds = iminds(ones(numstruct,1),:);
    iminds = iminds(:);
    structinds = repmat((1:numstruct)',[numim,1]);
    
    eps = 1e-20;
    numprojstruct = size(coord_axes,2);
    lpcs_vals = 1./numstruct.*ones(numim*numstruct,1);
        
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
    
    %% Generate random initial models from single model
    if run == 1
        disp('Generate random initial 3D models');
        randinds = randperm(numim);
        numimperstruct = floor(numim/numstruct);

        % Calculate projection norms squared
        disp('Calc approximated projection norms'); pause(0.05); tic;
        proj_est = reshape(projbasis*projcoeffs'.*maskimcol(:,ones(numproj,1)),[numpixsqrt,numpixsqrt,numproj]);
        temp = reshape(imrotate(proj_est,45,'bilinear','crop'), [numpix numproj]);
        projnorms = dot(temp,temp)';
        clear temp;
        toc;

        % Compute inner products
        disp('Calc inner products'); pause(0.05); tic;
        ips = comp_inner_prods_nocache(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans);
        toc;

        % For each seed, get random subset of images and do one round of
        % alignment, then do reconstruction to get starting model
        projinds = zeros(numim,1);
        rotinds = zeros(numim,1);
        transinds = zeros(numim,1);
        for s = 1:numstruct

            if s == numstruct
                currinds = randinds(numimperstruct*(s-1)+1:end);
            else
                currinds = randinds(numimperstruct*(s-1)+1:numimperstruct*s);
            end
            numcurrims = length(currinds);

            disp(['Calc SSDs for initial model ' num2str(s)]); pause(0.05);tic;
            [pinds,rinds,SSDs,tinds,scales] = comp_SSDs_fast_best_match_nocache(projnorms,projcoeffs,imcoeffs(currinds,:),ips,ctfinds(currinds),length(currinds),numctf,numprojstruct,numrot,searchtrans(:,currinds),imnorms(:,currinds),maxmem);
            projinds(currinds) = pinds;
            rotinds(currinds) = rinds;
            transinds(currinds) = tinds;
            toc;
            
            disp('Calc image weights and update projection templates'); tic;
            noisyims_g = gpuArray(noisyims(:,:,currinds));
            weights = zeros(numcurrims,1);
            for j = 1:numprojstruct
                inds = find(pinds == j);
                if length(inds > 0)
                    weights(inds) = scales(inds) ./ sum(scales(inds).^2);
                end
            end
            proj_est = update_templates2(noisyims_g,(1:numcurrims)',pinds,rinds,tinds,weights,rots,trans,numpixsqrt,numcurrims,numprojstruct);
            toc;

            disp(['Reconstruct initial model ' num2str(s)]); tic;
            fproj_est = zeros(size(proj_est));
            for j = 1:numprojstruct
                fproj_est(:,:,j) = fftshift(fft2(proj_est(:,:,j)));
            end
            takeoutinds = [];
            zeroim = zeros(numpixsqrt,numpixsqrt);
            % Take out directions that have no images to speed up reconstruction
            for j = 1:numprojstruct
                if isequal(fproj_est(:,:,j),zeroim)
                    takeoutinds = [takeoutinds; j];
                end
            end
            keepinds = 1:numprojstruct;
            keepinds(takeoutinds) = [];
            structure(:,:,:,s) = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,keepinds),data_axes(:,keepinds),ctfs(:,:,mod(keepinds-1,numctf)+1),mask,l_norm,l_smooth,iter_lim,stop_lim);
            structure(:,:,:,s) = lpf_at_res_sigma(structure(:,:,:,s),h.pixA,h.pixA*3,1); %%% TESTING LPF %%% 
            structure(:,:,:,s) = structure(:,:,:,s).*structmask;
            toc; 
            
            disp(['Get projections of initial model ' num2str(s)]); tic;
            fproj_struct = project_in_all_directions_w_ctf_par(structure(:,:,:,s),mask,coord_axes,ctfs);
            for j = 1:numprojstruct
                proj_struct(:,:,(s-1)*numprojstruct+j) = real(ifft2(ifftshift(fproj_struct(:,:,j)))).*maskim;
            end
            toc;
        end
        save('init_rand_structures.mat','structure','randinds');
        
        % Do PCA on the projections from all the structures
        disp('PCA of all projections'); tic;
        clear projbasis projcoeffs
        numproj = size(proj_struct,3);
        data = reshape(proj_struct,[size(proj_struct,1)*size(proj_struct,2),size(proj_struct,3)])';
        [coeffproj, scoreproj, latentproj] = pca(data);
        clear data;
        toc;

        % Setup projection subspace
        disp('Setting up projection subspace'); tic;
        latent_der = latentproj(2:end) - latentproj(1:end-1);
        latent_der2 = latent_der(2:end) - latent_der(1:end-1);
        latent_der2_avg = conv([1 1 1 ]./3,abs(latent_der2));
        numprojcoeffs = find(abs(latent_der2_avg) < t,1) + 3;
        if (latentproj(numprojcoeffs-2) == 0)
            numprojcoeffs = find(latentproj == 0,1);
        end
        clear latentproj latent_der latent_der2 latent_der2_avg
        disp(['Number of template basis elements: ' num2str(numprojcoeffs)]);
        meanproj = mean(proj_struct,3);
        projbasis = [meanproj(:), coeffproj(:,1:numprojcoeffs-1)];
        clear coeffproj meanproj;
        projcoeffs = [ones(numproj,1), scoreproj(:,1:numprojcoeffs-1)];
        clear scoreproj;    
        toc;
    
    end
    clear noisyims noisyims_g
    
    % Some final initializations
    onesprojcoeff = ones(numprojcoeffs,1,'int8');
    onesproj = ones(numproj,1,'int8');
    proj_est = reshape(projbasis*projcoeffs'.*maskimcol(:,onesproj),[numpixsqrt,numpixsqrt,numproj]);
    alphas = ones(numstruct,1)/numstruct;
    
    % Get first set of SSDs
    disp('Get SSDs from initial alignments'); tic;
    SSDs = zeros(numim*numstruct,1);
    ind = 1;
    if run > 1
        projinds = projinds(1:numstruct:end);
    end
    for i = 1:numim
        approxim = imbasis*imcoeffs(i,:)';
        for s = 1:numstruct
            currproj = imrotate(proj_est(:,:,projinds(i)+(s-1)*numprojstruct),rots(rotinds(i)),'bilinear','crop');
            dx = trans(transinds(i),1);
            dy = trans(transinds(i),2);
            currprojtrans = zeros(numpixsqrt,numpixsqrt,'single');
            if dy < 0
                if dx < 0
                    currprojtrans(1:end+dy,1:end+dx) = currproj(1-dy:end,1-dx:end);
                else
                    currprojtrans(1:end+dy,1+dx:end) = currproj(1-dy:end,1:end-dx);
                end
            else
                if dx < 0
                    currprojtrans(1+dy:end,1:end+dx) = currproj(1:end-dy,1-dx:end);
                else
                    currprojtrans(1+dy:end,1+dx:end) = currproj(1:end-dy,1:end-dx);
                end
            end

            SSDs(ind) = sum((approxim(:) - currprojtrans(:)).^2);
            ind = ind+1;
        end
    end   
    clear approxim currproj
    toc;
    
     %% Main loop
     for n = 1:maxnumiter
         
        itertime = tic;
        disp(' ');
        disp(['*****EM iteration ' num2str(n) '*****']);
        
        disp('E-Step');
  
        % Update latent probabilities and related things for M-Step
        disp('Update E-Step probs'); pause(0.05); tic;
        minSSDs = min(reshape(SSDs,[numstruct,numim]),[],1)';
        lpcs_vals = alphas(structinds) .* exp(- (SSDs - minSSDs(iminds(:)))./ (2*sigma1^2));
        clear minSSDs
        tempnorm = zeros(numim,1);
        for k = 1:length(lpcs_vals)
            tempnorm(iminds(k)) = tempnorm(iminds(k)) + lpcs_vals(k);
        end
        lpcs_vals = lpcs_vals ./ (tempnorm(iminds) + eps);
        clear tempnorm;
        sumoverim = zeros(numstruct,1);           
        for k = 1:length(lpcs_vals)
            sumoverim(structinds(k)) = sumoverim(structinds(k)) + lpcs_vals(k);
        end
        
        toc;
        
        % M-Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('M-Step');
        
        proj_last = proj_est;
        sigma1_2 = sigma1^2;
        sigma2_2 = sigma2^2;
        
        % Update mixing coefficients
        disp('Update mixing coeffs and noise standard deviations'); tic;
        alphas = sumoverim ./ numim;
        
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
        
        % Update projection direction and transformation estimates
        % Calculate projection norms squared
        disp('Update projection direction and transformation estimates');
        disp('Calc template norms'); pause(0.05); ssdtime = tic; tic;
        temp = reshape(imrotate(proj_est,45,'bilinear','crop'), [numpix numproj]);
        projnorms = dot(temp,temp)';
        clear temp;
        toc;
        
        % Compute inner products
        disp('Calc inner products'); pause(0.05); tic;
        ips = comp_inner_prods_nocache(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans);
        toc;
        
         % Calculate the SSDs
        disp('Calc SSDs'); pause(0.05);tic;
        [projinds,rotinds,SSDs,transinds] = comp_SSDs_fast_multi_ref(projnorms,projcoeffs,imcoeffs,ips,ctfinds,lpcs_vals,numim,numctf,numstruct,numproj,numrot,searchtrans,imnorms,maxmem);
        toc;
        clear ips;
        ssdtime = toc(ssdtime);
        ssdtimes(n) = ssdtime;
        
        % Calculate weighted average of rotated images using E-step weights (ALL ON GPU)
        disp('Calc weighted average of rotated images'); pause(0.05); tic;
        avgrotim = avg_rot_ims(imbasis,imcoeffs,projinds,rotinds,iminds,lpcs_vals,rots,numpix,numim,numproj,numrot,numimcoeffs,numpixsqrt,trans,transinds);
        toc;
        
        % Update template matrix
        disp('Update template matrix'); pause(0.05); tic;
        sumoverimperproj = zeros(numproj,1);
        for k = 1:length(lpcs_vals)
            sumoverimperproj(projinds(k)) = sumoverimperproj(projinds(k)) + lpcs_vals(k);
        end
        projbasis = double((sigma2_2*avgrotim+sigma1_2*reshape(proj_struct,[numpix,numproj]))*(projcoeffs*inv(sigma2_2*projcoeffs'*(sumoverimperproj(:,onesprojcoeff).*projcoeffs)+sigma1_2*projcoeffs'*projcoeffs)));
        toc;
        
        % Update template cofficients
        disp('Update template coeffs'); pause(0.05); tic;
        projcoeffs = double((inv(projbasis'*projbasis)*(projbasis'*(sigma2_2*avgrotim + sigma1_2*reshape(proj_struct,[numpix,numproj])))./(sigma2_2*sumoverimperproj(:,onesprojcoeff)' + sigma1_2))');
        clear avgrotim proj_struct;
        toc;
        
        % Update structure
        disp('Update templates and structures'); pause(0.05); tic;
        proj_est = reshape(projbasis*projcoeffs',[numpixsqrt,numpixsqrt,numproj]);
        fproj_est = zeros(size(proj_est));
        for j = 1:numproj
            fproj_est(:,:,j) = fftshift(fft2(proj_est(:,:,j)));
        end
        for s = 1:numstruct
            disp(['Reconstructing structure ' num2str(s)]);
            if n > 1
                structure(:,:,:,s) = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,(s-1)*numprojstruct+1:s*numprojstruct),data_axes,ctfs,mask,l_norm,l_smooth,iter_lim,stop_lim,structure(:,:,:,s));
            else
                structure(:,:,:,s) = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,(s-1)*numprojstruct+1:s*numprojstruct),data_axes,ctfs,mask,l_norm,l_smooth,iter_lim,stop_lim);
            end
            structure(:,:,:,s) = structure(:,:,:,s).*structmask;
        end
        clear fproj_est
        toc;
        
        % Update projections of structure
        disp('Update structure projections'); pause(0.05); tic;
        for s = 1:numstruct
            disp(['Projecting structure ' num2str(s)]);
            fproj_struct(:,:,(s-1)*numprojstruct+1:s*numprojstruct) = project_in_all_directions_w_ctf_par(structure(:,:,:,s),mask,coord_axes,ctfs);
        end
        proj_struct = zeros(size(proj_est));
        for j = 1:numproj
            proj_struct(:,:,j) = real(ifft2(ifftshift(fproj_struct(:,:,j)))).*maskim;
        end
        clear fproj_struct
        toc;
        
        % Masking for next round
        proj_est = proj_est .* maskim(:,:,onesproj);
        projbasis = projbasis .* maskimcol(:,onesprojcoeff);
        
        % Update translation search
        if numtrans > 1
            disp('Update translation search domain'); tic;
            searchtrans = get_trans_domains(iminds,transinds,lpcs_vals,trans,transwidth,transdelta,numim);
            toc;
        end
        
        % Check convergence
        disp('Check convergence'); tic;
        [done,err,pind] = check_convergence(proj_last,proj_est,convtol,numproj);
        if pind > 0
            disp(['Error: ' num2str(err) ' for proj ' num2str(pind)]);
        else
            disp('Converged!');
        end
        toc;
        
        % Timing stuff
        itertimes(n) = toc(itertime);
        disp(['Iteration time: ' num2str(itertimes(n)) ' seconds']);
        
        % Save iteration results
        if mod(n-1,1) == 0
            savename = ['run_' num2str(run) '_iter_' num2str(n)];
            save(savename,'-v7.3','structure','projbasis','projcoeffs','lpcs_vals','projinds','rotinds','transinds','scales','searchtrans','err');
        end
        
        if done == 1
            break;
        end
        
    end
    
    % Calculate reconstruction using final alignment params and original noisy images
    disp(' '); disp('Calc final estimate!'); disp('Load original particle images'); tic;
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
    
    disp('Update final projection templates'); tic;
    noisyims_g = gpuArray(noisyims);
    clear noisyims
    proj_est = update_templates2(noisyims_g,iminds,projinds,rotinds,transinds,lpcs_vals./(sumoverimperproj(projinds)+eps),rots,trans,numpixsqrt,numim,numproj);
    
    disp('Final reconstruction'); tic;
    fproj_est = zeros(size(proj_est));
    for j = 1:numproj
        fproj_est(:,:,j) = fftshift(fft2(proj_est(:,:,j)));
    end
    takeoutinds = [];
    zeroim = zeros(numpixsqrt,numpixsqrt);
    for j = 1:numproj
        if isequal(fproj_est(:,:,j),zeroim)
            takeoutinds = [takeoutinds; j];
        end
    end
    length(takeoutinds)
    keepinds = 1:numproj;
    keepinds(takeoutinds) = [];
    structure_final = zeros(size(structure));
    for s = 1:numstruct
        disp(['Reconstructing structure ' num2str(s)]);
        currinds = keepinds(keepinds >= (s-1)*numprojstruct+1 & keepinds <= s*numprojstruct);
        structure_final(:,:,:,s) = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,currinds),data_axes(:,currinds-(s-1)*numprojstruct),ctfs(:,:,mod(currinds-1,numctf)+1),mask,l_norm,l_smooth,iter_lim,stop_lim);
    end
    toc;
        
    % Display stats summary
    totaltime = toc(totaltime);
    disp(['Average time per serial iteration: ' num2str(mean(itertimes)) ' seconds']);
    disp(['Average time per serial SSD calc: ' num2str(mean(ssdtimes)) ' seconds']);
    disp(['Total wall time: ' num2str(totaltime) ' seconds']);
    
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
    savename = [savename '_fastmr_' num2str(numimcoeffs) '_' num2str(numprojcoeffs) '_' num2str(rotstep) 'd_' num2str(transmax) num2str(transdelta) num2str(transwidth) 't_' num2str(run) 'x'];
    if substep > 0
        savename = [savename '_sub' num2str(substep)];
    end
    if reconhalf
        savename = [savename '_h' num2str(reconstartind)];
    end
    save([pathout '/' savename '.mat'],'-v7.3','structure_final','structure','proj_struct','proj_est','alphas','lpcs_vals','projinds','rotinds','rots','SSDs','totaltime','itertimes','ssdtimes','n','sigma1','sigma2','projbasis','projcoeffs','searchtrans','transinds','trans','scales','numprojcoeffs','numimcoeffs','t','imfile','initprojfile','imreconfile','convtol');
    
    % Some clean up
    clear structure_final structure proj_est proj_last
end

% Recon with largest mask possible for fsc calcs
disp('Reconstruct one more time with largest mask possible');
mask = get_mask_struct_ncd([numpixsqrt numpixsqrt numpixsqrt],1); % reconstruct with largest mask possible
recon = zeros(numpixsqrt,numpixsqrt,numpixsqrt,numstruct);
[~,h] = ReadMRC(imreconfile,1,-1);
for s = 1:numstruct
    disp(['Reconstructing structure ' num2str(s)]);
    currinds = keepinds(keepinds >= (s-1)*numprojstruct+1 & keepinds <= s*numprojstruct);
    recon(:,:,:,s) = reconstruct_by_cg_w_ctf_par(fproj_est(:,:,currinds),data_axes(:,currinds-(s-1)*numprojstruct),ctfs(:,:,mod(currinds-1,numctf)+1),mask,l_norm,l_smooth,iter_lim,stop_lim);
    writeMRC(recon(:,:,:,s),h.pixA,['fmr_recon_' num2str(s) '.mrc']);
    if ~isequal(structmask,ones(size(structmask)))
        writeMRC(recon(:,:,:,s).*structmask,h.pixA,[pathout '/' 'fmr_recon_' num2str(s) '_masked.mrc']);
    end
end
% matlabpool close  % matlab2013
delete(gcp('nocreate')); 
save(savename,'-append','recon');

passed = 1;
