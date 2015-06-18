% Function to be used with fast_best_match code
% Calculates the sum of squared differences between projections and images
% and finds the best match for each image
% The function is edited to be compatible with caching data structure

function [projinds,rotinds,SSDs,transinds,scales] = comp_SSDs_fast_best_match(projnorms,projcoeffs,imcoeffs,ips,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem,...
    ipaddrs,login,path_rem,vars,sleeptime,path_res,printout,pathout)

% Initializations
projinds = -ones(numim,1);
rotinds = zeros(numim,1);
SSDs = zeros(numim,1);
transinds = zeros(numim,1);
scales = zeros(numim,1);
numst = size(searchtrans,1);
searchtrans = searchtrans';
currmem = monitor_memory_whos;
minscale = 1.0;
maxscale = 1.0;

% For each CTF class
for c = 1:numctf
    
    % Get the relevant indices and norms for the current CTF
    inds = c:numctf:numproj;
    numprojc = length(inds);
    onesprojc = ones(numprojc,1,'int8');
    projnormsc = projnorms(inds);
    currprojcoeffs = projcoeffs(inds,:);
    cis = find(ctfinds == c)';
    searchtransu = unique(searchtrans(cis,:),'rows');
    numstu = size(searchtransu,1);
   
    % For each set of translations
    fprintf('\nNumber of set of translations: %i\n', numstu);
    for st = 1:numstu
        currtrans = searchtransu(st,:);
        sinds = find(ismember(searchtrans(cis,:),currtrans,'rows'));
        
        % Determine number of images to process per batch for handling limited memory
        numsinds = length(sinds);
        numbatches = ceil((numprojc*numsinds*numrot*numst + numprojc*numrot*numst)*4/1048576 / (maxmem - currmem - 550));
        if numbatches < 1
            numbatches = numsinds;
        elseif numbatches > 1
            numbatches = numbatches + 1; %%%% Changed to deal with memory
        end
        batchsize = ceil(numsinds / numbatches);
        numbatches = ceil(numsinds / batchsize);
        %fprintf('\nNumber of images per batch to process: %i\n', numbatches);
        
        % For each batch of images
        for b = 1:numbatches
            
            % Get the indices of the images
            if b == numbatches
                curriminds = cis(sinds((b-1)*batchsize+1:end));
            else
                if (b*batchsize > length(sinds))
                    ind = b*batchsize
                    length(sinds)
                end
                curriminds = cis(sinds((b-1)*batchsize+1:b*batchsize));
            end
            numcurrim = length(curriminds);
            
            % Setup ssds matrix, rep proj norms, and get current imcoeffs
            ssds = inf(numprojc,numcurrim,numrot,numst,'single');
            currprojnorms = projnormsc(:,ones(numcurrim,1));
            ic = imcoeffs(curriminds,:)';
            
            % if using distributor, break data
            [ncluster ~] = find(ipaddrs==' ');
            ncluster = size(ncluster,2)+1;
            if ncluster == 1 % if distributor is not used
                % For each rotation
                for r = 1:numrot
                    % For each translation in the current search set
                    for t = 1:numst
                        % Check if translation exists
                        currt = currtrans(t);
                        if currt < 1
                            continue;
                        end
                        
                        % Set up the current image norms
                        currimnorms = imnorms(currt,curriminds);
                        currimnorms = currimnorms(onesprojc,:);
                        
                        % First calculate the inner products between
                        % projections and current images
                        %currips = currprojcoeffs*(ips.read_cached_array([0,0,r,currt])*ic);
                        currips = currprojcoeffs*(ips(:,:,r,currt)*ic);
                        
                        %                   % Calculate scale and adjust
                        s = currips ./ currprojnorms / 2;
                        s(s < minscale) = minscale;
                        s(s > maxscale) = maxscale;
                        
                        % Calculate the ssds between each projection and image
                        ssds(:,:,r,t) = currimnorms + s.^2.*currprojnorms - s.*currips;
                    end
                end
                
                % For each image in the batch
                pind = zeros(1,numcurrim);
                rind = zeros(1,numcurrim);
                tind = zeros(1,numcurrim);
                for i = 1:numcurrim
                    % Find the min ssd and get the indices of the parameters of
                    % the best batch
                    currssds = squeeze(ssds(:,i,:,:));
                    [SSDs(curriminds(i)),minind] = min(currssds(:));
                    [pind(i),rind(i),tind(i)] = ind2sub([numprojc,numrot,numst],minind);
                    projinds(curriminds(i)) = inds(pind(i));
                    rotinds(curriminds(i)) = rind(i);
                    transinds(curriminds(i)) = currtrans(tind(i));
                end
                
            else % if distributor is used
                addpath(fullfile(cd, '../src/dist-wrappers'));
                addpath(fullfile(cd, '../src/rshell-mat'));
                fprintf('\nCalculation using remotes\n');
                
                currfold = pwd; cd('../src/rshell-mat/');
                path_curr = pwd; path_curr = fixslash(path_curr);
                cd(currfold); cd('../src/best_match/');
                srcfold = pwd; srcfold = fixslash(srcfold);
                cd(currfold);
                pathsrc = srcfold;
                path_vars = pathout;
                
                path_cache = ips.window.cpath;
                cachevar = ips.window.vname;
                d = Distributor(login, path_rem, ipaddrs, path_vars, vars, path_cache, cachevar,...
                    path_curr, sleeptime, path_res, printout);
                in_split = struct('numst', numst, 'currtrans', currtrans, 'curriminds', curriminds, ...
                    'onesprojc', onesprojc, 'currprojcoeffs', currprojcoeffs, 'ic', ic, ...
                    'currprojnorms', currprojnorms, 'minscale', minscale, 'maxscale', maxscale, ...
                    'imnorms', imnorms, 'numprojc', numprojc, 'numcurrim', numcurrim, 'numrot', numrot, ...
                    'broken', ips.ibroken(), 'dimensions', ips.dimension(), 'ctype', ips.type(), ...
                    'ncluster', d.ncluster, 'path_vars', path_vars, 'vars', vars);
                in_merge = struct('ncluster', d.ncluster, 'numcurrim', numcurrim, 'path_res', path_res, ...
                    'vars', vars, 'numim', numim, 'curriminds', curriminds, ...
                    'numprojc', numprojc, 'numrot', numrot, 'numst', numst, 'currtrans', currtrans,...
                    'inds', inds);
                
                out = d.launch(@ssd_split, in_split, @ssd_wrap, @ssd_merge, in_merge);
                minindices = out.minindices;
                minvalues = out.minvalues;
                % For each image in the batch
                pind = zeros(1,numcurrim);
                rind = zeros(1,numcurrim);
                tind = zeros(1,numcurrim);
                for i = 1:numcurrim
                    [val, ind] = min(minvalues(:,i));
                    SSDs(curriminds(i)) = val;
                    minind = minindices(ind,i);
                    [pind(i),rind(i),tind(i)] = ind2sub([numprojc,numrot,numst],minind);
                    projinds(curriminds(i)) = inds(pind(i));
                    rotinds(curriminds(i)) = rind(i);
                    transinds(curriminds(i)) = currtrans(tind(i));
                end
            end
            
            % to calculate scales, need to sort by rind so that to have sequensial access to ips
            [s_rind, i_rind] = sort(rind);
            s_pind = pind(i_rind);
            s_tind = tind(i_rind);
            s_curriminds = curriminds(i_rind);
            for i = 1:numcurrim
                scales(s_curriminds(i)) = currprojcoeffs(s_pind(i),:) * ips(:,:,s_rind(i),currtrans(s_tind(i)))  *...
                    imcoeffs(s_curriminds(i),:)' / projnormsc(s_pind(i))/2;
            end
            
            clear ssds currssds currprojnorms ic currimnorms
        end
        progress_bar(st, numstu);
    end
    fprintf('\n');
end

scales(scales < minscale) = minscale;
scales(scales > maxscale) = maxscale;
end

function foldname = fixslash(foldname)
slash = foldname(end);
if (~isequal(slash, '\') && ~isequal(slash, '/'))
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % Windows
        foldname = [foldname '\'];
    else % Linux
        foldname = [foldname '/'];
    end
end
end

