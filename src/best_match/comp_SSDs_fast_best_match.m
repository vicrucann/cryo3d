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

% initialize distributor if needed
addpath(fullfile(cd, '../src/dist-wrappers'));
addpath(fullfile(cd, '../src/rshell-mat'));
path_vars = pathout;
currfold = pwd; 
cd('../src/rshell-mat/'); path_curr = pwd; path_curr = fixslash(path_curr);
cd(currfold);
debug = 1;

d = Distributor(login, path_rem, ipaddrs, path_vars, vars, path_curr, sleeptime, path_res, printout);
if (d.ncluster > 1)
    d.scp_cached_data(ips);
end

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
            currprojnorms = projnormsc(:,ones(numcurrim,1));
            ic = imcoeffs(curriminds,:)';
            
            %if d.ncluster == 1 % if distributor is not used
            if (debug)
                sum_currips = zeros(numrot, numst); % debug
                ssds = inf(numprojc,numcurrim,numrot,numst,'single');
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
                        sum_currips(r,t) = sum(sum(currips)); % debug, control sum!
                        
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
            end % debug   
            if debug % distributor debug, delete the whole section after done
                
                % =============================
                % split
                r_begin = zeros(1, d.ncluster);
                r_end = zeros(1, d.ncluster);
                dims = zeros(4, d.ncluster);
                numrot_ = ceil(numrot / d.ncluster);
                numrot_l = numrot - numrot_*(d.ncluster-1);
                for i=1:d.ncluster
                    rb = (i-1)*numrot_ + 1;
                    ds = ips.dimension();
                    if i~=d.ncluster
                        re = numrot_ * i;
                        ds(ips.ibroken())=numrot_;
                    else
                        re = numrot;
                        ds(ips.ibroken())=numrot_l;
                    end
                    r_begin(i) = rb;
                    r_end(i) = re;
                    dims(:,i) = ds;
                end
                
                % =============================
                % kernel
                vol = ips.window.volume(ips.ibroken());
                dimensions = ips.dimension();
                nvol = ceil(dimensions(ips.ibroken()) / vol);
                vol_end = dimensions(ips.ibroken()) - (nvol-1)*vol;
                dims = ips.window.volume;
                
                fprintf('The broken dimension of ips has size of %i or %i\n', vol, vol_end);
                fprintf('dims=%i %i %i %i\n', dims(1), dims(2), dims(3), dims(4));
                minindices = zeros(d.ncluster,numcurrim);
                minvalues = zeros(d.ncluster,numcurrim);
                sum_currips_di = zeros(numrot, numst);
                for i=1:d.ncluster
                    idxc = ceil(r_begin(i)/vol);
                    idxc_end = ceil(r_end(i)/vol);
                    fprintf('idxc=%i, idxc_end=%i\n', idxc, idxc_end);
                    file_dat = [ips.window.cpath ips.window.vname int2str(idxc) '.dat'];
                    mm = memmapfile(file_dat, 'Format', ips.type());
                    if (idxc == idxc_end)
                        vol=vol_end;
                        dims(ips.ibroken()) = vol_end;
                    end
                    ipsi = reshape(mm.Data, dims);
                    fprintf('Dat file is read\n');
                    
                    ssdi = inf(numprojc, numcurrim, r_end(i) - r_begin(i)+1, numst,'single');
                    fprintf('The calculation loop for r in range [%i %i] and t in range [%i %i]\n',...
                        r_begin(i), r_end(i), 1, numst);
                    
                    for r = r_begin(i):r_end(i)
                        for t = 1:numst
                            currt = currtrans(t);
                            if currt < 1
                                continue;
                            end
                            
                            % Set up the current image norms
                            currimnorms = imnorms(currt,curriminds);
                            currimnorms = currimnorms(onesprojc,:);
                            
                            % First calculate the inner products between
                            % projections and current images
                            idxc_ = ceil(r/vol);
                            if (idxc ~= idxc_) % memmapfile the next file
                                idxc = idxc_;
                                if (idxc == idxc_end)
                                    vol = vol_end;
                                    dims(ips.ibroken()) = vol_end;
                                end
                                fprintf('idxc=%i, r=%i, vol=%i\n', idxc, r, vol);
                                file_dat = [ips.window.cpath ips.window.vname int2str(idxc) '.dat'];
                                mm = memmapfile(file_dat, 'Format', ips.type());
                                ipsi = reshape(mm.Data, dims);
                                fprintf('Dat file is read\n');
                            end
                            
                            r_idx_loc = mod(r, vol);
                            if (r_idx_loc == 0)
                                r_idx_loc = vol;
                            end
                            
                            currips_di = currprojcoeffs*(ipsi(:, :, r_idx_loc, currt) * ic);
                            sum_currips_di(r,t) = sum(sum(currips_di));
                            if ~isequal(sum_currips(r,t), sum_currips_di(r,t))
                                fprintf('control sum is incorrect\n');
                            end
                            
                            % Calculate scale and adjust
                            s = currips_di ./ currprojnorms / 2;
                            s(s < minscale) = minscale;
                            s(s > maxscale) = maxscale;
                            % Calculate the ssds between each projection and image
                            ssdi(:, :, r - r_begin(i) + 1, t) = currimnorms + s.^2.*currprojnorms - s.*currips_di;
                        end
                    end
                    
                    fprintf('loop-1 terminated\n');
                    minidc = zeros(1,numcurrim);
                    minval = zeros(1,numcurrim);
                    for j = 1:numcurrim
                        currssdi = squeeze(ssdi(:,j,:,:));
                        [minval(j), minidc(j)] = min(currssdi(:));
                    end
                    fprintf('loop-2 terminated\n');
                    minindices(i,:) = minidc;
                    minvalues(i,:) = minval;    
                end
                if ~isequal(sum_currips, sum_currips_di)
                    fprintf('debug error!\n');
                end
                % =============================
                % merge
                minindices_glo = zeros(1,numcurrim);
                minvalues_glo = zeros(1,numcurrim);
                for i=1:numcurrim
                    [val, ind_loc] = min(minvalues(:,i));
                    minind_loc = minindices(ind_loc,i);
                    rot_size = ceil(numrot/d.ncluster);
                    if (ind_loc == d.ncluster)
                        rot_size = numrot - (d.ncluster-1)*rot_size;
                    end
                    [p_loc, r_loc, t_loc] = ind2sub([numprojc, rot_size, numst], minind_loc);
                    r_glo = r_loc + (ind_loc-1)*ceil(numrot/d.ncluster);
                    assert(r_glo <= numrot, 'r_glo calculation failed: out of range');
                    minindices_glo(i) = sub2ind([numprojc, numrot, numst], p_loc, r_glo, t_loc);
                    minvalues_glo(i) = val;
                end
                
                % =============================
                % out of distributor
                pind_di = zeros(1,numcurrim);
                rind_di = zeros(1,numcurrim);
                tind_di = zeros(1,numcurrim);
                for i = 1:numcurrim
                    val = minvalues_glo(i);
                    minind = minindices_glo(i);
                    %SSDs(curriminds(i)) = val; %val;
                    
                    [pind_di(i),rind_di(i),tind_di(i)] = ind2sub([numprojc,numrot,numst],minind);
                    if (pind_di(i) ~= pind(i) || rind_di(i) ~= rind(i) || tind_di(i) ~= tind(i))
                        fprintf('%i vs %i\n%i vs %i\n%i vs %i\n', pind_di(i), pind(i), rind_di(i), rind(i),...
                            tind_di(i), tind(i));
                        fprintf('numcurrim number: %i\n', i);
                    end
                    
                    projinds(curriminds(i)) = inds(pind(i));
                    rotinds(curriminds(i)) = rind(i);
                    transinds(curriminds(i)) = currtrans(tind(i));
                end
            end
            %else % if distributor is used
%                 in_split = struct('numst', numst, 'currtrans', currtrans, 'curriminds', curriminds, ...
%                     'onesprojc', onesprojc, 'currprojcoeffs', currprojcoeffs, 'ic', ic, ...
%                     'currprojnorms', currprojnorms, 'minscale', minscale, 'maxscale', maxscale, ...
%                     'imnorms', imnorms, 'numprojc', numprojc, 'numcurrim', numcurrim, 'numrot', numrot, ...
%                     'broken', ips.ibroken(), 'dimensions', ips.dimension(), 'ctype', ips.type(), 'volume', ips.window.volume,...
%                     'ncluster', d.ncluster, 'path_vars', path_vars, 'vars', vars);
%                 in_merge = struct('ncluster', d.ncluster, 'numcurrim', numcurrim, 'path_res', path_res, ...
%                     'vars', vars, 'numim', numim, 'curriminds', curriminds, ...
%                     'numprojc', numprojc, 'numrot', numrot, 'numst', numst, 'currtrans', currtrans,...
%                     'inds', inds);
%                 
%                 out = d.launch(@ssd_split, in_split, @ssd_wrap, @ssd_merge, in_merge);
%                 minindices = out.minindices;
%                 minvalues = out.minvalues;
%                 % For each image in the batch
%                 pind_di = zeros(1,numcurrim);
%                 rind_di = zeros(1,numcurrim);
%                 tind_di = zeros(1,numcurrim);
%                 for i = 1:numcurrim
%                     val = minvalues(i);
%                     minind = minindices(i);
%                     SSDs(curriminds(i)) = val; %val;                    
%                     
%                     [pind_di(i),rind_di(i),tind_di(i)] = ind2sub([numprojc,numrot,numst],minind);
%                     if (pind_di(i) ~= pind(i) || rind_di(i) ~= rind(i) || tind_di(i) ~= tind(i))
%                         fprintf('%i vs %i\n%i vs %i\n%i vs %i\n', pind_di(i), pind(i), rind_di(i), rind(i),...
%                             tind_di(i), tind(i));
%                         fprintf('numcurrim number: %i\n', i);
%                         error('indexing calculation failed');
%                     end
%                         
%                     projinds(curriminds(i)) = inds(pind(i));
%                     rotinds(curriminds(i)) = rind(i);
%                     transinds(curriminds(i)) = currtrans(tind(i));
%                 end
            %end
            
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

