% Function to be used with fast_best_match code
% Calculates the sum of squared differences between projections and images
% and finds the best match for each image

function [projinds,rotinds,SSDs,transinds,scales] = comp_SSDs_fast_best_match(projnorms,projcoeffs,imcoeffs,ips,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem)
%function [projinds,rotinds,SSDs,transinds,scales] = comp_SSDs_fast_best_match(projnorms,projcoeffs,imcoeffs,ips_cache, num_chunks,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem)

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
                
                %nc = ips.nchunks;
                %nr = ips.dimensions(ips.broken);
                %dr = ceil(nr/nc);
                %dr_last = nr-(nc-1)*dr;
                %chunk = read_cached_array_chunk(ips, 1, [ips.dimensions(1), ips.dimensions(2), dr, ips.dimensions(4)]);
                % For each rotation
                for r = 1:numrot
                    
                    % First calculate the inner products between
                    % projections and current images
                    
                    %currips = currprojcoeffs*(ips(:,:,r,currt)*ic);
                    
                    %idx_m = ceil(r/dr);
                    %idx_d = mod(r,dr);
                    %if (idx_d == 0)
                    %    idx_d = dr;
                    %end
                    % DEBUG: %fprintf('r=%i, currt=%i, idx_m=%i, idx_d=%i\n', r, currt, idx_m, idx_d);
                    %if ( idx_m > ceil((r-1)/dr) && r > 1)
                    %    if (idx_m == ips.nchunks)
                    %        dr = dr_last;
                    %    end
                    %    chunk = read_cached_array_chunk(ips, idx_m, [ips.dimensions(1), ips.dimensions(2), dr, ips.dimensions(4)]);
                    %end
                    %DEBUG: %if (~isequal(ips_curr,ips(:,:,r,currt))) error('Failed matrix equality test when testing caching functions.'); end
                    %currips = currprojcoeffs*(chunk(:,:,idx_d,currt) * ic);
                    
                    currips = currprojcoeffs*(read_cached_array(ips, [0,0,r,currt])*ic);
                    
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
            
            % sort by rind so that to have sequensial access to ips
            [s_rind, i_rind] = sort(rind);
            s_pind = pind(i_rind);
            s_tind = tind(i_rind);
            s_curriminds = curriminds(i_rind);
            %dr = ceil(nr/nc);
            %chunk = read_cached_array_chunk(ips, 1, [ips.dimensions(1), ips.dimensions(2), dr, ips.dimensions(4)]);
            for i = 1:numcurrim
                scales(s_curriminds(i)) = currprojcoeffs(s_pind(i),:)* read_cached_array(ips, [0, 0, s_rind, currtrans(s_tind) ]) *...
                    imcoeffs(s_curriminds(i),:)' / projnormsc(s_pind(i))/2; 
                % scales(curriminds(i)) = currprojcoeffs(pind,:)*ips(:,:,rind,currtrans(tind))*imcoeffs(curriminds(i),:)' / projnormsc(pind) / 2;
                % Calculate scale that gave the min ssd
                %r = s_rind(i);
                %idx_m = ceil(r/dr);
                %idx_d = mod(r,dr);
                %if (idx_d == 0)
                %    idx_d = dr;
                %end
                %if ( idx_m > ceil((r-1)/dr) && r > 1)
                %    if (idx_m == numChunks)
                %        dr = dr_last;
                %    end
                %    chunk = read_cached_array_chunk(ips, idx_m, [ips.dimensions(1), ips.dimensions(2), dr, ips.dimensions(4)]);
                %end
                %scales(s_curriminds(i)) = currprojcoeffs(s_pind(i),:)* chunk(:,:,idx_d,currtrans(s_tind(i))) *...
                %    imcoeffs(s_curriminds(i),:)' / projnormsc(s_pind(i))/2;  
            end
            
            clear ssds currssds currprojnorms ic currimnorms
        end 
        progress_bar(st, numstu);
    end
end
fprintf('\n');

scales(scales < minscale) = minscale;
scales(scales > maxscale) = maxscale;

