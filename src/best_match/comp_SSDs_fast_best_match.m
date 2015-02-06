% Function to be used with fast_best_match code
% Calculates the sum of squared differences between projections and images
% and finds the best match for each image

function [projinds,rotinds,SSDs,transinds,scales] = comp_SSDs_fast_best_match(projnorms,projcoeffs,imcoeffs,ips,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem)

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
                
                % For each rotation
                for r = 1:numrot
                    
                    % First calculate the inner products between
                    % projections and current images
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
            for i = 1:numcurrim
                
                % Find the min ssd and get the indices of the parameters of
                % the best batch
                currssds = squeeze(ssds(:,i,:,:));
                [SSDs(curriminds(i)),minind] = min(currssds(:));
                [pind,rind,tind] = ind2sub([numprojc,numrot,numst],minind);
                projinds(curriminds(i)) = inds(pind);
                rotinds(curriminds(i)) = rind;
                transinds(curriminds(i)) = currtrans(tind);
                
                % Calculate scale that gave the min ssd
                scales(curriminds(i)) = currprojcoeffs(pind,:)*ips(:,:,rind,currtrans(tind))*imcoeffs(curriminds(i),:)' / projnormsc(pind) / 2;
                
            end
            clear ssds currssds currprojnorms ic currimnorms
        end    
    end
end

scales(scales < minscale) = minscale;
scales(scales > maxscale) = maxscale;

