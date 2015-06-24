% Function to be used with fast_multi_ref.m to calculate SSDs

% Created by Nicha C. Dvornek, 01/2015
% Last modified 03/2015

function [projinds,rotinds,SSDs,transinds] = comp_SSDs_fast_multi_ref(projnorms,projcoeffs,imcoeffs,ips,ctfinds,lpcs_vals,numim,numctf,numstruct,numproj,numrot,searchtrans,imnorms,maxmem)

projinds = -ones(numim,1);
rotinds = zeros(numim,1);
SSDs = zeros(numim,1);
transinds = zeros(numim,1);
numst = size(searchtrans,1);
searchtrans = searchtrans';
currmem = monitor_memory_whos;
numdir = numproj/numstruct/numctf;
onesdir = ones(numdir,1,'int8');
lps = reshape(lpcs_vals,[numstruct,numim]);

for c = 1:numctf
    inds = c:numctf:numproj;
    projnormsc = projnorms(inds);
    currprojcoeffs = projcoeffs(inds,:);
    cis = find(ctfinds == c)';
    searchtransu = unique(searchtrans(cis,:),'rows');
    numstu = size(searchtransu,1);
    for st = 1:numstu
        currtrans = searchtransu(st,:);
        sinds = find(ismember(searchtrans(cis,:),currtrans,'rows'));
        % DETERMINE BATCH SIZE!! 
        numsinds = length(sinds);
        numbatches = ceil((numdir*numsinds*numrot*numst + numdir*numrot*numst)*4/1048576 / (maxmem - currmem - 550));
        if numbatches < 1
            numbatches = numsinds;
        elseif numbatches > 1
            numbatches = numbatches + 1; %%%% Changed to deal with memory
        end
        batchsize = ceil(numsinds / numbatches);      
        numbatches = ceil(numsinds / batchsize);
        for b = 1:numbatches
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
            ssds = inf(numdir,numcurrim,numrot,numst,'single');
            currprojnorms = projnormsc(:,ones(numcurrim,1)); 
            ic = imcoeffs(curriminds,:)';
            currlps = lps(:,curriminds);
            for t = 1:numst
                currt = currtrans(t);
                if currt < 1
                    continue;
                end
                currimnorms = imnorms(currt,curriminds);
                currimnorms = currimnorms(onesdir,:);
                for r = 1:numrot
                    % THE SSD CALCULATION
                    ssds(:,:,r,t) = currimnorms;
                    for s = 1:numstruct
                        currips = currprojcoeffs((s-1)*numdir+1:s*numdir,:)*(ips(:,:,r,currt)*ic);
                        ssds(:,:,r,t) = ssds(:,:,r,t) + currlps(s*onesdir,:).*(currprojnorms((s-1)*numdir+1:s*numdir,:) - currips);
                    end
                end
            end
            for i = 1:numcurrim
                currssds = squeeze(ssds(:,i,:,:));
                [SSDs(curriminds(i)),minind] = min(currssds(:));    
                [dind,rind,tind] = ind2sub([numdir,numrot,numst],minind);
                insertinds = (curriminds(i)-1)*numstruct+1:curriminds(i)*numstruct;
                projinds(insertinds) = inds(dind:numdir:end);
                rotinds(insertinds) = rind;
                transinds(insertinds) = currtrans(tind);
                SSDs(insertinds) = imnorms(transinds(insertinds),curriminds(i)) + projnorms(projinds(insertinds)) - projcoeffs(projinds(insertinds),:)*(ips(:,:,rind,currtrans(tind))*imcoeffs(curriminds(i),:)');
            end
            clear ssds currssds currprojnorms ic currimnorms
        end
    end
end