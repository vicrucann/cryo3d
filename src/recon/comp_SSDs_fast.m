% Function to be used with EM code
% Computes exp(SSD) between each image and template

% Created by Nicha C. Dvornek, 11/2013
% Last modified 03/2015

function [projinds,rotinds,iminds,lpcs_vals,SSDs,transinds] = comp_SSDs_fast(projnorms,projcoeffs,imcoeffs,ips,sigma1,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem)

projinds = -ones(numim,1);
rotinds = zeros(numim,1);
iminds = zeros(numim,1);
lpcs_vals = zeros(numim,1);
SSDs = zeros(numim,1);
transinds = zeros(numim,1);
numst = size(searchtrans,1);
searchtrans = searchtrans';
pdfbot = 2*sigma1^2;
indstart = 1;
currmem = monitor_memory_whos;

for c = 1:numctf
    inds = c:numctf:numproj;
    numprojc = length(inds);
    onesprojc = ones(numprojc,1,'int8');
    projnormsc = projnorms(inds);
    currprojcoeffs = projcoeffs(inds,:);
    cis = find(ctfinds == c)';
    searchtransu = unique(searchtrans(cis,:),'rows');
    numstu = size(searchtransu,1);
    for s = 1:numstu
        currtrans = searchtransu(s,:);
        sinds = find(ismember(searchtrans(cis,:),currtrans,'rows'));
        % DETERMINE BATCH SIZE
        numsinds = length(sinds);
        numbatches = ceil((numprojc*numsinds*numrot*numst + 2*numprojc*numrot*numst)*4/1048576 / (maxmem - currmem - 550));
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
            ssds = inf(numprojc,numcurrim,numrot,numst,'single');
            currprojnorms = projnormsc(:,ones(numcurrim,1)); 
            ic = imcoeffs(curriminds,:)';
            for t = 1:numst
                currt = currtrans(t);
                if currt < 1
                    continue;
                end
                currimnorms = imnorms(currt,curriminds);
                currimnorms = currimnorms(onesprojc,:);
                for r = 1:numrot
                    ssds(:,:,r,t) = currimnorms + currprojnorms - currprojcoeffs*(ips(:,:,r,currt)*ic);
                end
            end
            minssds = min(min(min(ssds,[],1),[],3),[],4);
            for i = 1:numcurrim
                currssds = squeeze(ssds(:,i,:,:));
                [pind,rtind,lval] = find(exp(-(currssds - minssds(i)).^2./pdfbot));
                indend = indstart+length(pind)-1;
                iminds(indstart:indend) = curriminds(i);
                projinds(indstart:indend) = inds(pind);
                [rind,tind] = ind2sub([numrot,numst],rtind);
                rotinds(indstart:indend) = rind;
                transinds(indstart:indend) = currtrans(tind);
                lpcs_vals(indstart:indend) = lval;
                ssdinds = sub2ind([numprojc,numrot,numst],pind,rind,tind);
                SSDs(indstart:indend) = currssds(ssdinds);
                indstart = indend + 1;
            end
            clear ssds currssds currprojnorms ic currimnorms
        end    
    end
end

outinds = projinds == -1;
projinds(outinds) = [];
rotinds(outinds) = [];
iminds(outinds) = [];
lpcs_vals(outinds) = [];
SSDs(outinds) = [];
transinds(outinds) = [];
