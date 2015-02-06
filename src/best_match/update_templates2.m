% Function to be used in EM code
% Updates projection templates from weighted averages of rotated/translated
% noisy images using latent probabilities

function proj_est = update_templates2(noisyims_g,iminds,projinds,rotinds,transinds,lprs_vals,rots,trans,numpixsqrt,numim,numproj)

proj_est = zeros(numpixsqrt,numpixsqrt,numproj);
g = gpuDevice;
neededmem = numpixsqrt*numpixsqrt*numim*8;
numbatches = ceil(neededmem / g.FreeMemory);
batchsize = ceil(numim / numbatches);
for b = 1:numbatches
    if b < numbatches
        binds = batchsize*(b-1)+1:batchsize*b;
    else
        binds = batchsize*(b-1)+1:numim;
    end
    numbinds = length(binds);
    validtrans = unique(transinds)';
    for t = validtrans
        dx = -trans(t,1);
        dy = -trans(t,2);
        transnoisyims_g = gpuArray.zeros(numpixsqrt,numpixsqrt,numbinds,'single');
        if dy < 0
            if dx < 0
                transnoisyims_g(1:end+dy,1:end+dx,:) = noisyims_g(1-dy:end,1-dx:end,binds);
            else
                transnoisyims_g(1:end+dy,1+dx:end,:) = noisyims_g(1-dy:end,1:end-dx,binds);
            end
        else
            if dx < 0
                transnoisyims_g(1+dy:end,1:end+dx,:) = noisyims_g(1:end-dy,1-dx:end,binds);
            else
                transnoisyims_g(1+dy:end,1+dx:end,:) = noisyims_g(1:end-dy,1:end-dx,binds);
            end
        end
        tinds = find(transinds == t & ismember(iminds,binds));
        rtinds = rotinds(tinds);
        rtindsu = unique(rtinds)';
        for r = rtindsu
            rinds = tinds(rtinds == r);
            [iindsu,~,tempimsinds] = unique(iminds(rinds));
            tempims = gather(imrotate(transnoisyims_g(:,:,iindsu-batchsize*(b-1)),-rots(r),'bilinear','crop'));
            for k = 1:length(rinds)
                proj_est(:,:,projinds(rinds(k))) = proj_est(:,:,projinds(rinds(k))) + tempims(:,:,tempimsinds(k))*lprs_vals(rinds(k));
            end
        end
    end
end

clear transnoisyims_g
