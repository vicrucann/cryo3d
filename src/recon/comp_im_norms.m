% Function to be used with subspaceEM/bestmatch code

% Created by Nicha C. Dvornek, 09/2013
% Last modified 03/2015

function [imnorms] = comp_im_norms(imbasis,imcoeffs,maskim,trans,maxmem)

numim = size(imcoeffs,1);
if nargin == 4
    batchsize = double(ceil(numim / 100));
    numbatches = ceil(numim / batchsize);
else
    currmem = monitor_memory_whos;
    numpix = size(imbasis,1);
    numbatches = ceil(numim*numpix*8/1048576 / (maxmem - currmem - 550 ));
    if numbatches < 1
        numbatches = numim;
    end
    batchsize = ceil(numim/numbatches);
end
numtrans = size(trans,1);
sizeim = size(maskim);
imnorms = zeros(numtrans,numim);

for b = 1:numbatches
    if b == numbatches
        inds = (b-1)*batchsize+1:numim;
    else
        inds = (b-1)*batchsize+1:b*batchsize;
    end        
    numbinds = length(inds);
    onesbinds = ones(numbinds,1);
    approxims = imbasis*imcoeffs(inds,:)';
    for t = 1:numtrans
        dx = trans(t,1);
        dy = trans(t,2);
        maskcurr = zeros(sizeim);
        if dy < 0
            if dx < 0
                maskcurr(1:end+dy,1:end+dx) = maskim(1-dy:end,1-dx:end);
            else
                maskcurr(1:end+dy,1+dx:end) = maskim(1-dy:end,1:end-dx);
            end
        else
            if dx < 0
                maskcurr(1+dy:end,1:end+dx) = maskim(1:end-dy,1-dx:end);
            else
                maskcurr(1+dy:end,1+dx:end) = maskim(1:end-dy,1:end-dx);
            end
        end
        maskcurrcol = maskcurr(:);
        imcurrs = approxims .* maskcurrcol(:,onesbinds);
        imnorms(t,inds) = dot(imcurrs,imcurrs);
    end
end
