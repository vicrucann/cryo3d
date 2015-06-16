function out = ssd_merge( in )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

ncluster = in.ncluster;
numcurrim = in.numcurrim;
numim = in.numim;
curriminds = in.curriminds;
numprojc = in.numprojc;
numrot = in.numrot;
numst = in.numst;
currtrans = in.currtrans;

minindices = zeros(ncluster,numcurrim);
minvalues = zeros(ncluster,numcurrim);
for i=1:ncluster
    load([in.path_res '/' 'result_' in.vars int2str(i) '.mat']);
    minindices(i,:) = minidc;
    minvalues(i,:) = minval;
end
 % For each image in the batch
 pind = zeros(1,numcurrim);
 rind = zeros(1,numcurrim);
 tind = zeros(1,numcurrim);
 SSDs = zeros(numim,1);
 
 projinds = -ones(numim,1);
 rotinds = zeros(numim,1);
 transinds = zeros(numim,1);
 for i = 1:numcurrim
     [val, ind] = min(minvalues(:,i));
     SSDs(curriminds(i)) = val;
     minind = minindices(ind,i);
     [pind(i),rind(i),tind(i)] = ind2sub([numprojc,numrot,numst],minind);
     projinds(curriminds(i)) = inds(pind(i));
     rotinds(curriminds(i)) = rind(i);
     transinds(curriminds(i)) = currtrans(tind(i));
 end
 out = struct('pind', pind, 'rind', rind, 'tind', tind, ...
     'projinds', projinds, 'rotinds', rotinds, 'transinds', transinds, 'SSDs', SSDs);
end

