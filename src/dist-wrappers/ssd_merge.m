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

minindices = zeros(ncluster,numcurrim);
minvalues = zeros(ncluster,numcurrim);
for i=1:ncluster
    load([in.path_res '/' 'result_' in.vars int2str(i) '.mat']);
    minindices(i,:) = minidc;
    minvalues(i,:) = minval;
end
out = struct('minindices', minindices, 'minvalues', minvalues);

%  for i = 1:numcurrim
%      [val, ind] = min(minvalues(:,i));
%      SSDs(curriminds(i)) = val;
%      minind = minindices(ind,i);
%      [pind(i),rind(i),tind(i)] = ind2sub([numprojc,numrot,numst],minind);
%      projinds(curriminds(i)) = in.inds(pind(i));
%      rotinds(curriminds(i)) = rind(i);
%      transinds(curriminds(i)) = in.currtrans(tind(i));
%  end
%  out = struct('pind', pind, 'rind', rind, 'tind', tind, ...
%      'projinds', projinds, 'rotinds', rotinds, 'transinds', transinds, 'SSDs', SSDs);
end

