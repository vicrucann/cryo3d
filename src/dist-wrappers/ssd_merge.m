function out = ssd_merge( in )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

ncluster = in.ncluster;
numcurrim = in.numcurrim;
numprojc = in.numprojc;
numrot = in.numrot;
numst = in.numst;

minindices = zeros(ncluster,numcurrim);
minvalues = zeros(ncluster,numcurrim);
%fprintf('Loading merged data\n');
%t_load = tic;
for i=1:ncluster
    load([in.path_res '/' 'result_' in.vars int2str(i) '.mat']);
    minindices(i,:) = minidc;
    minvalues(i,:) = minval;
end
%toc(t_load);

minindices_glo = zeros(1,numcurrim);
minvalues_glo = zeros(1,numcurrim);
for i=1:numcurrim
    [val, ind_loc] = min(minvalues(:,i));
    minind_loc = minindices(ind_loc,i);
    rot_size = ceil(numrot/ncluster);
    if (ind_loc == ncluster)
        rot_size = numrot - (ncluster-1)*rot_size;
    end
    [p_loc, r_loc, t_loc] = ind2sub([numprojc, rot_size, numst], minind_loc);
    %idc_glo = p_loc + numprojc * (r_loc*ind_loc - 1) + ind_loc * numrot * r_loc * numprojc * (t_loc - 1);
    r_glo = r_loc + (ind_loc-1)*ceil(numrot/ncluster);
    assert(r_glo <= numrot, 'r_glo calculation failed: out of range');
    minindices_glo(i) = sub2ind([numprojc, numrot, numst], p_loc, r_glo, t_loc);
    minvalues_glo(i) = val;
end
out = struct('minindices', minindices_glo, 'minvalues', minvalues_glo);

end


