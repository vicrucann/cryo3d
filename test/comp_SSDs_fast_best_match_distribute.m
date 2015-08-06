function [projinds,rotinds,SSDs,transinds,scales] = comp_SSDs_fast_best_match_distribute(...
    projnorms,projcoeffs,imcoeffs,ips,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem,...
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

d = Distributor(login, path_rem, ipaddrs, path_vars, vars, path_curr, sleeptime, path_res, printout);
if (d.ncluster > 1)
    d.scp_cached_data(ips);
end

for c = 1:numctf
    
end

scales(scales < minscale) = minscale;
scales(scales > maxscale) = maxscale;
end