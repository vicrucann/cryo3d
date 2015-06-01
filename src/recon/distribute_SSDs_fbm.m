function [projinds,rotinds,SSDs,transinds,scales] = ...
    distribute_SSDs_fbm(projnorms,projcoeffs,imcoeffs,ips,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem,...
    ipaddrs,login,ppath,varmat,sleeptime,resfold,printout)

[ncluster ~] = find(ipaddrs==' ');
ncluster = size(ncluster,2)+1;
remmat = 'comp_SSDs_fast_best_match_wrapper'; % if distributor is not used, an original function will be run from best_match
bashscript = fullfile(pwd,'dhead.sh');

% split data
fprintf('\n\nCalculation using remotes \n');
numim_ = ceil(numim / ncluster);
numim_l = numim - numim_*(ncluster-1); % size of last might be different
for i=1:ncluster
    if (i ~= ncluster)
        
    else
    end
    save([varmat int2str(i) '.mat'], 'projnorms_', 'projcoeffs_', 'imcoeffs_', 'ips_', 'ctfinds_', 'numim_',...
        'numctf_','numproj_','numrot_','searchtrans_','imnorms_','maxmem_');
end
system(['chmod u+x ' bashscript]);
if printout
    cmdStr = [bashscript ' ' login ' ' ppath ' ' ipaddrs ' ' remmat ' ' varmat ' ' int2str(sleeptime) ' ' resfold];
else
    cmdStr = [bashscript ' ' login ' ' ppath ' ' ipaddrs ' '...
        remmat ' ' varmat ' ' int2str(sleeptime) ' ' resfold '>' remmat '.log 2>&1'];
end
% perform the command
system(cmdStr);

% merge the results
projinds = -ones(numim,1);
rotinds = zeros(numim,1);
SSDs = zeros(numim,1);
transinds = zeros(numim,1);
scales = zeros(numim,1);
for i=1:ncluster
    load([resfold '/' 'result_' varmat int2str(i) '.mat']);
    if (i ~= ncluster)
        projinds(szx*(i-1)+1:szx*i) = projinds_;
        rotinds(szx*(i-1)+1:szx*i) = rotinds_;
        SSDs(szx*(i-1)+1:szx*i) = SSDs_;
        transinds(szx*(i-1)+1:szx*i) = transinds_;
        scales(szx*(i-1)+1:szx*i) = scales_;
    else
        projinds(szx*(i-1)+1:end) = projinds_;
        rotinds(szx*(i-1)+1:end) = rotinds_;
        transinds(szx*(i-1)+1:end) = transinds_;
        SSDs(szx*(i-1)+1:end) = SSDs_;
        scales(szx*(i-1)+1:end) = scales_;
    end
end