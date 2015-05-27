function [projinds,rotinds,SSDs,transinds,scales] = ...
    distribute_SSDs_fbm(projnorms,projcoeffs,imcoeffs,ips,ctfinds,numim,numctf,numproj,numrot,searchtrans,imnorms,maxmem,...
    ipaddrs,login,ppath,varmat,sleeptime,resfold,printout)

[ncluster ~] = find(ipaddrs==' ');
ncluster = size(ncluster,2)+1;
remmat = 'comp_SSDs_fast_best_match';