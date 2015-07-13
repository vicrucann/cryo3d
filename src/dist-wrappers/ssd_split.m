function output = ssd_split( in )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
numrot = in.numrot;
ncluster = in.ncluster;
dimensions = in.dimensions;

numrot_ = ceil(numrot / ncluster);
numrot_l = numrot - numrot_*(ncluster-1);

for i=1:ncluster
    r_begin = (i-1)*numrot_ + 1;
    dims = dimensions;
    if i~=ncluster
        r_end = numrot_ * i;
        dims(in.broken)=numrot_;
    else
        r_end = numrot;
        dims(in.broken)=numrot_l;
    end
    %fprintf('Saving split data\n');
    %t_save = tic;
    save([in.path_vars in.vars int2str(i) '.mat'], 'r_begin', 'r_end', 'dims', 'in');
    %toc(t_save);
end
output = 1;
end

