function out = ips_ssd_merge(in)

projinds_ = -ones(in.numim,1);
rotinds_ = zeros(in.numim,1);
SSDs_ = zeros(in.numim,1);
transinds_ = zeros(in.numim,1);
scales_ = ones(in.numim,1);

break_len = ceil(in.numim / in.ncluster);
%break_len_end = in.numim - break_len*(in.ncluster-1);

for i=1:in.ncluster
    load([in.path_res '/' 'result_' in.vars int2str(i) '.mat']);
    if (i~=in.ncluster)
        projinds_((i-1)*break_len+1:i*break_len) = projinds;
        rotinds_((i-1)*break_len+1:i*break_len) = rotinds;
        SSDs_((i-1)*break_len+1:i*break_len) = SSDs;
        transinds_((i-1)*break_len+1:i*break_len) = transinds;
        scales_((i-1)*break_len+1:i*break_len) = scales;
    else
        projinds_((i-1)*break_len+1:end) = projinds;
        rotinds_((i-1)*break_len+1:end) = rotinds;
        SSDs_((i-1)*break_len+1:end) = SSDs;
        transinds_((i-1)*break_len+1:end) = transinds;
        scales_((i-1)*break_len+1:end) = scales;
    end
end
out=struct('projinds', projinds_, 'rotinds', rotinds_, 'SSDs', SSDs_, ...
    'transinds', transinds_, 'scales', scales_);
end