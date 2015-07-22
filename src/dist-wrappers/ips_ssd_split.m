function output = ips_ssd_split(in)
output = 0;
imcoeffs = in.imcoeffs;
ctfinds = in.ctfinds;
%numim = in.numim;
searchtrans = in.searchtrans;
imnorms = in.imnorms;

break_len = ceil(size(imcoeffs,1) / in.ncluster);
break_len_end = size(imcoeffs,1) - break_len*(in.ncluster-1);

for i=1:in.ncluster
    if i~=in.ncluster
        in.imcoeffs = imcoeffs((i-1)*break_len+1:i*break_len,:);
        in.ctfinds = ctfinds((i-1)*break_len+1:i*break_len,:);
        in.numim = break_len;
        in.searchtrans = searchtrans(:,(i-1)*break_len+1:i*break_len);
        in.imnorms = imnorms(:,(i-1)*break_len+1:i*break_len);
    else
        in.imcoeffs = imcoeffs((i-1)*break_len+1:end,:);
        in.ctfinds = ctfinds((i-1)*break_len+1:end,:);
        in.numim = break_len_end;
        in.searchtrans = searchtrans(:,(i-1)*break_len+1:end);
        in.imnorms = imnorms(:,(i-1)*break_len+1:end);
    end
    save([in.path_vars in.vars int2str(i) '.mat'], 'in');
end
end