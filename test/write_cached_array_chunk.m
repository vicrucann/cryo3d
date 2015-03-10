function cacharr = write_cached_array_chunk(cacharr, chunk, idx_chunk)

if (cacharr.caching == 1)
    fname = [num2str(idx_chunk) '.dat'];
    fid = fopen([cacharr.path fname], 'Wb');
    fwrite(fid, chunk, cacharr.type);
    fclose(fid);
    if (idx_chunk ==  1)
        cacharr.data = chunk; % prepare data for the first use
    end
else
    batchsize = ceil(cacharr.dimensions(cacharr.broken) / cacharr.nchunks);
    if (idx_chunk < cacharr.nchunks)
        cacharr.data(:,:,batchsize*(idx_chunk-1)+1:batchsize*idx_chunk,:) = chunk;
    else
        cacharr.data(:,:,batchsize*(idx_chunk-1)+1:end,:) = chunk;
    end
end