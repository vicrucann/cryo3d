function chunk_x = read_cached_array(cacharr, indices)
% example of indices = [0, 0, 1, 1]; % zero stands for ':'

% check the indices are in the right dimension range
if (indices(3) > cacharr.dimensions(3) || indices(4) > cacharr.dimensions(4))
    fprintf('Exceeding index: %f out of %f and %f out of %f \n', indices(3), cacharr.dimensions(3), indices(4), cacharr.dimensions(4));
    fprintf('Other caching data: size(data) = [%i %i %i %i]', size(cacharr.data,1), size(cacharr.data, 2),...
        size(cacharr.data,3), size(cacharr.data,4));
    error('Index exceeds matrix dimensions.');
end

if (cacharr.caching == 1)
    nc = cacharr.nchunks; % number of chunks
    ix = cacharr.broken; % idx of broken dimension
    nx = cacharr.dimensions(ix); % total number of dimensions
    dx = ceil(nx/nc); % number of dimensions per chunk
    idx_chunk = ceil(indices(ix)/dx); % filename calculation
    idx_data = mod(indices(ix), dx); % current chunk data offset based on broken idx
    if (idx_data == 0)
        idx_data = dx;
    end
    if (idx_chunk > cacharr.currchunk && idx_chunk > 1)
        mm = memmapfile([cacharr.path num2str(idx_chunk) '.dat'], 'Format', cacharr.type);
        cacharr.currchunk = idx_chunk;
        if (idx_chunk == nc)
            dx = nx-(nc-1)*dx; % the last chunk might have different size in broken dimension
        end
        dims = [cacharr.dimensions(1), cacharr.dimensions(2), dx, cacharr.dimensions(4)];
        cacharr.data = reshape(mm.Data, dims);
    end
    
    chunk_x = cacharr.data(:,:,idx_data, indices(4)); % in order to create more general chunk reader, try to use eval function here
else
    chunk_x = cacharr.data(:,:,indices(3), indices(4));
end
