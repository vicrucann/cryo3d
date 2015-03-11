function cacharr = create_cached_array(size, path_cache, type, num_chunks, idx_broken, caching)

if ~exist(path_cache)
    mkdir(path_cache);
end

if (caching == 0)
    data = zeros(size, type);
else
    data = 0;
end

cacharr = struct('dimensions', size, 'path', [path_cache '\/'], 'type', type, 'nchunks', num_chunks, 'broken', idx_broken, ...
    'caching', caching, 'data', data, 'currchunk', 1);