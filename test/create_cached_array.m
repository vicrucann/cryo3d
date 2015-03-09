function cacharr = create_cached_array(size, path_cache, type, num_chunks)

if ~exist(path_cache)
    mkdir(path_cache);
end

cacharr = struct('dimensions', size, 'path', [path_cache '\'], 'type', type, 'nchunks', num_chunks);