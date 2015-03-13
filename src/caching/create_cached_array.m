function cacharr = create_cached_array(size, path_cache, type, num_chunks, idx_broken, caching)

if ~exist(path_cache)
    mkdir(path_cache);
else
    delete([path_cache '\/*.dat']);
    warning('Cache folder has been cleared from previous cache data.');
end

if caching == -1 % caching is determined automatically
    if isequal(type, 'single')
        reqmem = 4; % bytes for single
    else % have to fulfill to add more data types and their sizes
        reqmem = 8;
    end
    for i = 1:length(size)
        reqmem = reqmem*size(i); % total size of variable in bytes
    end
    reqmem = 1.3*reqmem; % assume it's 30% more than required to allow for other side variables
    
    archstr = computer('arch');
    archstr = archstr(1:3);
    if (isequal(archstr, 'win')) % if it's windows
        user = memory;
        if (user.MaxPossibleArrayBytes > reqmem)
            fprintf('No caching will be used, there is enough memory \n');
            caching = 0;
        else
            warning('Not enough memory: caching will be used. Processing time will be slower. ');
            caching = 1;
        end
    else % if linux
        [r w] = unix('free | grep Mem');
        stats = str2double(regexpr(w, '[0-9]*', 'match'));
        memsize = stats(1); % bytes
        freemem = stats(3) + stats(end); % bytes
        if (freemem > reqmem)
            fprintf('No caching will be used, there is enough memory \n');
            caching = 0;
        else
            warning('Not enough memory: caching will be used. Processing time will be slower. ');
            caching = 1;
        end
    end    
end

if (caching == 0)
    data = zeros(size, type);
else
    data = 0;
end

cacharr = struct('dimensions', size, 'path', [path_cache '\/'], 'type', type, 'nchunks', num_chunks, 'broken', idx_broken, ...
    'caching', caching, 'data', data, 'currchunk', 1);