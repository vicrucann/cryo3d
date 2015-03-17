% Part of caching function data structure
% Allows to avoid Matlab out of memory error by caching a large array into
% several files on hard disk and then reading the necessary chunks using
% memmapfile Matlab function
% Victoria Rudakova 2015, victoria.rudakova(at)yale.edu

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
    else
        reqmem = 8; % assume it's double otherwise
    end
    chunkmem = reqmem;
    for i = 1:length(size)
        reqmem = reqmem*size(i); % total size of variable in bytes
        if (i ~= idx_broken)
            chunkmem = chunkmem*size(i);
        else
            chunkmem = chunkmem*ceil(size(idx_broken)/num_chunks);
        end
    end
    
    reqmem = 1.5*(2*chunkmem+reqmem); % assume it's 20% more than required to allow for other side variables
    
    archstr = computer('arch');
    if (isequal(archstr(1:3), 'win')) % if it's windows
        user = memory;
        if (user.MaxPossibleArrayBytes > reqmem)
            fprintf('No caching will be used, there is enough memory \n');
            caching = 0;
        else
            warning('Not enough memory: caching will be used. Processing time will be slower. ');
            caching = 1;
        end
    elseif (isequal(archstr(1:5),'glnxa')) % if linux
        [r w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        %memsize = stats(1); % bytes
        freemem = stats(3) + stats(end); % bytes
        if (freemem > reqmem)
            fprintf('No caching will be used, there is enough memory \n');
            caching = 0;
        else
            warning('Not enough memory: caching will be used. Processing time will be slower. ');
            caching = 1;
        end
    else % mac?
        error('Unrecognized or unsupported architecture');
    end    
end

if (caching == 0)
    data = zeros(size, type);
else
    data = 0;
end

if (~strcmp(path_cache(end), '\') && ~strcmp(path_cache(end), '/'))
    path_cache = [path_cache '\/'];
end

cacharr = struct('dimensions', size, 'path', path_cache, 'type', type, 'nchunks', num_chunks, 'broken', idx_broken, ...
    'caching', caching, 'data', data, 'currchunk', 1);