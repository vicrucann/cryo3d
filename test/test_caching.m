% test caching functions
% given chunks of data, save them to cache folder
% then read chunk by chunk
% return the matrix chunk

num_chunks = 20;
type = 'single';
carr = create_cached_array([1000, 300, 360, 150], 'cache', type, num_chunks);
for i = 1:2
    tic;
    chunk = rand(1000,300,18,150, type); % break in 3rd dimension
    fprintf('done random matrix allocation\n');
    toc;
    
    tic;
    write_cached_array_chunk(carr, chunk, i);
    fprintf('done writing to cache of chunk %i\n', i);
    toc;
    
    tic;
    chunk_r = read_cached_array_chunk(carr, i, [1000, 300, 18, 150]);
    fprintf('done reading chunk %i from a cache\n', i);
    toc;
    fprintf('matrix equality: %i\n\n', isequal(chunk,chunk_r));    
end


    
    
