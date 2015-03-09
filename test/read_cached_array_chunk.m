function chunk = read_cached_array_chunk(cacharr, idx_chunk, dims)
% dims is vector of dimensions [dim1 x dim2 x dim3 x ...]

fprintf('reading from the cached file...');
tic;
mm = memmapfile([cacharr.path num2str(idx_chunk) '.dat'], 'Format', cacharr.type);
chunk = reshape(mm.Data, dims);
fprintf('done\n');
toc;
fprintf('\n');
