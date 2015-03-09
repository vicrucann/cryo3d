function passed = write_cached_array_chunk(cacharr, chunk, idx_chunk)

passed = 0;
fprintf('writing the cached file...');
tic;
fname = [num2str(idx_chunk) '.dat'];
fid = fopen([cacharr.path fname], 'Wb');
fwrite(fid, chunk, cacharr.type);
fclose(fid);
fprintf('done\n');
toc;
fprintf('\n');
passed = 1;