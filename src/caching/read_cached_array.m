% Part of caching function data structure
% Allows to avoid Matlab out of memory error by caching a large array into
% several files on hard disk and then reading the necessary chunks using
% memmapfile Matlab function
% Victoria Rudakova 2015, victoria.rudakova(at)yale.edu

function chunk_x = read_cached_array(cacharr, indices)
% example of indices = [0, 0, 1, 1]; % zero stands for ':'
% so it would be the same as cacharr(:, :, 1, 1);

% check indices have the right dimension
if (size(indices, 2) ~= size(cacharr.dimensions,2))
    fprintf('indices size: %i\n', size(indices, 2));
    error('Indices length is too large or too small, the read data might be not correct');
end

% check the indices are in the right dimension range
if (sum(indices(:) > cacharr.dimensions(:)) > 0)
%if (indices(3) > cacharr.dimensions(3) || indices(4) > cacharr.dimensions(4))
    fprintf('Exceeding indices:'); 
    fprintf('%i ', indices(:)); 
    fprintf('\nOutof: ');
    fprintf('%i ', cacharr.dimensions(:));
    fprintf('\n');
    fprintf('Other caching data: size(data) = ');
    fprintf('%i ', size(cacharr.data));
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
        %dims = [cacharr.dimensions(1), cacharr.dimensions(2), dx, cacharr.dimensions(4)];
        dims = cacharr.dimensions;
        dims(ix) = dx; % broken dimension has different size than original array
        cacharr.data = reshape(mm.Data, dims);
    end
    ind = indices;
    ind(ix) = idx_data;
    expr_ind = ind2str(ind);
    chunk_x = eval(['cacharr.data' expr_ind  ';']); % general, N-dimensional array
    %chunk_x = cacharr.data(:,:,idx_data, indices(4)); % in order to create more general chunk reader, try to use eval function here
else
    expr_ind = ind2str(indices);
    chunk_x = eval(['cacharr.data' expr_ind ';']);
    %chunk_x = cacharr.data(:,:,indices(3), indices(4));
end
end

% Converts array of indices of form [0 0 2 0] into string '(:,:,2,:)' for
% further eval usage
function expr = ind2str(indices)
expr = '(';
for i = 1:size(indices,2)
    if (indices(i) == 0) % take all elements - ':'    
        expr = strcat(expr, ':,');
    else % write down the index number
        expr = strcat(expr, [num2str(indices(i)) ',']);
    end
end
expr(end) = ')'; % get rid of the comma at the end and close the braket
end
