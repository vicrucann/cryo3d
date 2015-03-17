% Part of caching function data structure
% Allows to avoid Matlab out of memory error by caching a large array into
% several files on hard disk and then reading the necessary chunks using
% memmapfile Matlab function
% Victoria Rudakova 2015, victoria.rudakova(at)yale.edu

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
    indices = zeros(1, size(cacharr.dimensions, 2)+1);
    indices(cacharr.broken) = batchsize*(idx_chunk-1)+1; % obtain vector of form [0 0 0 1 1 0] -> (:,:,:,1:1,:)
    if (idx_chunk < cacharr.nchunks)
        indices(cacharr.broken+1) = batchsize*idx_chunk;
    else
        indices(cacharr.broken+1) = size(cacharr.data, cacharr.broken);
    end
    expr_ind = ind2str(indices);
    eval(['cacharr.data' expr_ind '=chunk;']);
%     if (idx_chunk < cacharr.nchunks)
%        cacharr.data(:,batchsize*(idx_chunk-1)+1:batchsize*idx_chunk,:,:) = chunk;
%     else
%        cacharr.data(:,batchsize*(idx_chunk-1)+1:end,:,:) = chunk;
%     end
end
end

% Converts array of indices of form [0 1 5 0] into string '(:,1:5,:)' for
% further eval function
function expr = ind2str(indices)
expr = '(';
f = 0;
for i = 1:size(indices,2)
    if (indices(i) == 0) % take all elements - ':'    
        expr = strcat(expr, ':,');
    else % write down the index number
        if (~f)
            expr = strcat(expr, [num2str(indices(i)) ':']);
            f = 1;
        else
            expr = strcat(expr, [num2str(indices(i)) ',']);
        end
    end
end
expr(end) = ')'; % get rid of the comma at the end and close the braket
end

