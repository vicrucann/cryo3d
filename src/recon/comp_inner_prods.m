% Function to be used in subspaceEM/bestmatch code
% Computes inner products between image basis and rotated/translated projection basis

% Created by Nicha C. Dvornek, 08/2013

% Last modified by Victoria Rudakova, 05/2015
% Function edited to save the ips variable to cache files - in order to
% avoid "out of memory" error when dealing with large dimensions
% ips_cache output variable contains the list of cache files
% by default the files are saved in the current directory in 'cache\' folder

function ips = comp_inner_prods(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans, caching, pathcache, ipaddrs)
%function ips_cache = comp_inner_prods(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans, caching)

projbasis3d_g = gpuArray(single(reshape(projbasis,[numpixsqrt, numpixsqrt, numprojcoeffs])));
imbasis_g = gpuArray(imbasis)';

if nargin == 8  % Only rotations ! the parameter number had changed, this block is obsolete
    
    ips_g = gpuArray.zeros(numimcoeffs,numprojcoeffs,numrot);
    for r = 1:numrot
        projbasisrot_g = imrotate(projbasis3d_g,rots(r),'bilinear','crop');
        ips_g(:,:,r) = imbasis_g * reshape(projbasisrot_g,[numpix,numprojcoeffs]);
    end    
    ips_g = permute(ips_g,[2 1 3]);
    ips_g = 2*ips_g;
    ips = gather(ips_g);
       
else            % Rotations + translations
    
    % Initializations, and determine if need to compute inner products in
    % batches due to limited space on gpu
        
    validtrans = unique(searchtrans(searchtrans > 0))';
    g = gpuDevice;
    
    % if the distributer will be used, numbatches = c*numservers; 
    % where c is a coeff [1..], normally would be 1 and more for more
    % requiring data sizes
    [ncluster ~] = find(ipaddrs==' ');
    ncluster = size(ncluster,2)+1;
    neededmem = numimcoeffs*numprojcoeffs*numrot*numtrans*8 + numpix*numprojcoeffs*8;
    if ncluster == 1    % no distributer used, just a local machine
        numbatches = ceil(neededmem / g.FreeMemory);
        if numbatches > 1
            numbatches = numbatches + 1;
        end
    else % in case if we use distributor
        numbatches = ncluster;
        caching=1;
        k=2;
        while (numbatches < ceil(neededmem / g.FreeMemory))
            numbatches = ncluster*k;
            k=k+1;
        end        
    end
    batchsize = ceil(numrot / numbatches);
    
    fprintf('Number of GPU batches: %i\n', numbatches);
    fprintf('Total memory size of all batches in Gb, less than: %i\n', floor(neededmem/1024^3));
    
    ips = CachedNDArray([numprojcoeffs,numimcoeffs,numrot,numtrans], 3, ...
    'type', 'single', 'var_name', 'ips', 'path_cache', pathcache,...
    'nchunks', numbatches, 'fcaching', caching, 'fdiscreet', 1);
    %ips_test = zeros(numprojcoeffs,numimcoeffs,numrot,numtrans,'single');
    %ips = Cacharr([numprojcoeffs,numimcoeffs,numrot,numtrans],...
    %    pathcache, 'single', numbatches, 3, caching, 'ips');
    %ips = zeros(numprojcoeffs,numimcoeffs,numrot,numtrans,'single');
    
    fprintf('Percent completed: ');
    
    for b = 1:numbatches
        
        % Setup for the current batch of rotations
        if b < numbatches
            currrots = rots(batchsize*(b-1)+1:batchsize*b);
        else
            currrots = rots(batchsize*(b-1)+1:end);
        end
        currnumrot = length(currrots);
        ips_g = gpuArray.zeros(numimcoeffs,numprojcoeffs,currnumrot,numtrans,'single');
        
        for r = 1:currnumrot
            
            % Rotate projection bases
            projbasisrot_g = imrotate(projbasis3d_g,currrots(r),'bilinear','crop');
            
            for t = validtrans
                
                % Translate rotated projection bases
                dx = trans(t,1);
                dy = trans(t,2);
                pbrottrans_g = gpuArray.zeros(numpixsqrt,numpixsqrt,numprojcoeffs,'single');
                if dy < 0
                    if dx < 0
                        pbrottrans_g(1:end+dy,1:end+dx,:) = projbasisrot_g(1-dy:end,1-dx:end,:);
                    else
                        pbrottrans_g(1:end+dy,1+dx:end,:) = projbasisrot_g(1-dy:end,1:end-dx,:);
                    end
                else
                    if dx < 0
                        pbrottrans_g(1+dy:end,1:end+dx,:) = projbasisrot_g(1:end-dy,1-dx:end,:);
                    else
                        pbrottrans_g(1+dy:end,1+dx:end,:) = projbasisrot_g(1:end-dy,1:end-dx,:);
                    end
                end
                
                % Calculate the inner products
                ips_g(:,:,r,t) = imbasis_g * reshape(pbrottrans_g,[numpix,numprojcoeffs]); 
                
            end
        end
        
        % Rearrange order and transfer back to host memory
        ips_g = permute(ips_g,[2 1 3 4]);
        ips_g = 2*ips_g;
        
        if b < numbatches
           ips(:,:,batchsize*(b-1)+1:batchsize*b,:) = gather(ips_g);
           %ips_test(:,:,batchsize*(b-1)+1:batchsize*b,:) = gather(ips_g);
        else
           ips(:,:,batchsize*(b-1)+1:end,:) = gather(ips_g);
           %ips_test(:,:,batchsize*(b-1)+1:end,:) = gather(ips_g);
        end
        
        %chunk = single(gather(ips_g));
        %ips.write_cached_array_chunk(chunk, b);
        clear chunk;
        
        perc = round(b/numbatches*100);
        fprintf('%u ', perc);
    end
        fprintf('\n');
    clear  pbrottrans_g g
    
end

clear ips_g projbasis3d_g imbasis_g projbasisrot_g
