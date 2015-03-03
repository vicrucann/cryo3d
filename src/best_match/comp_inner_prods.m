% Function to be used in EM code
% Computes inner products between image basis and rotated/translated projection basis

% Function edited to save the ips variable to cache files - in order to
% avoid "out of memory" error when dealing with large dimensions
% ips_cache output variable contains the list of cache files
% by default the files are saved in the current directory in 'cache\' folder

function ips = comp_inner_prods(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans)
%function [ips_cache, num_chunks] = comp_inner_prods(projbasis,imbasis,rots,numprojcoeffs,numrot,numimcoeffs,numpixsqrt,numpix,trans,searchtrans,numtrans)

projbasis3d_g = gpuArray(single(reshape(projbasis,[numpixsqrt, numpixsqrt, numprojcoeffs])));
imbasis_g = gpuArray(imbasis)';

%if ~exist('cache')
%    mkdir('cache');
%end
%ips_cache = 'cache\';
%type = 'single';

if nargin == 8  % Only rotations
    
    ips_g = gpuArray.zeros(numimcoeffs,numprojcoeffs,numrot);
    for r = 1:numrot
        projbasisrot_g = imrotate(projbasis3d_g,rots(r),'bilinear','crop');
        ips_g(:,:,r) = imbasis_g * reshape(projbasisrot_g,[numpix,numprojcoeffs]);
    end    
    ips_g = permute(ips_g,[2 1 3]);
    ips_g = 2*ips_g;
    ips = gather(ips_g);
    
%    fname = ['0.dat'];
%    fid = fopen([ips_cache fname], 'Wb');
%    fwrite(fid, ips, type); % this could be improved if this chunk info was already used
%    fclose(fid);
    
else            % Rotations + translations
    
    % Initializations, and determine if need to compute inner products in
    % batches due to limited space on gpu
    ips = zeros(numprojcoeffs,numimcoeffs,numrot,numtrans,'single');
    
    %fprintf('percent completed: ');
    
    validtrans = unique(searchtrans(searchtrans > 0))';
    g = gpuDevice;
    neededmem = numimcoeffs*numprojcoeffs*numrot*numtrans*8 + numpix*numprojcoeffs*8;
    numbatches = ceil(neededmem / g.FreeMemory);
    if numbatches > 1
        numbatches = numbatches + 1;
    end
    batchsize = ceil(numrot / numbatches);
    
    fprintf('Number of GPU batches: %i\n', numbatches);
    fprintf('Memory size of one batch: %i\n', batchsize);
    
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
        else
           ips(:,:,batchsize*(b-1)+1:end,:) = gather(ips_g);
        end
        
%         batch = gather(ips_g);
% gb = 1.5;
% type = 'single';
% stype = 4;
% dt = numtrans;
% npro = numprojcoeffs;
% nimg = numimcoeffs;
% dr = floor(gb*1024^3/(stype*nimg*npro*dt));

%         fname = [num2str(b) '.dat'];
%         fid = fopen([ips_cache fname], 'Wb');

%         fwrite(fid, chunk, type); % this could be improved if this chunk info was already used
%         fclose(fid);
%         
%         perc = b/numbatches*100;
%         fprintf('%i ', perc);
    end
    %fprintf('/n');
    clear  pbrottrans_g g
    
end

clear ips_g projbasis3d_g imbasis_g projbasisrot_g
