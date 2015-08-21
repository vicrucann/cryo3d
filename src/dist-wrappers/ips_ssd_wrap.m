function output = ips_ssd_wrap(file_mat, res_fname, ~, ~)

%% Inner product calcualtion

t_ips = tic;

fprintf('The mat file provided: %s\n', file_mat);
load(file_mat);

projbasis3d_g = gpuArray(single(reshape(in.projbasis,[in.numpixsqrt, in.numpixsqrt, in.numprojcoeffs])));
imbasis_g = gpuArray(in.imbasis)';
% Rotations + translations

% Initializations, and determine if need to compute inner products in
% batches due to limited space on gpu

validtrans = unique(in.searchtrans(in.searchtrans > 0))';
g = gpuDevice;

% if the distributer will be used, numbatches = c*numservers;
% where c is a coeff [1..], normally would be 1 and more for more
% requiring data sizes
neededmem = in.numimcoeffs*in.numprojcoeffs*in.numrot*in.numtrans*8 + in.numpix*in.numprojcoeffs*8;
numbatches = ceil(neededmem / g.FreeMemory);
if numbatches > 1
    numbatches = numbatches + 1;
end
batchsize = ceil(in.numrot / numbatches);

fprintf('Number of GPU batches: %i\n', numbatches);
fprintf('Total memory size of all batches in Gb, less than: %i\n', floor(neededmem/1024^3));

ips = CachedNDArray([in.numprojcoeffs,in.numimcoeffs,in.numrot,in.numtrans], 3, ...
    'type', 'single', 'var_name ', 'ips', 'work_path', in.pathcache,...
    'nchunks', numbatches, 'fcaching', in.caching, 'fdiscreet', 1);
%ips_test = zeros(numprojcoeffs,numimcoeffs,numrot,numtrans,'single');
%ips = Cacharr([numprojcoeffs,numimcoeffs,numrot,numtrans],...
%    pathcache, 'single', numbatches, 3, caching, 'ips');
%ips = zeros(in.numprojcoeffs,in.numimcoeffs,in.numrot,in.numtrans,'single');

fprintf('Percent completed: ');

for b = 1:numbatches
    
    % Setup for the current batch of rotations
    if b < numbatches
        currrots = in.rots(batchsize*(b-1)+1:batchsize*b);
    else
        currrots = in.rots(batchsize*(b-1)+1:end);
    end
    currnumrot = length(currrots);
    ips_g = gpuArray.zeros(in.numimcoeffs,in.numprojcoeffs,currnumrot,in.numtrans,'single');
    
    for r = 1:currnumrot
        
        % Rotate projection bases
        projbasisrot_g = imrotate(projbasis3d_g,currrots(r),'bilinear','crop');
        
        for t = validtrans
            
            % Translate rotated projection bases
            dx = in.trans(t,1);
            dy = in.trans(t,2);
            pbrottrans_g = gpuArray.zeros(in.numpixsqrt,in.numpixsqrt,in.numprojcoeffs,'single');
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
            ips_g(:,:,r,t) = imbasis_g * reshape(pbrottrans_g,[in.numpix,in.numprojcoeffs]);
            
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
    %clear chunk;
    
    perc = round(b/numbatches*100);
    fprintf('%u ', perc);
end
fprintf('\n');
clear  pbrottrans_g g


%save(res_fname, 'ips');
fprintf('ips calcualtion done\n');
toc(t_ips);

%% SSD calculation
t_ssd = tic;
% Initializations
projinds = -ones(in.numim,1);
rotinds = zeros(in.numim,1);
SSDs = zeros(in.numim,1);
transinds = zeros(in.numim,1);
scales = ones(in.numim,1);
numst = size(in.searchtrans,1);
searchtrans = in.searchtrans';
currmem = 0; %monitor_memory_whos;
minscale = 1.0;
maxscale = 1.0;

maxmem = 38000;

for c = 1:in.numctf
    
    % Get the relevant indices and norms for the current CTF
    inds = c:in.numctf:in.numproj;
    numprojc = length(inds);
    onesprojc = ones(numprojc,1,'int8');
    projnormsc = in.projnorms(inds);
    currprojcoeffs = in.projcoeffs(inds,:);
    cis = find(in.ctfinds == c)';
    %fprintf('size of ctfinds: [%d %d]\n', size(in.ctfinds,1), size(in.ctfinds,2));
    %fprintf('size of cis: [%d %d]\n', size(cis,1), size(cis,2)); %debug
    %fprintf('max and min: [%d %d]\n', max(cis), min(cis));
    fprintf('size of searchtrans: [%d %d]\n', size(searchtrans,1), size(searchtrans,2));
    
    searchtransu = unique(searchtrans(cis,:),'rows');
    numstu = size(searchtransu,1);
   
    % For each set of translations
    fprintf('\nNumber of set of translations: %i\n', numstu);
    for st = 1:numstu
        currtrans = searchtransu(st,:);
        sinds = find(ismember(searchtrans(cis,:),currtrans,'rows'));
        % Determine number of images to process per batch for handling limited memory
        numsinds = length(sinds);
        numbatches = ceil((numprojc*numsinds*in.numrot*numst + numprojc*in.numrot*numst)*4/1048576 / (maxmem - currmem - 550));
        if numbatches < 1
            numbatches = numsinds;
        elseif numbatches > 1
            numbatches = numbatches + 1; %%%% Changed to deal with memory
        end
        batchsize = ceil(numsinds / numbatches);
        numbatches = ceil(numsinds / batchsize);
        %fprintf('\nNumber of images per batch to process: %i\n', numbatches);
        
        % For each batch of images
        for b = 1:numbatches
            
            % Get the indices of the images
            if b == numbatches
                curriminds = cis(sinds((b-1)*batchsize+1:end));
            else
                if (b*batchsize > length(sinds))
                    ind = b*batchsize
                    length(sinds)
                end
                curriminds = cis(sinds((b-1)*batchsize+1:b*batchsize));
            end
            numcurrim = length(curriminds);
            
            % Setup ssds matrix, rep proj norms, and get current imcoeffs
            currprojnorms = projnormsc(:,ones(numcurrim,1));
            ic = in.imcoeffs(curriminds,:)';
            
            sum_currips = zeros(in.numrot, numst); % debug
            ssds = inf(numprojc,numcurrim,in.numrot,numst,'single');
            %ssds = CachedNDArray([numprojc,numcurrim,in.numrot,numst],...
            %    'single', 3, 'ssds', in.pathcache, 4, -1, 1, inf);
            
            % For each rotation
            %t_loop1 = tic;
            for r = 1:in.numrot
                % For each translation in the current search set
                for t = 1:numst
                    % Check if translation exists
                    currt = currtrans(t);
                    if currt < 1
                        continue;
                    end
                    
                    % Set up the current image norms
                    currimnorms = in.imnorms(currt,curriminds);
                    currimnorms = currimnorms(onesprojc,:);
                    
                    % First calculate the inner products between
                    % projections and current images
                    %currips = currprojcoeffs*(ips.read_cached_array([0,0,r,currt])*ic);
                    currips = currprojcoeffs*(ips(:,:,r,currt)*ic);
                    sum_currips(r,t) = sum(sum(currips)); % debug, control sum!
                    
                    %                   % Calculate scale and adjust
                    s = currips ./ currprojnorms / 2;
                    s(s < minscale) = minscale;
                    s(s > maxscale) = maxscale;
                    
                    % Calculate the ssds between each projection and image
                    ssds(:,:,r,t) = currimnorms + s.^2.*currprojnorms - s.*currips;
                end
            end
            %fprintf('loop1 is done\n');
            %toc(t_loop1);
            
            
            % For each image in the batch
            %t_loop2 = tic;
            pind = zeros(1,numcurrim);
            rind = zeros(1,numcurrim);
            tind = zeros(1,numcurrim);
            for i = 1:numcurrim
                % Find the min ssd and get the indices of the parameters of
                % the best batch
                currssds = squeeze(ssds(:,i,:,:));
                [SSDs(curriminds(i)),minind] = min(currssds(:));
                [pind(i),rind(i),tind(i)] = ind2sub([numprojc,in.numrot,numst],minind);
                projinds(curriminds(i)) = inds(pind(i));
                rotinds(curriminds(i)) = rind(i);
                transinds(curriminds(i)) = currtrans(tind(i));
            end
            clear ssds currssds currprojnorms ic currimnorms
        end
        progress_bar(st, numstu);
    end
   fprintf('\n');
end
fprintf('SSD calculation is done\n');
toc(t_ssd);
t_save = tic;
save(res_fname, 'projinds', 'rotinds', 'SSDs', 'transinds', 'scales');
fprintf('Wrapped data is saved\n');
toc(t_save);
output = 0;
end