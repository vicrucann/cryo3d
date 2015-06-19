function output = ssd_wrap( file_mat, res_fname, cache_vname, ncache )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

fprintf('The mat file provided: %s\n', file_mat);
load(file_mat);

%assert(sum(file_dat~='0')>0, 'Error while checking file_dat: format is not correct\n');
%fprintf('The dat file provided: %s\n', file_dat);
%mm = memmapfile(file_dat, 'Format', in.ctype);
%ipsi = reshape(mm.Data, dims);
fprintf('Dat file is read\n');

ssdi = inf(in.numprojc, in.numcurrim, r_end - r_begin, in.numst,'single');

fprintf('The calculation loop for r in range [%i %i] and t in range [%i %i]\n', r_begin, r_end, 1, in.numst);
for r = r_begin:r_end
    for t = 1:in.numst
        % Check if translation exists
        currt = in.currtrans(t);
        if currt < 1
            continue;
        end
         % Set up the current image norms
        currimnorms = in.imnorms(currt,in.curriminds);
        currimnorms = currimnorms(in.onesprojc,:);
        % First calculate the inner products between
        % projections and current images
        r_idx = r - r_begin + 1;
        
        currips = in.currprojcoeffs*(ipsi(:, :, r_idx, currt) * in.ic);
        
        % Calculate scale and adjust
        s = currips ./ in.currprojnorms / 2;
        s(s < in.minscale) = in.minscale;
        s(s > in.maxscale) = in.maxscale;
        % Calculate the ssds between each projection and image
        ssdi(:, :, r - r_begin + 1, t) = currimnorms + s.^2.*in.currprojnorms - s.*currips;
    end
end
fprintf('loop-1 terminated\n');
minidc = zeros(1,in.numcurrim);
minval = zeros(1,in.numcurrim);
for i = 1:in.numcurrim
    currssdi = squeeze(ssdi(:,i,:,:));
    [minval(i), minidc(i)] = min(currssdi(:));
end
fprintf('loop-2 terminated\n');
save(res_fname, 'minval', 'minidc');
fprintf('output variable saved\n');
end