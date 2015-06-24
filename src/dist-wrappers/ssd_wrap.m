function output = ssd_wrap( file_mat, res_fname, cache_vname, ncache )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

fprintf('The mat file provided: %s\n', file_mat);
load(file_mat);

vol = in.volume(in.broken);
nvol = ceil(in.dimensions(in.broken) / vol);
vol_end = in.dimensions(in.broken) - (nvol-1)*vol;
dims = in.volume;

fprintf('The broken dimension of ips has size of %i or %i\n', vol, vol_end);
fprintf('dims=%i %i %i %i\n', dims(1), dims(2), dims(3), dims(4));
idxc = ceil(r_begin/vol);
idxc_end = ceil(r_end/vol);
fprintf('idxc=%i, idxc_end=%i\n', idxc, idxc_end);
file_dat = [cache_vname int2str(idxc) '.dat'];
mm = memmapfile(file_dat, 'Format', in.ctype);

if (idxc == idxc_end)
    vol=vol_end;
    dims(in.broken) = vol_end;
end

ipsi = reshape(mm.Data, dims);
clear mm;
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
        idxc_ = ceil(r/vol);
        if (idxc ~= idxc_) % memmapfile the next file
            idxc = idxc_;
            if (idxc == idxc_end)
                vol = vol_end;
                dims(in.broken) = vol_end;
            end
            fprintf('idxc=%i, r=%i, vol=%i\n', idxc, r, vol);
            file_dat = [cache_vname int2str(idxc) '.dat'];
            mm = memmapfile(file_dat, 'Format', in.ctype);
            ipsi = reshape(mm.Data, dims);
            fprintf('Dat file is read\n');
        end
        
        r_idx_loc = mod(r, vol);
        if (r_idx_loc == 0)
            r_idx_loc = vol;
        end
        
        currips = in.currprojcoeffs*(ipsi(:, :, r_idx_loc, currt) * in.ic);
        
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
