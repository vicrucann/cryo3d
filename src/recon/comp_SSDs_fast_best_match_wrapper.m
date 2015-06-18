function comp_SSDs_fast_best_match_wrapper( fname, resfname )
%   To use as a part of cryo3d together with rshell-mat -
%   to parallelize ssds heavy computations
%   2015 Victoria Rudakova, vicrucann@gmail.com

fprintf('The data file provided: %s\n', fname);
load(fname);

ssdi = inf(numprojc, numcurrim, r_end - r_begin, numst,'single');

for r = r_begin:r_end
    for t = 1:numst
        % Check if translation exists
        currt = currtrans(t);
        if currt < 1
            continue;
        end
        % Set up the current image norms
        currimnorms = imnorms(currt,curriminds);
        currimnorms = currimnorms(onesprojc,:);
        
        % First calculate the inner products between
        % projections and current images
        currips = currprojcoeffs*(ipsi(:, :, r - r_begin + 1, currt) * ic);
        %currips = currprojcoeffs*(ips(:,:,r,currt)*ic);
        
        %                   % Calculate scale and adjust
        s = currips ./ currprojnorms / 2;
        s(s < minscale) = minscale;
        s(s > maxscale) = maxscale;
        % Calculate the ssds between each projection and image
        ssdi(:, :, r - r_begin + 1, t) = currimnorms + s.^2.*currprojnorms - s.*currips;
    end
end
save(resfname, 'ssdi');

end
