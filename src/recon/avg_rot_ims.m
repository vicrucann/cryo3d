% Function to be used in subspaceEM/bestmatch code
% Computes weighted average of rotated images for each projection template

% Created by Nicha C. Dvornek, 08/2013
% Last modified 03/2015


function avgrotim = avg_rot_ims(imbasis,imcoeffs,projinds,rotinds,iminds,lpcs_vals,rots,numpix,numim,numproj,numrot,numimcoeffs,numpixsqrt,trans,transinds)

avgrotim_g = gpuArray.zeros(numpix,numproj,'single');
imbasis3d_g = gpuArray(single(reshape(imbasis,[numpixsqrt numpixsqrt numimcoeffs])));
lpcs_vals_g = gpuArray(lpcs_vals);

if nargin == 13     % Only rotations
    
    for r = 1:numrot
        rotimbasis_g = reshape(imrotate(imbasis3d_g,-rots(r),'bilinear','crop'),[numpix numimcoeffs]);
        rinds = find(rotinds == r);
        lpcs_r_g = gpuArray.zeros(numim,numproj);
        inds2d = sub2ind([numim numproj], iminds(rinds), projinds(rinds));
        lpcs_r_g(inds2d) = lpcs_vals_g(rinds);
        currprojindsu = unique(projinds(rinds));
        currimindsu = unique(iminds(rinds));
        lpcs_r_g = lpcs_r_g(currimindsu,:);
        lpcs_r_g = lpcs_r_g(:,currprojindsu);
        avgrotim_g(:,currprojindsu) = avgrotim_g(:,currprojindsu) + rotimbasis_g*(imcoeffs(currimindsu,:)'*lpcs_r_g);
    end
    
else                % Rotations and translations
    
    % For each translation
    validtrans = unique(transinds)';
    for t = validtrans
        
        % Translate the image basis
        dx = -trans(t,1);
        dy = -trans(t,2);
        transimbasis_g = gpuArray.zeros(numpixsqrt,numpixsqrt,numimcoeffs,'single');
        if dy < 0
            if dx < 0
                transimbasis_g(1:end+dy,1:end+dx,:) = imbasis3d_g(1-dy:end,1-dx:end,:);
            else
                transimbasis_g(1:end+dy,1+dx:end,:) = imbasis3d_g(1-dy:end,1:end-dx,:);
            end
        else
            if dx < 0
                transimbasis_g(1+dy:end,1:end+dx,:) = imbasis3d_g(1:end-dy,1-dx:end,:);
            else
                transimbasis_g(1+dy:end,1+dx:end,:) = imbasis3d_g(1:end-dy,1:end-dx,:);
            end
        end
        
        % Get the relevant rotations that go with this translation
        tinds = find(transinds == t);
        rtinds = rotinds(tinds);
        rtindsu = unique(rtinds)';
        
        % For each of the relevant rotations
        for r = rtindsu
            
            % Rotate the translated image basis
            rotimbasis_g = reshape(imrotate(transimbasis_g,-rots(r),'bilinear','crop'),[numpix numimcoeffs]);
            
            % Get the relevant image and projection indices
            rinds = tinds(rtinds == r);
            curriminds = iminds(rinds);
            currprojinds = projinds(rinds); 
            [currimindsu,~,ciui] = unique(curriminds);
            [currprojindsu,~,cpui] = unique(currprojinds);

            % Set up latent probability matrix
            lpcs_r_g = gpuArray.zeros(length(currimindsu),length(currprojindsu),'single');
            inds2d = sub2ind([length(currimindsu),length(currprojindsu)],ciui,cpui);
            lpcs_r_g(inds2d) = lpcs_vals_g(rinds);
            
            % Add on the relevant weighted, rotated images to the running
            % average image for each projection template
            temp_g = rotimbasis_g*(imcoeffs(currimindsu,:)'*lpcs_r_g);
            clear rotimbasis_g lpcs_r_g
            avgrotim_g(:,currprojindsu) = avgrotim_g(:,currprojindsu) + temp_g;
            clear temp_g;

        end
        clear transimbasis_g
    end
    
end

avgrotim = gather(avgrotim_g);
clear avgrotim_g imbasis3d_g rotimbasis_g lpcs_vals_g lpcs_r_g