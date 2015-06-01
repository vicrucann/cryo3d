% Function to be used in subspaceEM/bestmatch code

% Created by Nicha C. Dvornek, 08/2013

function [done,err,pind] = check_convergence(proj_last,proj_est,convtol,numproj)

done = 1;
err = -1;
pind = -1;
for j = 1:numproj
    if sum(sum(((proj_last(:,:,j)-proj_est(:,:,j))).^2))/sum(sum((proj_last(:,:,j)).^2)) > convtol
        done = 0;
        err = sum(sum(((proj_last(:,:,j)-proj_est(:,:,j))).^2))/sum(sum((proj_last(:,:,j)).^2));
        pind = j;
        break;
    end
end