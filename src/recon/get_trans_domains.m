% Function to determine the current best translation and to update
% translation search indices

% Created by Nicha C. Dvornek, 02/2014


function searchtrans = get_trans_domains(iminds,transinds,lpcs_vals,trans,transwidth,transdelta,numim)

if transwidth == 0
    transdelta = 1;
end
searchtrans = -ones((2*(floor(transwidth/transdelta))+1)*(2*(floor(transwidth/transdelta))+1),numim);

% Make search faster
f = find(lpcs_vals > 0.1);
iminds_search = iminds(f);
transinds_search = transinds(f);
lpcs_vals_search = lpcs_vals(f);

for i = 1:numim
    
    currinds = find(iminds_search == i);
    if length(currinds) > 0
        [~,bestind] = max(lpcs_vals_search(currinds));
        center = trans(transinds_search(currinds(bestind)),:);
    else
        currinds = find(iminds == i);
        [~,bestind] = max(lpcs_vals(currinds));
        center = trans(transinds(currinds(bestind)),:);
    end
    
    xinds = trans(:,1) >= center(1)-transwidth & trans(:,1) <= center(1)+transwidth;
    mc = mod(center(1),transdelta);
    takeoutinds = mod(trans(:,1),transdelta) ~= mc;
    xinds(takeoutinds) = 0;
    yinds = trans(:,2) >= center(2)-transwidth & trans(:,2) <= center(2)+transwidth;
    mc = mod(center(2),transdelta);
    takeoutinds = mod(trans(:,2),transdelta) ~= mc;
    yinds(takeoutinds) = 0;
    tinds = find(xinds & yinds);
    searchtrans(1:length(tinds),i) = tinds;
end