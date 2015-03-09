function progress_bar(cit, tot, pb, ps)

if (nargin == 2)
    pb = 0.25;
    ps = 0.05;
end

perc_b = ceil(pb*tot);
perc_s = ceil(ps*tot);

if (mod(cit, perc_b) == 0)
    fprintf('%i', perc_b);
elseif (mod(cit, perc_s) == 0)
    fprintf('.');
end