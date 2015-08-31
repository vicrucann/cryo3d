function map = projection_map(p, q)
M = size(q,2);
map = ones(1,M);
for j = 1 : M
    k = find(p==q(j)); 
    if ~isempty(k)
        map(j) = k(1);
    end
end
end