function map = projection_map(p, q)
M = size(q,2);
map = zeros(1,M);
for j = 1 : M
    k = find(p==q(j)); 
    if isempty(k)
        fprintf('\n');
    end
    %map(j) = k(1);
end
end