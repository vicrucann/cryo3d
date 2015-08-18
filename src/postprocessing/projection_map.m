function map = projection_map(p, q)
map = zeros(size(q));
for j = 1 : M
    k = (p==q(j));
    map(k) = j;
end
end