function d = manhattan_similarity(a,b)
d= 0;
for i=1:size(a,2)
    d = d + abs(a(i)-b(i));
end
d = 1-d/size(a,2);