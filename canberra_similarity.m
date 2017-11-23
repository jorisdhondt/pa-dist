function d = canberra_similarity(a,b)
d= 0;
for i=1:size(a,2)
    if (a(i) ~= 0 || b(i) ~= 0)
       d = d + (abs(a(i) - b(i)))/(abs(a(i)) + abs(b(i))); 
    end
end
d = d/i;
d = 1-d;