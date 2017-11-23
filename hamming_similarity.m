function c = hamming_similarity(a,b)
at = size(find(a),2);
bt = size(find(b),2);
common = 0;
for i = 1:size(a,2)
    if(a(i) == 1 && b(i) == 1)
        common = common +1;
    end
end

c = common/(at^(0.5)*bt^(0.5));