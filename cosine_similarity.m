function [s] = cosine_similarity(a,b)
%tmp = a.Values;
%tmp2 = (b.Values)';
s = a*b';