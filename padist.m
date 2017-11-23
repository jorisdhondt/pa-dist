function Y = padist(X, distance, amountTerms, abs)
%               
%PADIST Pairwise distance between observations.
%   Y = PADIST(X) returns a vector Y containing the mcknow distances
%   between each pair of observations in the M-by-N data matrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is an
%   (M*(M-1)/2)-by-1 vector, corresponding to the M*(M-1)/2 pairs of
%   observations in X.
%
%   Y = PADIST(X, DISTANCE, AMOUNTTERMS, ABS) computes Y.  Choices are:       
%
%       DISTANCE:
%           'euclidean' - Euclidean distance
%       
%           'cosine'    - One minus the cosine of the included angle
%                         between observations (treated as vectors)
%       AMOUNTTERMS:
%           XXX         - The amount of terms used
%       ABS:
%           'normal'    - Take the highest XXX terms
%           'absolute'  - Take the highest XXX terms, after taking absolute
%                         values

[m, n] = size(X);
Y = zeros(1,m*(m-1)/2); %row matrix

if m < 2
   % Degenerate case, just return an empty of the proper size.
   Y = zeros(1,0);
   return;
end

% Create (I,J) defining all pairs of points
p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

for i=1:m*(m-1)/2
    Y(1,i) = distcalc(X(I(i),:), X(J(i),:), distance, amountTerms, abs);
end

% ----------------------------------------------
function d = distcalc(XI, XJ, dist, amount, absolute)

%DISTCALC Perform distance calculation for PDIST.
switch absolute
    
    case 'absolute' %d <- XI XJ
                %d = 1 - sum(XI.*XJ,2);
                
                [sortedXI,Isort_i] = sort(abs(XI));
                [sortedXJ,Isort_j] = sort(abs(XJ));
                subsp_ids = union(Isort_i(length(Isort_i)-(amount-1):length(Isort_i)),Isort_j(length(Isort_j)-(amount-1):length(Isort_j)));
                XI = XI(subsp_ids);
                XJ = XJ(subsp_ids);                
    
    case 'normal' %d <- XI XJ
                %d = 1 - sum(XI.*XJ,2);
                
                [sortedXI,Isort_i] = sort(XI);
                [sortedXJ,Isort_j] = sort(XJ);
                subsp_ids = union(Isort_i(length(Isort_i)-(amount-1):length(Isort_i)),Isort_j(length(Isort_j)-(amount-1):length(Isort_j)));
                XI = XI(subsp_ids);
                XJ = XJ(subsp_ids);
end

switch dist
    case 'cosine'
        XnormI = sqrt(sum(XI.^2, 2));
        XnormJ = sqrt(sum(XJ.^2, 2));
        deler = XnormI*XnormJ;
        d = 1 - (sum(XI.*XJ,2)/deler);
    case 'euclidean'
        d = sqrt(sum((XI-XJ).^2,2));
end