function ret = ismatrix(A)

s = size(A);
ret = length(s)==2 && all(s>=0);
