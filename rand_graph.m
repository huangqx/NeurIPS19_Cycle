function [A] = rand_graph(n,dmin,dmax)
%
A = zeros(n, n);
for i = 1 : n
    t = rand(1,1);
    d = floor(dmin*(1-t)+dmax*t);
    [s,order] = sort(rand(1,n));
    A(i, order(1:d)) = 1;
end
A = max(A,A');