function [cycleWeights] = cycle_weight_opt_cvx(J, consIds, eps2)
% Optimize the cycle weights to maximize the condition number 
A = full(J*J');
A = (A+A')/2;
[u,v] = eig(A);
v = diag(v)';
ids = find(v > 1e-5);
u = u(:, ids);
J = u'*J;
[dim, numC] = size(J);
%
cvx_begin
variable w(numC);
variable Hmin(dim, dim);
variable Hmax(dim, dim);
variable lmin;
variable lmax;
Hmin == semidefinite(dim);
Hmax == semidefinite(dim);
minimize lmax-lmin
subject to
Hmin == (J.*(ones(dim,1)*w'))*J' - lmin*eye(dim);
Hmax == lmax*eye(dim) - (J.*(ones(dim,1)*w'))*J';
sum(w) == 1;
w >= 0;
w(consIds) >= eps2;
cvx_end
%
cycleWeights = w;
