function [cycleWeights] = cycle_weight_opt_admm(J, consIds, eps2, Para_admm)
% Optimize the cycle weights to maximize the condition number 
J_proj = remove_null_space(J);
% Perform the ADMM optimizer to obtain the optimized weights
cycleWeights = admm_solver(J_proj, consIds, eps2, Para_admm);
%

function [w_opt] = admm_solver(J, consIds, eps2, Para)
%
[dim, dim_w] = size(J);
if dim > 100
    Para.exact_solver = 0;
end
w = ones(dim_w,1)/dim_w;
% primal variables
X_1 = (J.*(ones(dim,1)*w'))*J';
X_2 = X_1;
v = eig((X_1+X_1')/2);
s_1 = min(v);
s_2 = max(v);

% dual variables
y = 0;
Y_1 = zeros(dim, dim);
Y_2 = zeros(dim, dim);
z = zeros(dim_w,1);
mu = Para.mu_init;
delta = Para.mu_rho;

% Pre-compute the data structure for optimizing Y1 and Y2
dids = find(reshape(eye(dim), [1,dim*dim])); %also used later
J2 = kron(J, ones(dim,1)).*kron(ones(dim,1), J);
if Para.exact_solver == 1
    L = eye(dim*dim);
    
    L(dids,dids) = L(dids,dids) + ones(length(dids),length(dids));
    L = kron(eye(2),L);
    L = L + [J2;-J2]*[J2;-J2]';
    invL = inv(L);
end
%
for iter = 1: Para.num_iters
    % Optimize S1 and S2 via
    % S_1 = min_{S_1 \succeq 0} \|S_1 - Y_1 + mu*X_1\|_{F}^2
    % S_2 = min_{S_2 \succeq 0} \|S_2 - Y_2 + mu*X_2\|_{F}^2
    S_1 = sdp_cone_proj(Y_1 - mu*X_1);
    S_2 = sdp_cone_proj(Y_2 - mu*X_2);

    % Optimize Y1 and Y2 by minimizing the following objective function
    % min_{Y_1,Y_2}  \|1 - mu*s_1 - trace(Y_1)\|^2 
    %              + \|1 + mu*s_2 - trace(Y_2)\|^2
    %              + \|-S_1 + Y_1 + mu*X_1\|^2 
    %              + \|-S_2 + Y_2 + mu*X_2\|^2
    %              + \|diag(J^T(Y_2-Y_1)J) - z + y*1 + mu*w\|^2
    % This amounts to solve a linear system whose left-hand side is fixed
    fixed_vec_mul_uses = z - y*ones(dim_w,1) + mu*w;
    b1 = zeros(dim*dim, 1);
    b2 = b1;
    b1(dids) = b1(dids) + 1 + mu*s_1;
    b2(dids) = b2(dids) + 1 - mu*s_2;
    b1 = b1 + reshape(S_1 + mu*X_1, [dim*dim,1]) - J2*fixed_vec_mul_uses;
    b2 = b2 + reshape(S_2 + mu*X_2, [dim*dim,1]) + J2*fixed_vec_mul_uses;
    b = [b1;b2];
    if Para.exact_solver == 1
        x = invL*b;
    else
        x1_init = reshape(Y_1, [dim*dim,1]);
        x2_init = reshape(Y_2, [dim*dim,1]);
        x_init = [x1_init;x2_init];
        x = cgs(@(x)(afun(x, J2, dim, dids)), b, Para.cg_eps, Para.cg_iters, [],[], x_init);
    end
    Y_1 = reshape(x(1:(dim*dim)),[dim,dim]);
    Y_2 = reshape(x((dim*dim+1):(2*dim*dim)), [dim,dim]);

    % Optimize z via minimizing 
    % min_{z\geq 0} \|diag(J^T*(Y_2-Y_1)*J) - z + y*1 + mu*w -mu*eps2*1_{consids}\|^2
    fixed_vec_mul_uses2 = sum(J.*((Y_1-Y_2)*J))';
    z = -fixed_vec_mul_uses2 + y*ones(dim_w,1) - mu*w;
    z(consIds) = z(consIds) + eps2*mu;
    z = max(0, z);

    % Optimize y via minimizing
    % -2*mu*y + \|diag(J^T*(Y_2-Y_1)*J) - z + y*1 + mu*w\|^2
    y = (sum(z + fixed_vec_mul_uses2 + mu*w) - mu)/dim_w;
 
    % Update the primal variables 
    s_1 = s_1 + (1-trace(Y_1))/mu;
    s_2 = s_2 + (trace(Y_2)-1)/mu;
    X_1 = X_1 + (S_1-Y_1)/mu;
    X_2 = X_2 + (S_2-Y_2)/mu;
    w = w + (fixed_vec_mul_uses2 + z - y*ones(dim_w,1))/mu;
    
    % Update mu
    mu = mu*delta;
    %
    if mod(iter, 20) == 0
        fprintf('Finished iter = %d\n', iter);
    end
end
w_opt = w;

function [y] = afun(x, J2, dim, dids)
%
dim2 = dim*dim;
dim3 = 2*dim*dim;
x1 = x(1:dim2);
x2 = x((dim2+1):dim3);
y1 = x1;
y2 = x2;
y1(dids) = y1(dids) + sum(x1(dids));
y2(dids) = y2(dids) + sum(x2(dids));
%
tp = ((x1-x2)'*J2)';
tp2 = J2*tp;
%
y1 = y1 + tp2;
y2 = y2 - tp2;
%
y = [y1;y2];

function [X_pos] = sdp_cone_proj(X)
%
[u,v] = eig(X);
ids = find(diag(v)' > 0);
X_pos = u(:,ids)*v(ids,ids)*u(:,ids)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function removes the null space from the edge indictor space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J_proj] = remove_null_space(J)
A = full(J*J');
A = (A+A')/2;
[u,v] = eig(A);
v = diag(v)';
ids = find(v > 1e-5);
u = u(:, ids);
J_proj = u'*J;