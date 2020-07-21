This repository contains the code for generating a cycle-consistency basis from a undirected graph. Each cycle in the resulting cycle consistency basis basis has a weight.

Sample usage:
>> A = rand_graph(16,3,6); % This can be replaced by any graph you are interested
>> load('Para_admm.mat');
>> cycles = cycle_basis_generator(A, 64, Para_admm);

The code uses different solvers for solving the SDP problem. Recommendations: 
'small graphs' (Less than 50 edges): Para_admm.using_cvx = 1; 
'median graphs' (Less than 120 edges): Para_admm.using_cvx = 0; Para_admm.exact_solver = 1;
'large graphs' (Others): Para_admm.using_cvx = 0; Para_admm.exact_solver = 0;

