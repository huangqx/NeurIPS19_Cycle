function [cycles] = cycle_basis_generator(A, numcycles, Para_admm)
% This function computes a number of cycles from an undirected graph with 
% the target number of cycles
% Input arguments:
%   'A': A binary adjacency matrix that is symmetric and the non-zero
%        entries indicate the edges
%   'numcycles': the target number of edges, this number should be above
%                num_edges - num_vertices + 1
%   'Para_admm.using_cvx': if this flag is 1, then we use the cvx package
%                          to solve the induced optimization problem. If
%                          this flag is 0, then we use one of the two ADMM
%                          solvers.
%   'Para_admm.exact_solver': if this flag is 1, then we use the exact
%                             sovler to solve any induced linear systems.
%                             If this flag is 0, then we use the conjugate
%                             gradient solver to solve the induced linear
%                             system.
%   'Para_admm.mu_init': Used in the ADMM solver. The default value is 1.
%   'Para_admm.mu_rho': Used in the ADMM solver. The default value is 1.006
%   'Para_admm.num_iters': Used in the ADMM solver. The default value is
%                          600
%   'Para_admm.cg_eps': A parameter for the cg solver for linear systems. 
%                       1e-6 is the default value
%   'Para_admm.cg_iters': A parameter for the cg solver for linear systems.
%                         50 is the default value
%   Output argument:
%    'cycles': each cycle has three fields, i.e., weight that indicates the
%    importance of this cycle, 'cycle': an ordered list of vertices, and
%    'oriId', which indicates its index of the candidate cycles. 
numV = size(A, 1);
% Get the neighbors of each vertex
for id = 1 : numV
    neighborIds{id} = find(A(id,:) ~= 0);
end

% Compute candidate cycles via breadth-first spanning trees
parentIds = cell(1, numV);
depths = zeros(1, numV);
for id = 1 : numV
    [parentIds{id}, depths(id)] = bfs_search_tree(neighborIds, id);
end
%
[edgeMatrix, edges] = adj_mat_2_edge(A);
[min_depth, rootId] = min(depths);
%
% Compute candidate cycles
numE = size(edges, 2);
cand_cycles = cell(1, numV*(numE-numV+1));
constraintIds = [];
for vId = 1 : numV
    newcycles = cycles_from_tree(edges, parentIds{vId});
    startId = (numE-numV+1)*(vId-1) + 1;
    endId = (numE-numV+1)*vId;
    cand_cycles(startId:endId) = newcycles;
    if vId == rootId
        constraintIds = startId:endId;
    end
end

% Use the condition number to select cycles
% Remove null spaces
J = generate_jacobi_matrix(cand_cycles, edgeMatrix);
eps2 = 1/4/(numE-numV+1);
% Choose one of optimizers to optimize the cycle weights
if Para_admm.using_cvx == 1
    cycleWeights = cycle_weight_opt_cvx(J, constraintIds, eps2);
else
    cycleWeights = cycle_weight_opt_admm(J, constraintIds, eps2, Para_admm);
end
% Use the weights to guide cycle selection
[cycles] = cycle_sampling(cand_cycles, cycleWeights,...
    constraintIds, numcycles);

function [J] = generate_jacobi_matrix(cycles, edgeMatrix)
% Generate the matrix whose rows index through the edges and whose column
% index through the cycles of the underlying graph.
numE = max(max(edgeMatrix));
numC = length(cycles);
numE_total = 0;
for cId = 1 : numC
    numE_total = numE_total + length(cycles{cId});
end
rowsJ = zeros(1, numE_total);
colsJ = zeros(1, numE_total);
valsJ = zeros(1, numE_total);
off = 0;
for cId = 1 : numC
    cycle = cycles{cId};
    ids1 = cycle;
    n = length(cycle);
    ids2 = [cycle(n), cycle(1:(n-1))];
    for j = 1 : n
        id1 = ids1(j);
        id2 = ids2(j);
        off = off + 1;
        rowsJ(off) = edgeMatrix(id1, id2);
        colsJ(off) = cId;
        if id1 < id2
            valsJ(off) = 1;
        else
            valsJ(off) = -1;
        end
    end
end
J = sparse(rowsJ, colsJ, valsJ, numE, numC);

% Compute the path to the root
function [path2root] = path_2_root(parentIdvec, vid)
%
off = 1;
path2root(off) = vid;
while parentIdvec(vid) ~= -1
    off = off + 1;
    vid = parentIdvec(vid);
    path2root(off) = vid;
end

function [newcycles] = cycles_from_tree(edges, parentIdvec)
%
numE = size(edges, 2);
numV = length(parentIdvec);
numCycles = numE - numV + 1;
newcycles = cell(1, numCycles);
off = 0;
for eId = 1 : numE
    v1id = edges(1, eId);
    v2id = edges(2, eId);
    if parentIdvec(v1id) ~= v2id && parentIdvec(v2id) ~= v1id
        path1 = path_2_root(parentIdvec, v1id);
        path2 = path_2_root(parentIdvec, v2id);
        n1 = length(path1);
        n2 = length(path2);
        cycle = [path1, zeros(1, n2-1)];
        for i = 2 : n2
            cycle(n1 + (i-1)) = path2(n2+1-i);
        end
        off = off + 1;
        newcycles{off} = cycle;
    end
end

function [parentIds, depth] = bfs_search_tree(neighborIds, rootNodeId)
%
numV = length(neighborIds);
parentIds = zeros(1, numV);
parentIds(rootNodeId) = -1;
visitedlist = zeros(1, numV);
visitedlist(1) = rootNodeId;
fringe_start = 1;
fringe_end = 2;
depth = 0;
while fringe_start < fringe_end
    off = fringe_end;
    depth = depth + 1;
    for i = fringe_start:(fringe_end-1)
        vid = visitedlist(i);
        nids = neighborIds{vid};
        for j = 1 : length(nids)
            nid = nids(j);
            if parentIds(nid) == 0
                parentIds(nid) = vid;
                visitedlist(off) = nid;
                off = off + 1;
            end
        end
    end
    fringe_start = fringe_end;
    fringe_end = off;
end

% Compute edge matrix
function [edgeMatrix, edges] = adj_mat_2_edge(A)
%
numV = size(A, 1);
[rows, cols, vals] = find(A);
ids = find(rows < cols);
rows = rows(ids);
cols = cols(ids);
%
edges = [rows, cols]';
numE = size(edges, 2);
edgeMatrix = zeros(numV, numV);
for eId = 1 : numE
    edgeMatrix(edges(1, eId), edges(2, eId)) = eId;
    edgeMatrix(edges(2, eId), edges(1, eId)) = eId;
end
