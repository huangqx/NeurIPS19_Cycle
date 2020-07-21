function [cycles] = cycle_basis_generator(A, numcycles, using_cvx)
% This function computes a number of cycles from an undirected graph with 
% the target number of cycles
numV = size(A, 1);
%
for id = 1 : numV
    neighborIds{id} = find(A(id,:) ~= 0);
end
%
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
%
J = generate_jacobi_matrix(cand_cycles, edgeMatrix);
eps2 = 1/4/(numE-numV+1);
if using_cvx
    cycleWeights = cycle_weight_opt_cvx(J, constraintIds, eps2);
else
    cycleWieghts = cycle_weight_opt_admm(J, constraintIds, eps2);
end
%
[cycles] = cycle_sampling(cand_cycles, cycleWeights,...
    constraintIds, numcycles);
%
% A = J*diag(cycleWeights)*J';
% ids = zeros(2,length(cycles));
% for id = 1 : length(cycles)
%     ids(1,id) = cycles{id}.oriId;
%     ids(2,id) = cycles{id}.weight;
% end
% A1 = J(:,ids(1,:))*diag(ids(2,:))*J(:,ids(1,:))';
% h = 10;

function [J] = generate_jacobi_matrix(cycles, edgeMatrix)
%
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
