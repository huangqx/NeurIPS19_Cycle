function [cycles] = cycle_sampling(cand_cycles, cycleweights,...
    consIds, numCycles)
%
for id = 1 : length(consIds)
    cId = consIds(id);
    cycles{id}.weight = cycleweights(cId);
    cycles{id}.cycle = cand_cycles{cId};
    cycles{id}.oriId = cId;
end
flags = zeros(1, length(cand_cycles));
flags(consIds) = 1;
remainingIds = find(flags == 0);
[sampleCycleIds, adjustedweights] = sub_sampling(cycleweights(remainingIds),...
    numCycles-length(consIds));
%
offset = length(consIds);
for id = 1 : length(sampleCycleIds)
    cId = remainingIds(sampleCycleIds(id));
    cycles{offset+id}.weight = adjustedweights(id);
    cycles{offset+id}.cycle = cand_cycles{cId};
    cycles{offset+id}.oriId = cId;
end

%
function [ids, newWeights] = sub_sampling(inputWeights, numCycles)
%
for p = 1:-0.01:0
    tp = power(inputWeights/max(inputWeights),p);
    if sum(tp) >= numCycles
        checkprobability = tp*numCycles/sum(tp);
        adjustedWeights = power(inputWeights, 1-p);
        ids = find(rand(1,length(tp)) < checkprobability');
        newWeights = adjustedWeights(ids)';
        break;
    end
end
newWeights = newWeights*sum(inputWeights)/sum(newWeights);