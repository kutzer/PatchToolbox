function cycles = findUndirectedGraphCycles(G)
% Input: G is the adjacency matrix of an undirected graph
% Output: cycles is a cell array containing the largest set of non-overlapping cycles

numNodes = size(G, 1);
allCycles = {};
visited = false(1, numNodes);

% Find all cycles using DFS
for i = 1:numNodes
    if ~visited(i)
        path = [];
        allCycles = [allCycles; dfsFindCycles(G, i, visited, path, i)];
    end
end

% Sort cycles by length (descending)
cycleLens = cellfun(@length, allCycles);
[~, idx] = sort(cycleLens, 'descend');
allCycles = allCycles(idx);

% Select largest set of non-overlapping cycles
cycles = {};
usedNodes = false(1, numNodes);

for i = 1:length(allCycles)
    cycle = allCycles{i};
    if all(~usedNodes(cycle))  % Check if any node is already used
        cycles = [cycles; {cycle}];
        usedNodes(cycle) = true;
    end
end

end

function cycles = dfsFindCycles(G, node, visited, path, startNode)
% Helper function to find cycles using DFS
visited(node) = true;
path = [path, node];
cycles = {};

neighbors = find(G(node, :));
for i = 1:length(neighbors)
    neighbor = neighbors(i);
    if neighbor == startNode && length(path) > 2
        % Found a cycle
        cycles = [cycles; {path}];
    elseif ~visited(neighbor)
        new_cycles = dfsFindCycles(G, neighbor, visited, path, startNode);
        cycles = [cycles; new_cycles];
    end
end
end
