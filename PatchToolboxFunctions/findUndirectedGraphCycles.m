function cycles = findUndirectedGraphCycles(G)
% FINDUNDIRECTEDGRAPHCYCLES finds the largest set of non-overlapping cycles
% in a directed graph. 
%   cycles = findUndirectedGraphCycles(G)
%
%   Input(s)
%       G - undirected graph object defining directional connectivity 
%           between nodes
%
%   Output(s)
%          cycles - N-element cell array containing ordered node indices
%                   defining closed cycles in the directedGraph
%
%   See also graph findDirectedGraphCycles
%
%   M. Kutzer, 13Sep2024, USNA

%% Check input(s)
narginchk(1,1);

switch lower( class(G) )
    case 'graph'
        % Good input
    otherwise
        error('Graph input must be specified as an undirected graph.');
end

%% Find cycles
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
