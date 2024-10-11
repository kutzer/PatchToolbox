function [cycles,subCycles] = findGraphCycles(directedGraph, startNode)
% FINDGRAPHCYCLES finds closed cycles on a directed graph given a starting
% node.
%   cycles = findGraphCycles(directedGraph, startNode)
%
%   Input(s)
%       directedGraph - directed graph object defining directional
%                       connectivity between nodes
%           startNode - scalar integer specifying the starting node for
%                       finding cycles in the directed graph.
%
%          cycles - N-element cell array containing ordered node indices
%                   defining closed cycles in the directedGraph
%       subCycles - M-element cell array containing ordered node indices
%                   defining the shortest subset of cycles that cover all
%                   nodes
%
%   See also digraph
%
%   M. Kutzer, 13Sep2024, USNA

debugTXT = false;

%% Check input(s)
narginchk(2,2);

switch lower( class(directedGraph) )
    case 'digraph'
        % Good input
    otherwise
        error('Graph input must be specified as a directed graph.');
end

if numel(startNode) ~= 1
    error('Start node must be scalar.');
end

if startNode < 1 || startNode > numnodes(directedGraph)
    error('Start node must be a value between 1 and %d for the specified graph.',...
        numnodes(directedGraph));
end

%% Initialize variables
% Initialize all nodes as unvisited
visited = false(1, numnodes(directedGraph));
% Initialize an empty current path
path = [];
% Initialize an empty set of cycles found
cycles = {};

%% Define cycles
% NOTE: This uses a recursive depth-first search (DFS) function to find all
%        cycles starting and ending at startNode
cycles = dfsCycles(directedGraph, startNode, startNode, visited, path, cycles);

%% Post-process cycles
% Identify trivial cycles
% -> Node A to Node B to Node A
tfRemove = cellfun(@(x)(numel(x) == 3),cycles);

if debugTXT
    % Display all found cycles
    fprintf('\n');
    fprintf('Removing "trivial" cycles found:\n');
    for i = find(tfRemove)
        fprintf('cycles{%d} = [ ',i);
        fprintf('%d ',cycles{i});
        fprintf(']\n');
    end
    fprintf('\n');
end

% Remove trivial cycles
cycles(tfRemove) = [];

% Sort the array from longest to shortest
nNodesPerCycle = cellfun(@numel,cycles);
[~,idx] = sort(nNodesPerCycle);
cycles = cycles(idx);

if debugTXT
    % Display all found cycles
    fprintf('Closed cycles found:\n');
    for i = 1:numel(cycles)
        fprintf('cycles{%d} = [ ',i);
        fprintf('%d ',cycles{i});
        fprintf(']\n');
    end
    fprintf('\n');
end

%% Find shortest subset of cycles that cover all nodes
% Initialize binary indicating unique cycles
tfCycles = false(1,numel(cycles));

% Find all nodes contained
nodesALL = [];
for i = 1:numel(cycles)
    nodesALL = unique([nodesALL,cycles{i}]);
end

idxLast = numel(cycles)-1;
if numel(nodesALL) ~= numel( unique(cycles{end}))
    fprintf('***SPECIAL CASE\n');
    idxLast = numel(cycles);
end

if debugTXT
    fprintf('nodesALL: ')
    fprintf('%d ',nodesALL)
    fprintf('\n');
end

for k = 2:idxLast
    nodes = [];
    idxCombo = nchoosek(1:idxLast,k);
    for i = 1:numel(idxCombo)
        nodes = unique([ cycles{idxCombo(i,:)} ]);

        if debugTXT
            fprintf('idxCombo = [ ');
            fprintf('%d ',idxCombo(i,:));
            fprintf(']\n');
            fprintf('   nodes: ')
            fprintf('%d ',nodes)
            fprintf('\n');
        end

        if numel(nodes) == numel(nodesALL)
            tfCycles(idxCombo(i,:)) = true;
            break
        end
    end

    if any(tfCycles)
        break
    end
end

% Isolate subcycles
subCycles = cycles(tfCycles);

if debugTXT
    % Display all found cycles
    fprintf('\n');
    fprintf('Minimum set of closed, shortest cycles found:\n');
    for i = 1:numel(subCycles)
        fprintf('subCycles{%d} = [ ',i);
        fprintf('%d ',subCycles{i});
        fprintf(']\n');
    end
end

% END FUNCTION
end

% Internal functions (unique workspace)
function [cycles,visited] = dfsCycles(graph, currentNode, startNode, visited, path, cycles)
% DFSCYCLES performs a depth-first search to find all cycles starting and
% ending as a specified node

% Mark current node as visited
visited(currentNode) = true;
path = [path, currentNode];  % Add the current node to the path

% Get successors (outgoing neighbors) of the current node
neighbors = successors(graph, currentNode);

% Explore each neighbor
for i = 1:length(neighbors)
    nextNode = neighbors(i);
    if nextNode == startNode
        % A cycle is found if we reach the startNode again
        cycles{end+1} = [path, startNode];  % Append cycle to the list
    elseif ~visited(nextNode)
        % Recursively explore the neighbor if it hasn’t been visited
        [cycles,visited] = ...
            dfsCycles(graph, nextNode, startNode, visited, path, cycles);
    end
end

% Backtrack: mark the current node as unvisited and remove from the path
visited(currentNode) = false;
end
