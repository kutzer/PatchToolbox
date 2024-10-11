function cycles = findUndirectedGraphCycles(G)
% FINDUNDIRECTEDGRAPHCYCLES finds the set of non-overlapping cycles
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
N = size(G, 1);
visited = false(1,N);

cycles = {};
while ~all(visited,'all')
    s = find(~visited,1,'first');
    cycles{end+1} = dfsearch(G,s);
    visited(cycles{end}) = true;
end