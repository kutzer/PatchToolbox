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

debugTXT = false;

%% Check input(s)
narginchk(1,1);

switch lower( class(G) )
    case 'graph'
        % Good input
    otherwise
        error('Graph input must be specified as an undirected graph.');
end

%% Find cycles
N = numnodes(G);
visited = false(1,N);

if debugTXT
    nStr0 = sprintf('%d',N);
    nStr0 = numel(nStr0);
    str0 = repmat('-',1,nStr0);
    
    str = sprintf('%d',1:9);
    str = sprintf('%s ',str);
    str = repmat(str,1,ceil(N/10));
    str = str(1:N);
    fprintf('%s-%s\n',str0,str);
    
    fprintf(['%0',num2str(nStr0),'d-'],0);
    fprintf('%d',visited);
    fprintf('\n');
end

cycles = {};
while ~all(visited,'all')
    s = find(~visited,1,'first');
    cycles{end+1} = dfsearch(G,s);
    visited(cycles{end}) = true;

    if debugTXT
        fprintf(['%0',num2str(nStr0),'d-'],numel(cycles));
        fprintf('%d',visited);
        fprintf('\n');
    end
end