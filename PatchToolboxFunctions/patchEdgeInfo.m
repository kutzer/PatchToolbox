function edgeInfo = patchEdgeInfo(ptch)
% PATCHEDGEINFO calculates the unique edges associated with a patch object
% and associates each edge index with the face(s) it defines.
%   edgeInfo = PATCHEDGEINFO(ptch) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object.
%
%       edgeInfo.Edges      - vertex index pairs defining each unique edge 
%       edgeInfo.Faces      - edge indices associated with each face
%       edgeInfo.Vertices   - edge indices that contain each vertex
%       edgeInfo.Directions - direction of each face's edge (-1 or 1)
%           [v1, v2] vertices are in the correct order  ( 1)
%           [v2, v1] vertices are in the opposite order (-1)
%
%   M. Kutzer, 01Jun2020, USNA

% Updates:
%   10Oct2024 - Updated to account for NaN values in face rows

%% Parse inputs
try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Define all edges
%{
% METHOD 1 - Assumes no NaN values witin faces
% Repeat the first vertex index 
f_wrap = [f, f(:,1)];
% Define vertex and face mapping
edges = [];
for i = 1:size(f,2)
    edges = [edges; [f_wrap(:,i),f_wrap(:,i+1)]];
end

%% Define faces associated with each edge
% Define face indices
n = size(f,1);          % Total number of faces 
f_idx = transpose(1:n); % Index of each face
faces = repmat(f_idx,size(f,2),1);
%}

% MethOD 2 - Allows for NaN values within faces
edges = [];
faces = [];
nEdgesPerFace = 0;      % Maximum edges per face
n = size(f,1);          % Total number of faces 
for i = 1:n
    tf = ~isnan( f(i,:) );
    faceEdges = f(i,tf);
    faceEdges(2,:) = faceEdges([2:end,1]);
    edges = [edges; faceEdges.'];
    faces = [faces; repmat(i,nnz(tf),1)];
    nEdgesPerFace = max([nEdgesPerFace,nnz(tf)]);
end

%% Find unique edges
[sEdges,sIdx] = sort(edges,2);
%[uEdges,idxUnique,idxMap] = unique(sEdges,'rows');
[uEdges,~,idxMap] = unique(sEdges,'rows');

%% Define edge "directions"
bin = sIdx(:,1) == 1;
directions(bin) = 1;   % [v1, v2] vertices are in the correct order
directions(~bin) = -1; % [v2, v1] verteces are in the opposite order

%% Map unique edges to face indices and assign directions to face's edges 
h = waitbar(0,'Mapping unique edges to face indices...','Name','patchEdgeInfo.m');
n = size(f,1);          % Total number of faces 
for face = 1:n
    bin = faces == face;
    eFaces(face,:) = ...
        [reshape(idxMap(bin)    ,1,[]), nan(1,nEdgesPerFace-nnz(bin),1)];
    fDirections(face,:) = ...
        [reshape(directions(bin),1,[]), nan(1,nEdgesPerFace-nnz(bin),1)];
    % TODO - make sure the vertices are appropriately ordered to preserve
    %        face normal.
    waitbar(face/n,h);
end
delete(h);

%% Map vertices to unique edge indices
h = waitbar(0,'Mapping unique edges to vertex indices...','Name','patchEdgeInfo.m');
m = size(v,1);          % Total number of vertices
for vertex = 1:m
    bin = any(uEdges == vertex,2);
    %fprintf('%6d ',find(bin'));
    %fprintf('\n');
    v2eIdx = reshape(find(bin),1,[]);
    nV2eIdx = numel(v2eIdx);
    if vertex > 1 && nV2eIdx > size(eVertices,2)
        eVertices(:,nV2eIdx) = 0;
    end
    eVertices(vertex,1:nV2eIdx) = reshape(find(bin),1,[]);
    waitbar(vertex/m,h);
end
delete(h);

% Replace zeros with nan
tf = eVertices == 0;
eVertices(tf) = nan;

%% Package output
edgeInfo.Edges = uEdges;
edgeInfo.Faces = eFaces;
edgeInfo.Vertices = eVertices;
edgeInfo.Directions = fDirections;
