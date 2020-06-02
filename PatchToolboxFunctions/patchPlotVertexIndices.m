function txt = patchPlotVertexIndices(varargin)
% PATCHPLOTVERTEXINDICES adds text indicating the vertex index value 
% located at the position of the vertex of the specified patch object.
%   txt = PATCHPLOTVERTEXINDICES(ptch) defines the patch as either patch 
%   object.
%
%   txt = PATCHPLOTVERTEXINDICES(ptch,idx) allows the user to specify the
%   index values.
%
%   M. Kutzer, 02Jun2020, USNA

%% Check input(s)
narginchk(1,2);
ptch = varargin{1};

try
    v = ptch.Vertices;
    if size(v,2) < 3
        % Make vertices 3D
        v(:,3) = 0;
    end
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

switch nargin
    case 1
        idx = 1:size(v,1);
    case 2
        idx = varargin{2};
    otherwise
        error('Unexpected number of inputs.');
end

%% Get parent of patch object 
mom = get(ptch,'Parent');

%% Create temporary axes
figTMP = figure('Visible','off');
axsTMP = axes('Parent',figTMP);
hold(axsTMP,'on');

%% Isolate vertices
X = v(idx,:);

%% Add text
for i = 1:size(X,1)
    str = sprintf('v_{%d}',idx(i));
    txt(i,1) = text(X(i,1),X(i,2),X(i,3),str,'Parent',axsTMP,...
        'Tag','Patch Vertex Index Labels');
end

%% Migrate text to parent and delete temporary figure
set(txt,'Parent',mom);
delete(figTMP);