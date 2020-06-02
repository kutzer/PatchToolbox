function txt = patchPlotFaceIndices(varargin)
% PATCHPLOTFACEINDICES adds text indicating the face index value located at
% the centroid of each face of the specified patch object.
%   txt = PATCHPLOTFACEINDICES(ptch) defines the patch as either patch 
%   object.
%
%   txt = PATCHPLOTFACEINDICES(ptch,idx) allows the user to specify the
%   index values.
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   M. Kutzer, 02Jun2020, USNA

%% Check input(s)
narginchk(1,2);
ptch = varargin{1};

try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

switch nargin
    case 1
        idx = 1:size(f,1);
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

%% Calculate centroids
X = patchFaceCentroid(ptch,idx);

%% Add text
for i = 1:size(X,1)
    str = sprintf('f_{%d}',idx(i));
    txt(i,1) = text(X(i,1),X(i,2),X(i,3),str,'Parent',axsTMP,...
        'Tag','Patch Face Index Labels');
end

%% Migrate text to parent and delete temporary figure
set(txt,'Parent',mom);
delete(figTMP);