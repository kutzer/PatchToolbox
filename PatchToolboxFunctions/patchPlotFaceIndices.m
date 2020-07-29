function txt = patchPlotFaceIndices(varargin)
% PATCHPLOTFACEINDICES adds text indicating the face index value located at
% the centroid of each face of the specified patch object.
%   txt = PATCHPLOTFACEINDICES(ptch) defines the patch as either patch 
%   object.
%
%   txt = PATCHPLOTFACEINDICES(ptch,idx) allows the user to specify the
%   index values.
%
%   txt = PATCHPLOTFACEINDICES(ptch,idx,offset) allows the user to specify
%   an offset between the patch face and text label. 
%
%   Inputs:
%         ptch - patch object or structured array containing the fields 
%                "Vertices" and "Faces".
%          idx - N-element array containing face index value(s) to display. 
%                If idx = [] the function will display face labels.
%       offset - offset between the face centroid and the center of the
%                text for the face label (along the outward pointing unit
%                normal). This value defaults to 0.
%
%   Outputs:
%       txt - Nx1 array containing text object handles for each face label.
%       
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh.
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
        offset = 0;
    case 2
        idx = varargin{2};
        offset = 0;
    case 3
        idx = varargin{2};
        offset = varargin{3};
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

%% Calculate face unit normals
N = patchFaceNormal(ptch,idx);

%% Combine centroids and normals using offset 
X_offset = X + N*offset;

%% Add text
for i = 1:size(X_offset,1)
    str = sprintf('f_{%d}',idx(i));
    txt(i,1) = text(X_offset(i,1),X_offset(i,2),X_offset(i,3),str,...
        'Parent',axsTMP,'Tag','Patch Face Index Labels',...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

%% Migrate text to parent and delete temporary figure
set(txt,'Parent',mom);
delete(figTMP);