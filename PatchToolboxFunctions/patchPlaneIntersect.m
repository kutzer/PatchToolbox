function iInfo = patchPlaneIntersect(ptch,abcd,ZERO)
% PATCHPLANEINTERSECT finds the intersection between a patch and a plane
%   tbd = patchPlaneIntersect(ptc,abcd)
%   ___ = patchPlaneIntersect(ptc,abcd,ZERO)
%
%   Input(s)
%       ptch - patch object or structured array with fields "Faces" and
%              "Vertices"
%       abcd - 4-element array defining the coefficients for a 3D plane
%       ZERO - [OPTIONAL] scalar positive value close to zero. Default
%              value is 1e-8.
%
%   Output(s)
%      tbd
%
%   M. Kutzer, 10Oct2024, USNA

%% Check input(s)
narginchk(2,3);

if nargin < 3
    ZERO = 1e-8;
end

if numel(abcd) ~= 4
    error('Plane must be defined using a 4-element array.');
end
abcd = reshape(abcd,1,[]);

% Parse inputs
try
    eInfo = patchEdgeInfo(ptch);
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Find edges that intersect with plane
iInfo.X = [];
iInfo.tfFaces = [];
for i = 1:size(eInfo.Edges,1)
    edgePts = ptch.Vertices(eInfo.Edges(i,:),:).';
    pnt = intersectPlaneSegment(abcd,edgePts);

    if ~isempty(pnt)
        iInfo.X(:,end+1) = pnt;
        iInfo.tfFaces(:,end+1) = any(eInfo.Faces == i,2);
    end
end
