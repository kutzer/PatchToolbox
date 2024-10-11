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

debug = true;

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
    % Edge info
    eInfo = patchEdgeInfo(ptch);
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end
% Patch vertices (3xM array)
Xv = ptch.Vertices.';

%% Debug plot
if debug
    fig = figure('Name','patchPlaneIntersect.m');
    axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
    view(axs,3);

    ptc = patch(ptch,'FaceColor','b','FaceAlpha',0.2,'EdgeColor','none');
    
    Xlim(1,:) = xlim(axs);
    Xlim(2,:) = ylim(axs);
    Xlim(3,:) = zlim(axs);

    s = norm( diff(Xlim,1,2) )/2;
    X = mean(Xlim,2);
    plotPlane(axs,abcd,X,s);
end

%% Find edges that intersect with plane
% Initialize variables
Xint = [];      % 3xN array of intersection points
tfFaces = [];   % FxN logical array representing faces adjacent to point
tfEdges = [];   % ExN logical array representing edges adjacent to point
adjXint = [];   % NxN logical array defining adjacency between intersects

nEdges = size(eInfo.Edges,1);
nFaces = size(eInfo.Faces,1);
for i = 1:nEdges
    % Extract end-points of Edge_i (3x2 array)
    edgePts = Xv(:,eInfo.Edges(i,:));

    % Define the point(s) of intersection with the slicing plane
    %   pnt - []  -> No intersect
    %   pnt - 3x1 -> Segment intersect plane
    %   pnt - 3x2 -> Segment on plane
    [pnt,tfEndPnt] = intersectPlaneSegment(abcd,edgePts,ZERO);
    
    % Process intersection by type
    % -> Intersection(s) are the edge end-point(s)
    if any(tfEndPnt)

        % Point(s) of intersection is/are end-points of the edge
        jj = 0;
        for j = 1:numel(tfEndPnt)

            % Skip point if it is not on the plane
            if ~tfEndPnt(j)
                continue;
            end

            % Append intersection point
            jj = jj+1;
            Xint(:,end+1) = pnt(:,jj);
            
            % Define vertex index
            idxVert = eInfo.Edges(i,j);
            
            % Define adjacent edge indices
            idxEdge = eInfo.Vertices(idxVert,:);
            % -> Remove NaN values
            idxEdge = idxEdge(~isnan(idxEdge));

            % Update tfEdges
            tfEdges(:,end+1) = false(nEdges,1);
            tfEdges(idxEdge,end) = true;

            % Define adjacent faces
            tfFaces(:,end+1) = false(nFaces,1);
            for k = idxEdge
                tfFaces(:,end) = tfFaces(:,end) | any(eInfo.Faces == k,2);
            end

        end

        continue
    end

    % -> Intersection is not the edge end-point
    switch size(pnt,2)
        case 0  % pnt - []  -> No intersect
            % Do nothing
            continue
        case 1  % pnt - 3x1 -> Segment intersect plane
            
            % Append intersection point
            Xint(:,end+1);

            % Define edge index
            idxEdge = i;

            % Update tfEdges
            tfEdges(:,end+1) = false(nEdges,1);
            tfEdges(idxEdge,end) = true;

            % Define adjacent faces
            tfFaces(:,end+1) = false(nFaces,1);
            for k = idxEdge
                tfFaces(:,end) = tfFaces(:,end) | any(eInfo.Faces == k,2);
            end

        case 2  % pnt - 3x2 -> Segment on plane
            % UNEXPECTED CASE!
            % This should be addressed in:
            %   "Intersection(s) are the edge end-point(s)"
            tfEndPnt
            pnt
            error('Unexpected case');
    end

end
