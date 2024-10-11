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
    eInfo = patchEdgeInfo(ptch);
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

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
iInfo.X = [];
iInfo.tfFaces = [];
for i = 1:size(eInfo.Edges,1)
    edgePts = ptch.Vertices(eInfo.Edges(i,:),:).';
    pnt = intersectPlaneSegment(abcd,edgePts);
    tfFaces = any(eInfo.Faces == i,2);

    if debug
        plt_e(i) = plot3(axs,edgePts(1,:),edgePts(2,:),edgePts(3,:),'-k');
    end
    
    switch size(pnt,2)
        case 1
            if ~isempty(iInfo.tfFaces)
                if any( any(iInfo.tfFaces,2) & tfFaces )
                    % Repeat face
                    % TODO - address this again
                    %continue
                end
            end

            iInfo.X(:,end+1) = pnt;
            iInfo.tfFaces(:,end+1) = tfFaces;

            if debug
                set(plt_e(i),'Color','g');
                plot3(axs,pnt(1,:),pnt(2,:),pnt(3,:),'*r');
            end

        case 2
            if ~isempty(iInfo.tfFaces)
                if any( any(iInfo.tfFaces,2) & tfFaces )
                    % Repeat face
                    % TODO - address this again
                    %continue
                end
            end

            iInfo.X(:,end+1) = pnt(:,1);
            iInfo.X(:,end+1) = pnt(:,2);
            iInfo.tfFaces(:,end+1) = tfFaces;
            iInfo.tfFaces(:,end+1) = tfFaces;

            if debug
                set(plt_e(i),'Color','m');
                plot3(axs,pnt(1,:),pnt(2,:),pnt(3,:),'*r');
            end
    end

end
