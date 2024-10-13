function [Xints,iInfo] = patchPlaneIntersect(ptch,abcd,ZERO)
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
%       Xints - N-element cell array containing the vertices of the unique 
%               paths of intersection
%       iInfo - structured array containing intersection info
%
%   TODO - finalize slicing of polyshape and differentiate between regions
%          and holes
%
%   M. Kutzer, 10Oct2024, USNA

debug = false;

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

%% Patch vertices (3xM array)
Xv = ptch.Vertices.';

%% Debug plot
if debug
    fig = figure('Name','patchPlaneIntersect.m');
    axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
    view(axs,3);

    ptc = patch(ptch,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');

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

    % DEBUG
    if debug
        switch size(pnt,2)
            case 0
                plt_e(i) = plot3(axs,...
                    edgePts(1,:),edgePts(2,:),edgePts(3,:),'-k');
            case 1
                if nnz(tfEndPnt) == 0
                    plt_e(i) = plot3(axs,...
                        edgePts(1,:),edgePts(2,:),edgePts(3,:),'-c',...
                        'LineWidth',1.5);
                else
                    plt_e(i) = plot3(axs,...
                        edgePts(1,:),edgePts(2,:),edgePts(3,:),'-m',...
                        'LineWidth',1.5);
                end
            case 2
                plt_e(i) = plot3(axs,...
                    edgePts(1,:),edgePts(2,:),edgePts(3,:),'-y',...
                    'LineWidth',1.5);
        end
    end

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
            Xint(:,end+1) = pnt;

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

% DEBUG
if debug
    fprintf('tfFaces = \n');
    for i = 1:nFaces
        fprintf('\t[');
        fprintf('%d',tfFaces(i,:));
        fprintf(']\n');
    end

    fprintf('tfEdges = \n');
    for i = 1:nEdges
        fprintf('\t[');
        fprintf('%d',tfEdges(i,:));
        fprintf(']\n');
    end
end

%% Remove redundant points of intersection
% Keep only points with unique edge adjacency -----------------------------
[~,idxKeep,idxSame] = unique(tfEdges.','rows');
idxKeep = idxKeep.';
idxSame = idxSame.';

Xint = Xint(:,idxKeep.');

nKeep = numel(idxKeep);
tfEdgesTMP = false(nEdges,nKeep);
tfFacesTMP = false(nFaces,nKeep);
for i = 1:nKeep
    tfEdgesTMP(:,i) = any( tfEdges(:,idxSame == i), 2 );
    tfFacesTMP(:,i) = any( tfFaces(:,idxSame == i), 2 );

    % DEBUG
    %{
    if debug
        fprintf('tfFaces(:,%d) = any(___,2) \n',i);
        for j = 1:nFaces
            fprintf('\t[ ');
            fprintf('%d ',tfFaces(j,idxSame == i));
            fprintf('] = [ %d ]\n',tfFacesTMP(j,i));
        end

        fprintf('tfEdges(:,%d) = any(___,2) \n',i);
        for j = 1:nEdges
            fprintf('\t[ ');
            fprintf('%d ',tfEdges(j,idxSame == i));
            fprintf('] = [ %d ]\n',tfEdgesTMP(j,i));
        end
    end
    %}
end
tfEdges = tfEdgesTMP;
tfFaces = tfFacesTMP;

if debug
    fprintf('%s',repmat('v',1,60));
    fprintf('\nKeep only points with unique edge adjacency\n')
    fprintf('tfFaces = \n');
    for i = 1:nFaces
        fprintf('\t[ ');
        fprintf('%d ',tfFaces(i,:));
        fprintf(']\n');
    end

    fprintf('tfEdges = \n');
    for i = 1:nEdges
        fprintf('\t[ ');
        fprintf('%d ',tfEdges(i,:));
        fprintf(']\n');
    end
    fprintf('%s',repmat('^',1,60));
    fprintf('\n');
end
% -------------------------------------------------------------------------

% Keep only points with unique face adjacency -----------------------------
[~,idxKeep,idxSame] = unique(tfFaces.','rows');
idxKeep = idxKeep.';
idxSame = idxSame.';

Xint = Xint(:,idxKeep.');

nKeep = numel(idxKeep);
tfEdgesTMP = false(nEdges,nKeep);
tfFacesTMP = false(nFaces,nKeep);
for i = 1:nKeep
    tfEdgesTMP(:,i) = any( tfEdges(:,idxSame == i), 2 );
    tfFacesTMP(:,i) = any( tfFaces(:,idxSame == i), 2 );

    % DEBUG
    %{
    if debug
        fprintf('tfFaces(:,%d) = any(___,2) \n',i);
        for j = 1:nFaces
            fprintf('\t[ ');
            fprintf('%d ',tfFaces(j,idxSame == i));
            fprintf('] = [ %d ]\n',tfFacesTMP(j,i));
        end

        fprintf('tfEdges(:,%d) = any(___,2) \n',i);
        for j = 1:nEdges
            fprintf('\t[ ');
            fprintf('%d ',tfEdges(j,idxSame == i));
            fprintf('] = [ %d ]\n',tfEdgesTMP(j,i));
        end
    end
    %}
end
tfEdges = tfEdgesTMP;
tfFaces = tfFacesTMP;
if debug
    fprintf('%s',repmat('v',1,60));
    fprintf('\nKeep only points with unique face adjacency\n')
    fprintf('tfFaces = \n');
    for i = 1:nFaces
        fprintf('\t[ ');
        fprintf('%d ',tfFaces(i,:));
        fprintf(']\n');
    end

    fprintf('tfEdges = \n');
    for i = 1:nEdges
        fprintf('\t[ ');
        fprintf('%d ',tfEdges(i,:));
        fprintf(']\n');
    end
    fprintf('%s',repmat('^',1,60));
    fprintf('\n');
end
% -------------------------------------------------------------------------

% DEBUG
if debug
    % Plot and label intersections
    plt = plot3(axs,Xint(1,:),Xint(2,:),Xint(3,:),'.g','MarkerSize',15);
    for i = 1:size(Xint,2)
        txt(i) = text(axs,Xint(1,i),Xint(2,i),Xint(3,i),sprintf('%d',i),...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'BackgroundColor','w');
    end
end

%% Populate adjacency
N = size(Xint,2);
adjXint = false(N,N);   % NxN adjacency between intersects
for i = 1:N
    for j = (i+1):N
        % Intersections share an edge
        if any( tfEdges(:,i) & tfEdges(:,j), 1 )
            adjXint(i,j) = true;
            adjXint(j,i) = true;
        end

        % Intersections share a face
        if any( tfFaces(:,i) & tfFaces(:,j), 1 )
            adjXint(i,j) = true;
            adjXint(j,i) = true;
        end
    end
end

%% Find cycles
G = graph(adjXint);
cycles = findUndirectedGraphCycles(G);

%% Package outputs
f = ptch.Faces;
nFaces = size(f,1);
mFaces = size(f,2);
f0 = f;
f0(isnan(f)) = 0;
f0_srt = sort(f0,2);

nCycles = numel(cycles);
for i = 1:nCycles
    Xints{i} = Xint(:,cycles{i});
    
    % Check for special-case
    % -> Intersection cycle corresponds to a face
    if numel(cycles{i}) <= mFaces
        f_i = reshape(cycles{i},1,[]);
        if numel(f_i) < mFaces
            f_i(mFaces) = 0;
        end
        f_i = sort(f_i);
        f_i = repmat(f_i,nFaces,1);

        tf = all(f0_srt == f_i,2);
        if nnz(tf) == 1
            f_i = f(tf,:);
            f_i = f_i(f_i > 0);
            Xints{i} = Xint(:,f_i);
        end
    end

    % TODO - package vertex indices 

    % DEBUG
    if debug
        pltInt = plot3(axs,...
            Xints{i}(1,[1:end,1]),...
            Xints{i}(2,[1:end,1]),...
            Xints{i}(3,[1:end,1]),'-g','LineWidth',1.5);
    end
end

%% Package temporary output
iInfo.Xint = Xint;
iInfo.tfFaces = tfFaces;
iInfo.tfEdges = tfEdges;
