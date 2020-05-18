function [patchIntersect,pathIntersect,slicedPatch] = patchEllipsoidIntersect(ptch,efit)
% PATCHELLIPSOIDINTERSECT finds the intersect between a patch and
% ellipsoid...
%   [patchIntersect,pathIntersect] = PATCHELLIPSOIDINTERSECT(ptch,efit) 
%   defines the input patch (ptch) as either a structured array containing  
%   fields "Faces" and "Vertices" or a patch object and the
%   intersecting/slicing ellipsoid using efit:
%
%       efit - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%
%   The function returns patch/ellipsoid intersect information 
%   (patchIntersect) and the path of intersection (pathIntersect)
%
%       patchIntersect - structured array containing the following fields
%           patchIntersect.Vertices - Nx3 array containing N vertices of 
%                                     the input patch object 
%           patchIntersect.Faces    - Mx3 array containing M 3-element face
%                                     index values associated with the
%                                     faces that intersect the ellipsoid
%
%       pathIntersect - k-element structured array containing closed path
%       intersect information for each path of intersection between the STL
%       and ellipsoid.
%           pathIntersect(i).Points      - Hx3 array containing points of
%                                          intersection
%           pathIntersect(i).Directions  - Hx3 array containing the
%                                          instaneous direction of the path
%                                          at the corresponding point of
%                                          intersection
%           pathIntersect(i).FaceNormals - Hx3 array containing the
%                                          normal of the patch face at the
%                                          corresponding point of
%                                          intersection
%           pathIntersect(i).SurfNormals - Hx3 array containing the
%                                          normal of the ellipsoind at the
%                                          corresponding point of
%                                          intersection
%
%   [patchIntersect,pathIntersect,slicedPatch] = PATCHELLIPSOIDINTERSECT...
%   returns interior and exterior patch information including new faces and
%   vertices created from the patch.
%           slicedPatch.Interior.Vertices
%           slicedPatch.Interior.Faces
%           slicedPatch.Exterior.Vertices
%           slicedPatch.Exterior.Faces
%
%   References (NOT CURRENTLY USED)
%       [1] P. Klein, "On the Ellipsoid and Plane Intersection Equation,"
%       Applied Mathematics, 2012, 3, 1634-1640.
%       https://file.scirp.org/pdf/AM20121100009_89014420.pdf
%
%   M. Kutzer, 07May2019, USNA

%% Debug/plot flags
plotsOn = false;

%% Parse inputs
try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Check inputs
if size(v,2) ~= 3
    error('Vertices must be specified as an Mx3 array.');
end

if size(f,2) ~= 3
    warning('This function assumes the patch is created using a triangular mesh!');
    error('Faces must be specified as an Nx3 array.');
end

%% Initialize Outputs
pathIntersect = [];
patchIntersect = [];

%% Define position/orientation of ellipsoid relative to patch frame
H_i2o = eye(4);
H_i2o(1:3,1:3) = efit.Rotation;
H_i2o(1:3,4) = efit.Center;

%% Transform vertices to the ellipsoid frame
v_tmp = transpose(v);
v_tmp(4,:) = 1;

v_i = (H_i2o^(-1)) * v_tmp;
v_i(4,:) = [];
v_i = transpose(v_i);

%% Define ellipsoid function
% frac{x^2}{a^2} + frac{y^2}{b^2} + frac{z^2}{c^2} = 1
a = efit.PrincipalRadii(1);
b = efit.PrincipalRadii(2);
c = efit.PrincipalRadii(3);
fcnEllipsoid = @(x) (x(:,1).^2)./(a^2) + (x(:,2).^2)./(b^2) + (x(:,3).^2)./(c^2);

% Represent semi-axes as a diagonal
D = diag(efit.PrincipalRadii);

%% Initialize patch output working variable
pIntersect.Vertices = v;
pIntersect.Faces = [];

%% Calculate intesects
S = [0, 1; 1, 1];
invS = S^(-1);

intFace = [];
intPoints = [];
intVertices = [];
intDirections = [];
intFaceNormals = [];
intSurfNormals = [];
intConnectivity = [];
intOnEllipsoid = [];
nFaces = size(f,1);
for idxFace = 1:nFaces
    % Get current face
    face = f(idxFace,:);
    % Get vertices of face relative to ellipsoid frame
    p(1,:) = v_i(face(1),:); % D
    p(2,:) = v_i(face(2),:); % E
    p(3,:) = v_i(face(3),:); % F
    % Define outward pointing unit normal of the face
    n = cross(p(2,:) - p(1,:), p(3,:) - p(1,:));
    n = n./norm(n);
    % Find where points lie relative to the surface of the ellipsoid
    % -> onEllipsoid < 1 - Inside the surface
    % -> onEllipsoid = 1 - On the surface
    % -> onEllipsoid > 1 - Outside the surface
    onEllipsoid = fcnEllipsoid(p);
    
    % any(onEllipsoid < 1) % -> Faces inside or intersecting the ellipsoid
    % any(onEllipsoid > 1) % -> Faces outside or intersecting the ellipsoid
    if any(onEllipsoid < 1) && any(onEllipsoid > 1) % -> Face intersects the ellipsoid
        % Append intersecting face
        pIntersect.Faces(end+1,:) = face;
        
        % Find smooth intersection
        % TODO - Investigate the correct application of [1] and/or find an
        %   alternative resource... Something is VERY wrong with the
        %   solution presented. See SCRIPT_SliceEllipsoid for work in
        %   progress.
        
        % For now, let's just do a line segment intersect approach...
        typeIntersects = '';
        nIntersects = 0;
        % Check vertex combinations for an intersection
        for ii = 1:3
            for jj = 2:3
                % Consider vertex combinations: [1,2], [1,3], [2,3]
                if ii < jj
                    % Use combination
                else
                    % Skip combination
                    continue
                end
                % Fit a parametric line to a combination of vertices
                Q = [p(ii,:); p(jj,:)].';
                M = Q * invS;
                % Find candidate intersects
                m11 = M(1,1);
                m12 = M(1,2);
                m21 = M(2,1);
                m22 = M(2,2);
                m31 = M(3,1);
                m32 = M(3,2);
                % Solve for s-values resulting in an intersection
                %   solve( ((m11*s + m12)^2)/(a^2) + ((m21*s + m22)^2)/(b^2) + ((m31*s + m32)^2)/(c^2) == 1, s)
                s(1) = -(a*b*c*(a^2*b^2*m31^2 + a^2*c^2*m21^2 - a^2*m21^2*m32^2 + 2*a^2*m21*m22*m31*m32 - a^2*m22^2*m31^2 + b^2*c^2*m11^2 - b^2*m11^2*m32^2 + 2*b^2*m11*m12*m31*m32 - b^2*m12^2*m31^2 - c^2*m11^2*m22^2 + 2*c^2*m11*m12*m21*m22 - c^2*m12^2*m21^2)^(1/2) + a^2*b^2*m31*m32 + a^2*c^2*m21*m22 + b^2*c^2*m11*m12)/(a^2*b^2*m31^2 + a^2*c^2*m21^2 + b^2*c^2*m11^2);
                s(2) = -(a^2*b^2*m31*m32 - a*b*c*(a^2*b^2*m31^2 + a^2*c^2*m21^2 - a^2*m21^2*m32^2 + 2*a^2*m21*m22*m31*m32 - a^2*m22^2*m31^2 + b^2*c^2*m11^2 - b^2*m11^2*m32^2 + 2*b^2*m11*m12*m31*m32 - b^2*m12^2*m31^2 - c^2*m11^2*m22^2 + 2*c^2*m11*m12*m21*m22 - c^2*m12^2*m21^2)^(1/2) + a^2*c^2*m21*m22 + b^2*c^2*m11*m12)/(a^2*b^2*m31^2 + a^2*c^2*m21^2 + b^2*c^2*m11^2);
                
                % Check if intersects occur on segment
                nSegIntersects = 0;
                for kk = 1:2
                    if isreal(s(kk))
                        if s(kk) >= 0 && s(kk) <= 1
                            % Intersection occurs
                            
                            % Save intersect information
                            typeIntersects = sprintf('%s[%d,%d]',typeIntersects,ii,jj);
                            nIntersects = nIntersects + 1;
                            nSegIntersects = nSegIntersects + 1;
                            
                            % Define index of face
                            intFace(end+1,:) = idxFace;
                            % Define point of intersection
                            pInt = M*[s(kk);1];
                            intPoints(end+1,:) = pInt.';
                            % Define vertex indices associated with 
                            % intersecting segment
                            [intVertices(end+1,:),idxORDER] = sort( [face(ii),face(jj)] );
                            % Keep onEllipsoid info
                            intOnEllipsoid(end+1,:) = [onEllipsoid(ii),onEllipsoid(jj)];
                            intOnEllipsoid(end,:) = intOnEllipsoid(end,idxORDER);
                            
                            % Define "direction" of intersection
                            % -> Define outward pointing normal of
                            %    ellipsoid
                            nPlane = [(2*pInt(1))/(a^2),(2*pInt(2))/(b^2),(2*pInt(3))/(c^2)];
                            nPlane = nPlane./norm(nPlane);
                            % -> Define "direction" information associated with intersection
                            direction = cross(n,nPlane);
                            direction = direction./norm(direction);
                            intDirections(end+1,:) = direction;
                            intFaceNormals(end+1,:) = n;
                            intSurfNormals(end+1,:) = nPlane;
                        end
                    end
                end
                
                % Check for special condition(s)
                if nSegIntersects > 1
                    % Multiple intersections with the same segment
                    fprintf('Face %d, Total segment intersects = %d\n',idxFace,nSegIntersects);
                    warning('This segment intersects twice! Something crazy happened!');
                end
            end
        end
        
        % Define connectivity
        if nIntersects == 2
            % Two point connectivity (easy)
            % -> Two of the three face edges each intersects the ellipsoid
            % once.
            % P1 -> P2
            intConnectivity(end+1,:) = size(intPoints,1);
            % p2 -> P1
            intConnectivity(end+1,:) = size(intPoints,1)-1;
        elseif nIntersects < 2
            % Less than two intersects detected... Something went wrong.
            % -> This should only happen if/when a vertex or single edge of
            % a face sits right on the ellipsoid.
            warning('This face intersects less than two times! Something crazy happened!');
        else
            % Something besides two-point connectivity!
            % -> This corresponds to the "Multiple intersections with the
            % same segment" warning above.
            switch typeIntersects
                case '[1,2][1,2][1,3][2,3]'
                    % P1 -> P3
                    intConnectivity(end+1,:) = size(intPoints,1)-1;
                    % P2 -> P4
                    intConnectivity(end+1,:) = size(intPoints,1);
                    % P3 -> P1
                    intConnectivity(end+1,:) = size(intPoints,1)-3;
                    % P4 -> P2
                    intConnectivity(end+1,:) = size(intPoints,1)-2;
                case '[1,2][1,3][2,3][2,3]'
                    % P1 -> P3
                    intConnectivity(end+1,:) = size(intPoints,1)-1;
                    % P2 -> P4
                    intConnectivity(end+1,:) = size(intPoints,1);
                    % P3 -> P1
                    intConnectivity(end+1,:) = size(intPoints,1)-3;
                    % P4 -> P2
                    intConnectivity(end+1,:) = size(intPoints,1)-2;
                case '[1,2][1,3][1,3][2,3]'
                    % P1 -> P2
                    intConnectivity(end+1,:) = size(intPoints,1)-2;
                    % P2 -> P1
                    intConnectivity(end+1,:) = size(intPoints,1)-3;
                    % P3 -> P4
                    intConnectivity(end+1,:) = size(intPoints,1);
                    % P4 -> P3
                    intConnectivity(end+1,:) = size(intPoints,1)-1;
                otherwise
                    % Plot specific circumstance
                    str = sprintf('Face %d',idxFace);
                    figTMP = figure('Name',str);
                    axsTMP = axes('Parent',figTMP);
                    hold(axsTMP,'on');
                    view(axsTMP,3);
                    %daspect(axsTMP,[1 1 1]);
                    % Show the face associated
                    ptcTMP = [];
                    ptcTMP.Vertices = p;
                    ptcTMP.Faces = 1:3;
                    ptcTMP = patch(ptcTMP);
                    set(ptcTMP,'FaceColor','g','EdgeColor','k');
                    % Label vertices
                    for vi = 1:3
                        text(axsTMP,p(vi,1),p(vi,2),p(vi,3),sprintf('V_{%d}',vi));
                    end
                    drawnow;
                    % Get current axes limits
                    xxlim = xlim(axsTMP);
                    yylim = ylim(axsTMP);
                    zzlim = zlim(axsTMP);
                    % Plot the ellipsoid
                    efitTMP = efit;
                    efitTMP.Rotation = eye(3);
                    efitTMP.Center = zeros(1,3);
                    plotEllipsoid(axsTMP,efitTMP);
                    % Plot the points of intersection
                    pntTMP = intPoints(end-(nIntersects-1):end,:);
                    plot3(axsTMP,pntTMP(:,1),pntTMP(:,2),pntTMP(:,3),'*k');
                    for pi = 1:size(pntTMP,1)
                        text(axsTMP,pntTMP(pi,1),pntTMP(pi,2),pntTMP(pi,3),sprintf('P_{%d}',pi));
                    end
                    % Reset axes limits to specified face
                    title(axsTMP,typeIntersects);
                    xlim(axsTMP,xxlim);
                    ylim(axsTMP,yylim);
                    zlim(axsTMP,zzlim);
                    drawnow;
                    fprintf('Face %d, Total face intersects = %d\n',idxFace,nIntersects);
                    fprintf('Intersect types: "%s"\n',typeIntersects);
                    warning('This face does not intersect twice! Something crazy happened!');
                    %pause
            end
        end
    end
end

%% Define new set of faces with intersect info
if nargout > 2
    warning('THIS METHOD IS INCOMPLETE!');
    % Define interior and exterior patches
    slicedPatch.Interior = patchEllipsoidInterior(ptch,efit);
    slicedPatch.Exterior = patchEllipsoidExterior(ptch,efit);
    % Initialize added faces
    interFaces = [];            % Faces to be appended to existing faces if/when rendering interior patch
    exterFaces = [];            % Faces to be appended to existing faces if/when rendering exterior patch
    %newVertices = intPoints;   % Vertices to be appended to existing vertices
    for idxFace = unique(intFace).'
        % Get vertices of face
        face = f(idxFace,:);
        % Isolate instances of intersections with this face
        binFace = intFace == idxFace;
        % -> Intersections
        intPointsTMP = intPoints(binFace,:);
        % -> Vertices
        intVertsTMP = intVertices(binFace,:);
        % -> Interior/exterior info
        intOnEllipsoidTMP = intOnEllipsoid(binFace,:);
        % -> Define "face index" array & Identify the vertex that is used twice
        intFaceIDX  = zeros(size(intVertsTMP));
        idxTwice = [];
        for ii = 1:numel(face)
            binTMP = intVertsTMP == face(ii);
            intFaceIDX(binTMP) = ii;
            if sum(sum( binTMP )) == 2
                idxTwice = ii;
            end
        end
        if isempty(idxTwice)
            face
            intVertsTMP
            sum(binFace)
            intFaceIDX
            warning('NO TWO!!!');
            continue;
        end
        
        % Define new vertices
        % -> NOTE: This assumes the vertices are appended onto the previous
        %          vertex list!
        newVerts = reshape( find(binFace), [], 1 ) + size(v,1);
        
        % -> Sort faces to preserve normal
        for ii = 1:2
            [intFaceIDX(ii,:),idxTMP] = sort(intFaceIDX(ii,:));
            intVertsTMP(ii,:) = intVertsTMP(ii,idxTMP);
            intOnEllipsoidTMP(ii,:) = intOnEllipsoidTMP(ii,idxTMP);
        end
        % -> Sort face pairs
        switch idxTwice
            case 3  % Vertex 3 is used twice
                [~,idxTMP] = sort(intFaceIDX(:,1));
            otherwise
                [~,idxTMP] = sort(intFaceIDX(:,2));
        end
        %intPointsTMP      = intPointsTMP(idxTMP,:);
        intVertsTMP       = intVertsTMP(idxTMP,:);
        intOnEllipsoidTMP = intOnEllipsoidTMP(idxTMP,:);
        intFaceIDX        = intFaceIDX(idxTMP,:);
        newVerts          = newVerts(idxTMP,:);
        newVerts
        % Create new faces
        switch idxTwice
            case 1
                facesNEW{1} = [face(1),newVerts(1),newVerts(2)];
                facesNEW{2} = [...
                    newVerts(1),face(2),face(3);...
                    newVerts(1),face(3),newVerts(2)];
                % Seperate interior and exterior faces
                if intOnEllipsoidTMP(1,1) < 1
                    interFaces = [interFaces; facesNEW{1}];
                    exterFaces = [exterFaces; facesNEW{2}];
                else
                    exterFaces = [exterFaces; facesNEW{1}];
                    interFaces = [interFaces; facesNEW{2}];
                end
            case 2
                facesNEW{1} = [face(2),newVerts(2),newVerts(1)];
                facesNEW{2} = [...
                    newVerts(1),face(3),face(1);...
                    newVerts(2),face(3),newVerts(1)];
                % Seperate interior and exterior faces
                if intOnEllipsoidTMP(1,2) < 1
                    interFaces = [interFaces; facesNEW{1}];
                    exterFaces = [exterFaces; facesNEW{2}];
                else
                    exterFaces = [exterFaces; facesNEW{1}];
                    interFaces = [interFaces; facesNEW{2}];
                end
            case 3
                facesNEW{1} = [face(3),newVerts(1),newVerts(2)];
                facesNEW{2} = [...
                    newVerts(1),face(1),face(2);...
                    newVerts(1),face(2),newVerts(2)];
                % Seperate interior and exterior faces
                if intOnEllipsoidTMP(1,2) < 1
                    interFaces = [interFaces; facesNEW{1}];
                    exterFaces = [exterFaces; facesNEW{2}];
                else
                    exterFaces = [exterFaces; facesNEW{1}];
                    interFaces = [interFaces; facesNEW{2}];
                end
        end
    end
    
    % Concatenate arrays
    slicedPatch.Interior.Faces = [slicedPatch.Interior.Faces; interFaces];
    slicedPatch.Exterior.Faces = [slicedPatch.Exterior.Faces; exterFaces];
    
    size(slicedPatch.Interior.Vertices)
    size(intPoints)
    slicedPatch.Interior.Vertices = [slicedPatch.Interior.Vertices; intPoints];
    slicedPatch.Exterior.Vertices = [slicedPatch.Exterior.Vertices; intPoints];
end

%% Remove redundant intersections
% Define unique set of face vertices
% -> intVertices - Unique vertex indices
% -> idx         - ordered index values of unique vertex indices
% -> idxConnect  - index values mapping unique indices back to original
%                  array
[intVertices,idx,idxConnect] = unique(intVertices,'Rows');
% Down-select intersections and intersect directions
intPoints = intPoints(idx,:);
intDirections = intDirections(idx,:);
intFaceNormals = intFaceNormals(idx,:);
intSurfNormals = intSurfNormals(idx,:);
intOnEllipsoid = intOnEllipsoid(idx,:);
% Define connectivity between unique point intersections
% -> Connectivity of calculated intersections
cellConnectivity = cell(numel(idx),1);
cellConnectivity(:,1) = {[]};
% -> Faces containing specified vertex pair
cellFace = cell(numel(idx),1);
cellFace(:,1) = {[]};
for i = 1:numel(idxConnect)
    cellConnectivity{ idxConnect(i) }(end+1) = idxConnect( intConnectivity(i) );
    cellFace{ idxConnect(i) }(end+1) = intFace(i);
end

%% Define connectivity to find "cycles" or "loops"
% USEFUL TOOLS:
% conncomp - Connected graph components
% dfsearch - Depth-first graph search

% Define adjacency
A = zeros(numel(cellConnectivity),numel(cellConnectivity));
for i = 1:numel(cellConnectivity)
    for j = 1:numel(cellConnectivity{i})
        A(i,cellConnectivity{i}(j)) = 1;
        A(cellConnectivity{i}(j),i) = 1;
    end
end

% Define graph
G = graph(A);
% Define connected components
bins = conncomp(G);%,'Type','weak');

%% Review the degree of the graph
% NOTE: Nominally all nodes should be associated with a degree of 2
degG = degree(G);
u_degG = unique(degG);
fprintf('Unique degree values:\n');
updateAdjacency = false;
for i = 1:numel(u_degG)
    % Find all instances of degree value
    bin_degG = degG == u_degG(i);
    % Define nodes associated with degree value
    node_degG = find(bin_degG);
    % Define total number of nodes associated with degree value
    n_degG = sum(bin_degG);
    % Display result
    fprintf('\tDegree: %d, Occurances: %d\n',u_degG(i),n_degG);
    
    % Display information associated with nodes no associated with degree 2
    if u_degG(i) ~= 2
        % Plot original graph
        if plotsOn
            figGRPH = figure('Name','patchEllipsoidIntersect, Graph Ordering');
            axsGRPH(1) = subplot(1,2,1,'Parent',figGRPH,'Tag','Original Graph');
            axsGRPH(2) = subplot(1,2,2,'Parent',figGRPH,'Tag','Corrected Graph');
            hold(axsGRPH(1),'on');
            hold(axsGRPH(2),'on');
            title(axsGRPH(1),'Original Graph');
            title(axsGRPH(2),'Corrected Graph');
            plot(axsGRPH(1),G);
            drawnow;
        end
        
        updateAdjacency = true;
        for j = 1:numel(node_degG)
            % Get node
            nodeINITc = node_degG(j);
            
            % Create figure to display faces/intersects
            %{
            str = sprintf('Node %d',nodeINITc);
            figTMP = figure('Name',str);
            axsTMP = axes('Parent',figTMP);
            hold(axsTMP,'on');
            view(axsTMP,3);
            pltTMP = plot3(axsTMP,intPoints(nodeINITc,1),intPoints(nodeINITc,2),intPoints(nodeINITc,3),'*m');
            %}
            
            fprintf('\tNode: %d [%d,%d], Faces: ',nodeINITc,intVertices(nodeINITc,1),intVertices(nodeINITc,2));
            for jj = 1:numel(cellFace{nodeINITc})
                idxFace = cellFace{nodeINITc}(jj);
                fprintf('%d',idxFace);
                if jj < numel(cellFace{nodeINITc})
                    fprintf(', ');
                else
                    fprintf('\n');
                end
                
                %{
                % Show the face associated
                face = f(idxFace,:);
                p(1,:) = v_i(face(1),:); % D
                p(2,:) = v_i(face(2),:); % E
                p(3,:) = v_i(face(3),:); % F
                % Patch the face
                ptcTMP = [];
                ptcTMP.Vertices = p;
                ptcTMP.Faces = 1:3;
                ptcTMP = patch(ptcTMP);
                set(ptcTMP,'FaceColor','r','EdgeColor','k');
                % Label vertices
                %{
                for vi = 1:3
                    text(axsTMP,p(vi,1),p(vi,2),p(vi,3),sprintf('V_{%d}',vi));
                end
                %}
                drawnow;
                %}
                
            end
            
            % Correct adjacency?
            distConnection = [];
            A_row = A(nodeINITc,:);
            nodesA_row = find(A_row);
            for k = 1:numel(nodesA_row)
                
                nodeGOALc = nodesA_row(k);
                distConnection(k) = norm(intPoints(nodeINITc,:) - intPoints(nodeGOALc,:));
                
                %{
                pltTMP = plot3(axsTMP,intPoints(nodeGOALc,1),intPoints(nodeGOALc,2),intPoints(nodeGOALc,3),'*k');
                pltTMP = plot3(axsTMP,...
                    [intPoints(nodeINITc,1),intPoints(nodeGOALc,1)],...
                    [intPoints(nodeINITc,2),intPoints(nodeGOALc,2)],...
                    [intPoints(nodeINITc,3),intPoints(nodeGOALc,3)],':m','LineWidth',2);
                %}
                
                fprintf('\t\t%d. Adjacency: %d (%.2f) [%d,%d], Faces: ',k,nodeGOALc,distConnection(k),...
                    intVertices(nodeGOALc,1),intVertices(nodeGOALc,2));
                
                for jj = 1:numel(cellFace{nodeGOALc})
                    idxFace = cellFace{nodeGOALc}(jj);
                    fprintf('%d',idxFace);
                    if jj < numel(cellFace{nodeGOALc})
                        fprintf(', ');
                    else
                        fprintf('\n');
                    end
                    
                    %{
                    % Show the face associated
                    face = f(idxFace,:);
                    p(1,:) = v_i(face(1),:); % D
                    p(2,:) = v_i(face(2),:); % E
                    p(3,:) = v_i(face(3),:); % F
                    % Patch the face
                    ptcTMP = [];
                    ptcTMP.Vertices = p;
                    ptcTMP.Faces = 1:3;
                    ptcTMP = patch(ptcTMP);
                    set(ptcTMP,'FaceColor','g','EdgeColor','k','FaceAlpha',0.5);
                    % Label vertices
                    %{
                    for vi = 1:3
                        text(axsTMP,p(vi,1),p(vi,2),p(vi,3),sprintf('V_{%d}',vi));
                    end
                    %}
                    drawnow
                    %}
                    
                end
                
                %{
                % Get current axes limits
                lims = axis(axsTMP);
                % Plot the ellipsoid
                efitTMP = efit;
                efitTMP.Rotation = eye(3);
                efitTMP.Center = zeros(1,3);
                plotEllipsoid(axsTMP,efitTMP);
                % Reset axes limits to specified face
                axis(axsTMP,lims);
                drawnow;
                %}
                
            end
            
            % Remove connections associated with two largest distances
            if numel(distConnection) == 4
                % Remove incorrect adjacencies
                [~,rmvIDX] = sort(distConnection);
                rmvIDX(1:2) = [];
                fprintf('\t\tRemoving:\t');
                for k = rmvIDX
                    nodeGOALc = nodesA_row(k);
                    A(nodeINITc,nodeGOALc) = 0;
                    A(nodeGOALc,nodeINITc) = 0;
                    fprintf('%d. %d ',k,nodeGOALc);
                end
                fprintf('\n');
                % Create corrected adjacency
                A(nodesA_row(rmvIDX(1)),nodesA_row(rmvIDX(2))) = 1;
                A(nodesA_row(rmvIDX(2)),nodesA_row(rmvIDX(1))) = 1;
                fprintf('\t\tConnecting:\t%d. %d %d. %d\n',rmvIDX(1),nodesA_row(rmvIDX(1)),rmvIDX(2),nodesA_row(rmvIDX(2)));
            end
            
        end
    end
end

if updateAdjacency
    % Define graph
    G = graph(A);
    % Plot the graph
    if plotsOn
        plot(axsGRPH(2),G);
    end
    % Define connected components
    bins = conncomp(G);%,'Type','weak');
end


%% Cycle through each connected component in the graph
for i = 1:max(bins)
    % Find the first node for the given connected component
    idxStart = find(bins == i,1,'first');
    % Find the ordered set using depth first search
    vORDER = dfsearch(G,idxStart);
    % TODO - there may be branches and/or cyclesloops that are not closed.
    % We need to consider some sort of check for these.
    %   NOTE: Branches and/or additional cycles should result in "zig-zags"
    %         in the path.
    % Check if cycle is closed
    % -> Create adjacency for cycle
    A_cycle = zeros(size(A));
    for idxINIT = 1:numel(vORDER)
        nodeINIT = vORDER(idxINIT);
        if idxINIT < numel(vORDER)
            idxGOAL = idxINIT+1;
        else
            idxGOAL = 1;
        end
        nodeGOAL = vORDER(idxGOAL);
        A_cycle(nodeINIT,nodeGOAL) = A(nodeINIT,nodeGOAL);
        A_cycle(nodeGOAL,nodeINIT) = A(nodeGOAL,nodeINIT);
    end
    % -> Check for closed cycle
    for idxINIT = 1:numel(vORDER)
        nodeINIT = vORDER(idxINIT);
        if idxINIT < numel(vORDER)
            idxGOAL = idxINIT+1;
        else
            idxGOAL = 1;
        end
        nodeGOAL = vORDER(idxGOAL);
        if A(nodeINIT,nodeGOAL) ~= 1
            nConnections = sum(A_cycle);
            idxBAD = find(nConnections == 1);
            unique_nConnection = unique(nConnections);
            txt = '[';
            for idxConnect = 1:numel(unique_nConnection)
                txt = sprintf('%s,%d',txt,unique_nConnection(idxConnect));
            end
            txt(end+1) = ']';
            txt2 = '[';
            for idxTMP = 1:numel(idxBAD)
                txt2 = sprintf('%s,%d',txt2,idxBAD(idxTMP));
            end
            txt2(end+1) = ']';
            warning(['Connected Component %d:\n',...
                '\tNode %d is not connected to Node %d\n',...
                '\tUnique Connections: %s\n',...
                '\tRows associated with 1 connection: %s'],...
                i,idxINIT,idxGOAL,txt,txt2);
        else
            %{
            nConnections = sum(A_cycle);
            unique_nConnection = unique(nConnections);
            txt = '[';
            for idxConnect = 1:numel(unique_nConnection)
                txt = sprintf('%s,%d',txt,unique_nConnection(idxConnect));
            end
            txt(end+1) = ']';
            fprintf(['Connected Component %d:\n',...
                '\tNode %d *is* connected to Node %d\n',...
                '\tUnique Connections: %s\n'],...
                i,idxINIT,idxGOAL,txt);
            %}
        end
    end
    
    % Define the ordered set of points and directions
    pathPoints = intPoints(vORDER,:);
    pathDirections = intDirections(vORDER,:);
    pathFaceNormals = intFaceNormals(vORDER,:);
    pathSurfNormals = intSurfNormals(vORDER,:);
    
    % Transform intersection points/directions back to the world frame
    % -> Append 1's for positions
    pathPoints = pathPoints.';
    pathPoints(4,:) = 1;
    % -> Append 0's for directions
    pathDirections = pathDirections.';
    pathDirections(4,:) = 0;
    pathFaceNormals = pathFaceNormals.';
    pathFaceNormals(4,:) = 0;
    pathSurfNormals = pathSurfNormals.';
    pathSurfNormals(4,:) = 0;
    % -> Transform to world frame
    pathPoints = H_i2o * pathPoints;
    pathDirections = H_i2o * pathDirections;
    pathFaceNormals = H_i2o * pathFaceNormals;
    pathSurfNormals = H_i2o * pathSurfNormals;
    
    % Package output
    pathIntersect(i).Points = pathPoints(1:3,:).';
    pathIntersect(i).Directions = pathDirections(1:3,:).';
    pathIntersect(i).FaceNormals = pathFaceNormals(1:3,:).';
    pathIntersect(i).SurfNormals = pathSurfNormals(1:3,:).';
end

%

%% Package output
patchIntersect = pIntersect;
