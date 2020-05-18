function SA = patchSurfaceArea(ptch)
% PATCHSURFACEAREA returns the surface area of a patch defined using a 
% triangular mesh.
%   V = PATCHSURFACEAREA(ptch) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object.
%
%   References:
%       [1] Weisstein, Eric W. "Triangle Area." From MathWorld-A Wolfram 
%       Web Resource. http://mathworld.wolfram.com/TriangleArea.html
%
%   See also patchOrient patchVolume patchCentroid patchInertia
%
%   M. Kutzer, 16May2019, USNA

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

%% Calculate Volume
nFaces = size(f,1);
SA = 0;
for i = 1:nFaces
    % Get current face
    face = f(i,:);
    % Get vertices of face
    D = v(face(1),:);
    E = v(face(2),:);
    F = v(face(3),:);
    
    % Calculate area of triangle
    % -> From Eq. (10)
    G = E-D;
    H = F-D;
    N = cross(G,H);
    A = (1/2)*norm(N);
    
    % Calculate surface area
    SA = SA + A;
end