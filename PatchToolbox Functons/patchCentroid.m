function [X,V] = patchCentroid(ptch)
% PATCHCENTROID returns the volumetric centroid of a patch defined using 
% a triangular mesh assuming faces are defined in a counterclockwise sense 
% about the outward normal of each facet.
%   X = PATCHCENTROID(ptch) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object and
%   returns a 1x3 array containing the centroid position.
%
%   [X,V] = PATCHCENTROID(ptch) returns both the position of the centroid
%   and the volume of the patch.
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   References:
%       [1] A. Dobrovolshis, "Inertia of Any Polyhedron," ICARUS 124,
%       Article No. 0243, pp. 698-704, 1996.
%
%   See also patchOrient patchVolume
%
%   M. Kutzer, 02May2019, USNA

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

%% Calculate Center
nFaces = size(f,1);
V = 0;
R = zeros(1,3);
for i = 1:nFaces
    % Get current face
    face = f(i,:);
    % Get vertices of face
    D = v(face(1),:);
    E = v(face(2),:);
    F = v(face(3),:);
    
    % Calculate Delta R
    DeltaR = (D+E+F)./4;
    
    % Calculate DeltaV
    % Using (4) -----------------------------------------------------------
    G = E-D;
    H = F-D;
    N = cross(G,H);
    DeltaV = dot( D./3, N./2 );
    % ---------------------------------------------------------------------
    
    % Using (5) -----------------------------------------------------------
    %DeltaV = dot(D,cross(E,F))./6; % scalar triple product
    % ---------------------------------------------------------------------
    
    % Calculate R (prior to dividing by volume)
    R = R + DeltaV*DeltaR;
    % Calculate V
    V = V + DeltaV;
end

% Divide by volume to find center
R = R./V;

% Set Output
X = transpose(R);
