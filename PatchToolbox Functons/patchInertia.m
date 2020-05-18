function I = patchInertia(ptch,rho)
% PATCHINERTIA returns the inertia tensor of a patch defined using a
% triangular mesh assuming uniform density and faces that are defined in
% a counterclockwise sense about the outward normal of each facet.
%   I = PATCHINERTIA(ptch) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object and
%   returns a 3x3 array containing the inertia tensor. Note that the
%   density (rho) is assumed to be 1 when not provided.
%
%   I = PATCHINERTIA(ptch,rho) allows the user to specify the density rho.
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh.
%
%   NOTE: This function returns the inertia tensor. When comparing to the
%   "positive" inertia matrix returned by Solidworks, the off-diagonal
%   terms will have the opposite sign.
%
%   References:
%       [1] A. Dobrovolshis, "Inertia of Any Polyhedron," ICARUS 124,
%       Article No. 0243, pp. 698-704, 1996.
%
%   M. Kutzer, 07May2019, USNA

%% Parse inputs
try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Set default(s)
if nargin < 2
    rho = 1;
end

%% Check inputs
if size(v,2) ~= 3
    error('Vertices must be specified as an Mx3 array.');
end

if size(f,2) ~= 3
    warning('This function assumes the patch is created using a triangular mesh!');
    error('Faces must be specified as an Nx3 array.');
end

%% Calculate products of intertia
% Define the total number of faces
nFaces = size(f,1);
% Initialize products of inertia
P = zeros(3,3);
% Calculate products of inertia
% 1 - x
% 2 - y
% 3 - z
for j = 1:3
    for k = 1:3
        if j > k
            % Ignore yx, zy, zx
            continue
        else
            % Calculate
            for i = 1:nFaces
                % Get current face
                face = f(i,:);
                % Get vertices of face
                D = v(face(1),:);
                E = v(face(2),:);
                F = v(face(3),:);
                
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
                
                % Calculate DeltaP
                DeltaP = ((rho*DeltaV)/20)*...
                    (2*D(j)*D(k) + 2*E(j)*E(k) + 2*F(j)*F(k) ...
                    + D(j)*E(k) + D(k)*E(j) ...
                    + D(j)*F(k) + D(k)*F(j) ...
                    + E(j)*F(k) + E(k)*F(j));
                % Calculate P_{j,k}
                P(j,k) = P(j,k) + DeltaP;
            end
        end
    end
end

%% Calculate intertial tensor
I(1,1) = P(2,2) + P(3,3);
I(2,2) = P(1,1) + P(3,3);
I(3,3) = P(1,1) + P(2,2);

I(2,3) = -P(2,3);
I(1,3) = -P(1,3);
I(1,2) = -P(1,2);

I(3,2) = -P(2,3);
I(3,1) = -P(1,3);
I(2,1) = -P(1,2);

%% Comparison Note!!!
% Solidworks returns the "positive" inertia matrix whose off-diagonal terms
% are the inverse sign of those returned using this algorithm [2].
%
% Note that the inertia tensor definition follows [1].
%
% [2] D. Saba, "Solidworks 2013 shows the 'positive' matrix of inertia," 
% https://forum.solidworks.com/thread/84904, Apr. 2015.

%{
% Convert to "positive" to match Solidworks
K = eye(3) - (ones(3,3) - eye(3));
I = K.*I;
%}
