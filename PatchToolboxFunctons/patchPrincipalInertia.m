function [H_i2o,I_i,efit] = patchPrincipalInertia(ptch,rho)
% PATCHPRINCIPALINERTIA returns the coordinate transformation defining the 
% position and orientation of the principal axes of rotation (i.e.
% principal directions) relative to the reference frame of the patch
% object, and the principal moments of inertia. These properties are
% calculated for a patch object defined using a triangular mesh assuming  
% uniform density and faces that are defined in a counterclockwise sense 
% about the outward normal of each facet.
%   H_i2o = PATCHPRINCIPALINERTIA(ptch) defines the patch as either a 
%   structured array containing fields "Faces" and "Vertices" or a patch 
%   object and returns a 4x4 rigid body transformation describing the 
%   position and orientation of the principal frame relative to the 
%   reference frame of the patch. For a single input, the density (rho) is
%   assumed to be 1 when not provided. 
%
%   H_i2o = PATCHPRINCIPALINERTIA(ptch,rho) allows the user to specify the
%   density rho.
%
%   [H_i2o,I_i] = PATCHPRINCIPALINERTIA(___) returns the 4x4 rigid body
%   transform describing the position and orientation of the principal 
%   frame and the inertia tensor relative to the principal frame. 
%
%   [H_i2o,I_i,efit] = PATCHPRINCIPALINERTIA(___) returns a structured
%   array representation of the ellipsoid with the same inertia tensor of
%   the patch. The efit returned can be used with the PLOTELLIPSOID
%   function.
%
%       efit - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh.
%
%   See also patchOrient patchCentroid patchVolume patchInertia
%
%   M. Kutzer, 17May2019, USNA

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

%% Calculate center of mass and volume
[X,V] = patchCentroid(ptch);
M = rho*V;

% Parse centroid
x = X(1);
y = X(2);
z = X(3);

%% Calculate inertia
I = patchInertia(ptch,rho);

%% Center inertia at center of mass
% Parallel axis theorem
I_c = I - M * ...
    [y^2 + z^2,      -x*y,      -x*z;...
          -x*y, x^2 + z^2,      -y*z;...
          -x*z,      -y*z, x^2 + y^2];

%% Calculate eigenvectors and eigenvalues
[R,I_i] = eig(I_c);

%% Define rigid body transformation
H_i2o = eye(4);
H_i2o(1:3,1:3) = R;
H_i2o(1:3,4) = X;

%% Define ellipsoid structured array
efit.Center = reshape(X,1,3);
efit.Rotation = R;
efit.PrincipalRadii = diag(I_i);
          