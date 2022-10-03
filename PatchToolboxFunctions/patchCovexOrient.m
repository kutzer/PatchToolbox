function [ptc,msg] = patchCovexOrient(ptc)
% PATCHCONVEXORIENT updates the face definition of a patch providing a
% triangle-based mesh of the convex polygon resulting from the convex hull.
% For convex patch objects, this will provide faces with outward facing
% normals.
%   ptc = patchCovexOrient(ptc)
%
%   Input(s)
%       ptc - patch object or structured array with fields 'Vertices' and
%            'Faces'
%           ptc.Vertices - Nx3 array containing 3D vertex coordinates
%           ptc.Faces    - Mx3 array containing vertex indices describing
%                          triangular faces
%
%   Output(s)
%       ptc - patch object or structured array with fields 'Vertices' and
%            'Faces'
%           ptc.Vertices - Nx3 array containing 3D vertex coordinates
%           ptc.Faces    - Kx3 array containing vertex indices describing
%                          triangular faces with outward facing normals.
%                          Note, if the specified patch is convex, K = M.
%       msg - character array describing unexpected changes to patch
%             object.
%
%   M. Kutzer, 03Oct2022, USNA

%% Check input(s)
narginchk(1,1);

if numel(ptc) > 1
    error('A single patch object must be specified.');
end

try
    v = ptc.Vertices;
    f = ptc.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Calculate faces of convex hull
msg = '';

try
    ptc.Faces = convhull(v(:,1),v(:,2),v(:,3));
catch ME
    msg = ME.message;
end

if size(ptc.Faces,1) ~= size(f,1)
    msg = 'Patch is not convex.';
end