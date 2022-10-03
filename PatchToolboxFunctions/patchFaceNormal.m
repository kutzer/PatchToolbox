function [n_hat,n] = patchFaceNormal(ptc,idx)
% PATCHFACENORMAL calculates the normal to specified faces of a patch
% object/struct.
%   n_hat = PATCHFACENORMAL(ptc) returns an Nx3 array containing the unit
%   normals of each face of the patch object/struct.
%
%   n_hat = PATCHFACENORMAL(ptc,idx) allows the user to specify an Mx1
%   array of face indices, and the function returns an Mx3 array containing
%   the unit normals to the specified faces.
%
%   [n_hat,n] = PATCHFACENORMAL(___) also returns the normal.
%
%   Input(s)
%       ptc - patch object or structured array with fields 'Vertices' and
%            'Faces'
%           ptc.Vertices - Nx3 array containing 3D vertex coordinates
%           ptc.Faces    - Mx3 array containing vertex indices describing
%                          triangular faces
%       idx - [OPTIONAL] face index or indices for calculating the normal.
%             If idx is not specified, normals are calculated for all
%             faces.
%
%   Output(s)
%       n_hat - Kx3 array defining the unit normal for each face specified
%               in idx, or all faces if idx is not specified (i.e. K = M).
%       n     - Kx3 array defining the normal for each face
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   M. Kutzer, 12Jun2020, USNA

% Update(s)
%   03Oct2022 - Updated documentation

%% Check/parse inputs
narginchk(1,2);

try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

if nargin < 2
    M = size(ptc.Faces,1);
    idx = 1:M;
end

%% Calculate normals
n = nan(numel(idx),3);
n_hat = n;
for i = 1:numel(idx)
    % Get vertices
    v = ptc.Vertices(ptc.Faces(idx(i),:),:);
    % Calculate normal
    n(i,:) = cross( (v(3,:) - v(2,:)), (v(1,:) - v(2,:)) );
    n_hat(i,:) = n(i,:)./norm(n(i,:));
end

