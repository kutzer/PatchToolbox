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
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   M. Kutzer, 12Jun2020, USNA

%% Check/parse inputs
narginchk(1,2);

if nargin < 2
    M = size(ptc.Faces,1);
    idx = 1:M;
end

%% Calculate normals
n = nan(numel(idx),3);
n_hat = n;
for i = 1:numel(idx)
    % Get vertices
    v = p.Vertices(p.Faces(idx(i),:),:);
    % Calculate normal
    n(i,:) = cross( (v(3,:) - v(2,:)), (v(1,:) - v(2,:)) );
    n_hat(i,:) = n./norm(n);
end

