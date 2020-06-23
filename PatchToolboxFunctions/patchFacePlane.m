function abcd = patchFacePlane(ptc,idx)
% PATCHFACEPLANE calculates the coefficients for the plane(s) containing
% specified face(s) of a patch object/struct. Each plane is defined such
% that a*x + b*y + c*z + d = 0.
%   abcd = PATCHFACEPLANE(ptc) returns an Nx4 array containing the
%   coefficients of the plane of each face of the patch object/struct.
%
%   abcd = PATCHFACEPLANE(ptc,idx) allows the user to specify an Mx1 array
%   of face indices, and the function returns an Mx4 array containing the
%   unit normals to the specified faces.
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
abcd = nan(numel(idx),4);
for i = 1:numel(idx)
    % Get vertices
    v = ptc.Vertices(ptc.Faces(idx(i),:),:);
    % Calculate normal
    n = cross( (v(3,:) - v(2,:)), (v(1,:) - v(2,:)) );
    n_hat = n./norm(n);
    % Calculate d
    d = -n_hat*v(2,:).';
    % Append plane coefficients
    abcd(i,:) = [n_hat,d];
end