function X = patchFaceCentroid(varargin)
% PATCHFACECENTROID calculates the centroid(s) a designated faces of a
% patch object. 
%   X = PATCHFACECENTROID(ptch,idx) calculates the centroid(s) of the faces
%   of a patch object specified by face index values specified in "idx." 
%   "ptch" defines the patch as either a structured array containing fields
%   "Faces" and "Vertices" or a patch object and returns an Nx3 array 
%   containing the centroid positions of the designated face indices.
%
%   X = PATCHFACECENTROID(ptch) calculates the centroids of all faces of
%   the patch object.
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   M. Kutzer, 02Jun2020, USNA

%% Check input(s)
narginchk(1,2);
ptch = varargin{1};

try
    v = ptch.Vertices;
    if size(v,2) < 3
        % Make vertices 3D
        v(:,3) = 0;
    end
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

switch nargin
    case 1
        idx = 1:size(f,1);
    case 2
        idx = varargin{2};
    otherwise
        error('Unexpected number of inputs.');
end

%% Calculate centroid of designated faces
X = [];
for i = idx
    % Get face vertices
    ff = f(i,:);
    v_w = v(ff,:)
    % Calculate body-fixed coordinate frame for the face vertices
    %   -> Vertices in the body-fixed frame will only have x/y coordinates
    %   with vertex 2 located at the origin.
    z = cross(v_w(3,:) - v_w(2,:), v_w(1,:) - v_w(2,:));
    y = v_w(3,:) - v_w(2,:);
    z_hat = z./norm(z);
    y_hat = y./norm(y);
    x_hat = cross(y_hat,z_hat);
    d = v_w(2,:);
    % Combine transformation
    H_o2w = eye(4);
    H_o2w(1:3,:) = [x_hat; y_hat; z_hat; d].';
    % Transform to body-fixed frame
    v_w(:,4) = 1;
    v_o = invSE(H_o2w)*v_w.';
    % Define polyshape and calculate centroid
    p = polyshape(v_o(1,:),v_o(2,:));
    [x,y] = centroid(p);
    c_o = [x; y; 0; 1];
    c_w = H_o2w*c_o;
    % Append centroid to output
    X(end+1,:) = c_w(1:3).';
end

