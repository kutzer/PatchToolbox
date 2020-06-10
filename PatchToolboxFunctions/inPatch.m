function [in,on] = inPatch(p,varargin)
% INPATCH finds the points located inside or on the edge/face of a patch
% object. 
%   in = inPatch(p,V) returns an Nx1 binary array "in" indicating which
%   vertices contained in the Nx3 array "V" is contained within the patch 
%   object/struct defined in "p."
%
%   in = inPatch(p,x,y,z) specifies the vertices using distinct arrays for
%   the x, y, and z components.
%
%   [in,on] = inPatch(___) also returns an Nx1 binary array "on" indicating
%   which vertices lie on a face or edge of the patch object/struct.
%
%   References:
%       Latombe, J.C. "Robot Motion Planning," 1991.
%
%   M. Kutzer, 10Jun2020, USNA

%% Check & parse inputs
narginchk(2,4);

if nargin == 3
    error('This function accepts 2 or 4 input arguments.');
end

if nargin == 2
    V = varargin{1};
else
    for i = 1:3
        V(:,i) = reshape(varargin{i},[],1);
    end
end

%% Define face planes
M = size(p.Faces,1);
for i = 1:M
    % Get vertices
    v = p.Vertices(p.Faces(i,:),:);
    % Calculate normal
    n = cross( (v(3,:) - v(2,:)), (v(1,:) - v(2,:)) );
    n_hat = n./norm(n);
    % Calculate d
    d = -n_hat*v(2,:).';
    % Append plane coefficients
    abcd(i,:) = [n_hat,d];
end

%% Find points in the 