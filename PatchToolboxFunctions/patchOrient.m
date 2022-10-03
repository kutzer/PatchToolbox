function ptch = patchOrient(ptch,n)
% PATCHORIENT updates the face definition of a patch object to ensure that
% faces are defined in a counterclockwise sense about the outward normal 
% of each face.
%   ptch = PATCHORIENT(ptch,n) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object.
%   Normals are defined as an Nx3 array "n".
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   References:
%       [1] http://wwwf.imperial.ac.uk/~rn/centroid.pdf
%
%   M. Kutzer, 02May2019, USNA

%% Parse inputs
narginchk(2,2);

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

if size(n,2) ~= 3
    error('Normals must be specified as an Nx3 array.');
end

if size(f,1) ~= size(n,1)
    error('The number of faces must match the number of provided normals.');
end

%% Check unit normal direction
for i = 1:size(f,1)
    for j = 1:3
        p(j,:) = v(f(i,j),:);
    end
    v1 = p(2,:) - p(1,:);
    v2 = p(3,:) - p(1,:);
    
    n_now = n(i,:);
    n_chk = cross(v1,v2);
    
    if dot(n_now,n_chk) < 0
        fprintf('Flipping %d\n',i);
        f(i,:) = flip(f(i,:));
    end
end

%% Update patch
ptch.Faces = f;