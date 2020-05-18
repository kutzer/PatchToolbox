function pOut = patchEllipsoidInterior(ptch,efit)
% PATCHELLIPSOIDINTERIOR finds the portion of a patch that lies within the
% interior of an ellipsoid.
%   pOut = PATCHELLIPSOIDINTERIOR(ptch,efit) defines the input 
%   patch (ptch) as either a structured array containing fields "Faces" 
%   and "Vertices" or a patch object and the slicing ellipsoid using efit:
%
%       efit - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%
%   The function returns patch information (pOut) 
%
%       pOut - structured array containing the following fields
%           pOut.Vertices  - Nx3 array containing N vertices of 
%                            the input patch object 
%           pOut.Faces     - Mx3 array containing M 3-element face
%                            index values associated with the
%                            faces that are inside of the ellipsoid
%
%   M. Kutzer, 11Sep2019, USNA

% TODO - Consider calculating/adding "sliced" faces

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

%% Initialize Outputs
pOut = [];

%% Define position/orientation of ellipsoid relative to patch frame
H_i2o = eye(4);
H_i2o(1:3,1:3) = efit.Rotation;
H_i2o(1:3,4) = efit.Center;

%% Transform vertices to the ellipsoid frame
v_tmp = transpose(v);
v_tmp(4,:) = 1;

v_i = (H_i2o^(-1)) * v_tmp;
v_i(4,:) = [];
v_i = transpose(v_i);

%% Define ellipsoid function
% frac{x^2}{a^2} + frac{y^2}{b^2} + frac{z^2}{c^2} = 1
a = efit.PrincipalRadii(1);
b = efit.PrincipalRadii(2);
c = efit.PrincipalRadii(3);
fcnEllipsoid = @(x) (x(:,1).^2)./(a^2) + (x(:,2).^2)./(b^2) + (x(:,3).^2)./(c^2);

% Represent semi-axes as a diagonal
D = diag(efit.PrincipalRadii);

%% Initialize patch output working variable
pTMP.Vertices = v;
pTMP.Faces = [];

%% Evaluate all vertices 
e_i = fcnEllipsoid(v_i);

%% Remove faces that are not within the ellipsoid
idxVerts = find( e_i >= 1 ).'; % All points outside ellipsoid
for idxVert = idxVerts
    binALL = f == idxVert;
    binROW = any(binALL,2); % Does any face contain this vertex?
    bin    = any(binROW,1); % Do any of the faces contain this vertex?
    
    if bin
        % Remove face(s)
        f(binROW,:) = [];
        
        % Exit if we ran out of faces
        if isempty(f)
            break;
        end
    end
    
end
pTMP.Faces = f;

%% Package Output(s)
pOut = pTMP;
