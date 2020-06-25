function [adj,faceAdj] = patchFaceAdjacency(ptch)
% PATCHFACEADJACENCY returns the undirected graph adjacency identifying
% connections between faces with shared edges. 
%   adj = PATCHFACEADJACENCY(ptch) defines the patch as either a structured
%   array containing fields "Faces" and "Vertices" or a patch object.
%
%   [___,faceAdj] = PATCHFACEADJACENCY(___) additionally returns an Nx3
%   array containing face index values.
%
%   Inputs:
%       ptch - patch object or structure
%           ptch.Faces    - Nx3 array containing vertex indices
%           ptch.Vertices - Mx3 array containing vertex coordinates
%
%   Outputs:
%       adj - NxN sparse matrix containing a "1" if the corresponding face
%             indices share an edge
%       faceAdj - Nx3 array containing the face indices that are adjacent
%                 to the face defined by the row index
%           faceAdj(i,:) - the index values of the three faces adjacent to
%                          face i
%
%   NOTE: This function assumes the patch is defined using a triangular
%   mesh
%
%   M. Kutzer, 01Jun2020, USNA

%% Check input(s)
narginchk(1,1);

try
    v = ptch.Vertices;
    if size(v,2) < 3
        % Make vertices 3D
        v(:,3) = 0;
    end
    f = ptch.Faces;
    if size(f,2) > 3
        warning('This function assumes the patch is defined using a triangular mesh.');
    end
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Find faces that share vertices
N = size(f,1);
m = size(f,2); % This *should* be 3 (if the mesh is triangular)

adj = sparse(N,N);
fAdj = uint32(zeros(N,m)); % TODO - consider 8-bit, 16-bit, 32-bit switch based on number of faces
if size(f,1) > (2^32)
    error('More faces than the max number for uint32! Update this function.');
end

if N > 500
    useWaitbar = true;
else
    useWaitbar = false;
end

if useWaitbar
    wb = waitbar(0,'Finding face adjacency...','Name',mfilename);
end

for i = 1:N
    % Find shared vertices
    bin = false(N,m);
    for j = 1:m
        f_check = repmat(f(i,j),N,m);
        bin = bin | (f == f_check);
    end
    % Sum to find pairs of shared vertices
    % TODO - This step becomes more complicated for meshes that are not
    %        triangular
    sBIN = sum(bin,2);
    
    % Define edge adjacencies
    % -> Two shared vertices on a trianglular mesh - edges are shared.
    faceAdj(i,:) = reshape( uint32(find(sBIN == 2)),1,[]); % TODO - consider 8-bit, 16-bit, 32-bit switch based on number of faces
    adj(i,faceAdj(i,:)) = 1;
    
    if useWaitbar
        prc = i/N;
        waitbar(prc,wb,sprintf('Finding face adjacency (%10.6f%% complete)...',prc*100));
    end 
end
if useWaitbar
    delete(wb);
end


