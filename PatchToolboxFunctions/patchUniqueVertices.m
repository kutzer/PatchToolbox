function ptchUnique = patchUniqueVertices(ptch)
% PATCHUNIQUEVERTICES reduces the size of the patch object by eliminating
% redundant vertices.
%   ptchUnique = PATCHUNIQUEVERTICES(ptch) defines the patch as either a 
%   structured array containing fields "Faces" and "Vertices" or a patch 
%   object. 
%
%   M. Kutzer, 21May2019, USNA

%% Parse/check inputs
try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Define unique vertices
[vUnique,~,idxMap] = unique(v,'Rows');
idxMap = -idxMap; % Make index values negative to keep them seperate from old index values

h = waitbar(0,'Mapping faces to unique vertices...','Name','patchUniqueVertices.m');
fUnique = f;
n = numel(idxMap);
dn = round(n/100);
for i = 1:n
    bin = fUnique == i;
    fUnique(bin) = idxMap(i);
    if mod(i,dn) == 0
        waitbar(i/n,h);
    end
end
delete(h);

%% Define output
ptchUnique.Vertices = vUnique;
ptchUnique.Faces = -fUnique;

