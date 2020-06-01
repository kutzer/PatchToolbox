function ptchUnique = patchUniqueFaces(ptch)
% PATCHUNIQUEFACES reduces the size of the patch object by eliminating
% redundant faces.
%   ptchUnique = PATCHUNIQUEFACES(ptch) defines the patch as either a
%   structured array containing fields "Faces" and "Vertices" or a patch
%   object.
%
%   M. Kutzer, 01Jun2020, USNA

% TODO - confirm uniform face-normal for redundant faces

%% Parse/check inputs
try
    v = ptch.Vertices;
    f = ptch.Faces;
catch
    error('Patch must be defined with "Vertices" and "Faces".');
end

%% Define unique faces
fprintf('Checking %d original faces...',size(f,1));
fSort = sort(f,2);
%[fSortUnique,idxUnique,idxMap] = unique(fSort,'Rows');
[~,idxUnique,~] = unique(fSort,'Rows');
fprintf('%d unique faces found...',numel(idxUnique));

%% Define output
ptchUnique.Vertices = v;
ptchUnique.Faces = f(idxUnique,:);
fprintf('%d faces removed.\n',size(f,1) - size(ptchUnique.Faces,1));
