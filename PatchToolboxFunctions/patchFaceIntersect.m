function [ptc1_out,ptc2_out] = patchFaceIntersect(ptc1,ptc2,idx1,idx2)
% PATCHFACEINTERSECT
%   [ptc1_out,ptc2_out] = PATCHFACEINTERSECT(ptc1,ptc2,idx1,idx2)
%
%   M. Kutzer, 12Jun2020, USNA

%% Check inputs & parse inputs
narginchk(4,4);

ptc{1} = ptc1;
ptc{2} = ptc2;
idx{1} = idx1;
idx{2} = idx2;

%% Get faces and vertices 
for i = 1:2
    if numel(idx{i}) > 1
        error('This function only works with single index values.');
    end
    
    f{i} = ptc{i}.Faces(idx{i},:);
    v{i} = ptc{i}.Vertices(f{i},:);
    
    
%% Fit planes
abcd1 = patchFacePlane(ptc1,idx1);
abcd2 = patchFacePlane(ptc2,idx2);