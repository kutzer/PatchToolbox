function txt = patchPlotVertexInfo(ptch)
% PATCHPLOTVERTEXINFO overlays the individual vertex number
% for each vertex of a provided patch object.
%   txt = PATCHPLOTVERTEXINFO(ptch) overlays the vertex number for each
%   vertex of a patch object ptch. An array of text handles is returned.
%
%   M. Kutzer, 07May2019, USNA

% Get parent of patch object
mom = get(ptch,'Parent');

% Get vertices and faces of patch object
v = ptch.Vertices;
f = ptch.Faces;

% Create temporary axes
figTMP = figure('Visible','off');
axsTMP = axes('Parent',figTMP);
hold(axsTMP,'on');

for i = 1:size(v,1)
    v_now = v(i,:);
    
    str = sprintf('%d',i);
    
    txt(i,1) = text(v_now(1),v_now(2),v_now(3),str,'Parent',axsTMP);
end

set(txt,'Parent',mom);
delete(figTMP);